import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit


def find_fits_file():
    patterns = [
        "data/*.fits",
        "data/*.fit",
        "*.fits",
        "*.fit"
    ]

    files = []

    for pattern in patterns:
        files.extend(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(
            "No FITS file found. Put the Kepler light curve FITS file in data/ or the repo root."
        )

    return files[0]


def load_kepler_light_curve(path):
    with fits.open(path) as hdul:
        table = hdul[1].data
        columns = table.columns.names

        if "TIME" not in columns:
            raise ValueError("FITS table does not contain a TIME column.")

        time = np.asarray(table["TIME"], dtype=float)

        if "PDCSAP_FLUX" in columns:
            flux = np.asarray(table["PDCSAP_FLUX"], dtype=float)
            flux_source = "PDCSAP_FLUX"
        elif "SAP_FLUX" in columns:
            flux = np.asarray(table["SAP_FLUX"], dtype=float)
            flux_source = "SAP_FLUX"
        else:
            raise ValueError("FITS table does not contain PDCSAP_FLUX or SAP_FLUX.")

        if "PDCSAP_FLUX_ERR" in columns:
            flux_error = np.asarray(table["PDCSAP_FLUX_ERR"], dtype=float)
        elif "SAP_FLUX_ERR" in columns:
            flux_error = np.asarray(table["SAP_FLUX_ERR"], dtype=float)
        else:
            flux_error = np.full_like(flux, np.nan)

        if "QUALITY" in columns:
            quality = np.asarray(table["QUALITY"], dtype=int)
        else:
            quality = np.zeros_like(time, dtype=int)

    raw_count = len(time)

    finite_mask = np.isfinite(time) & np.isfinite(flux)
    quality_mask = quality == 0
    positive_mask = flux > 0

    clean_mask = finite_mask & quality_mask & positive_mask

    cleaned = pd.DataFrame({
        "time_bkjd": time[clean_mask],
        "flux": flux[clean_mask],
        "flux_error": flux_error[clean_mask],
        "quality": quality[clean_mask]
    })

    flagged = pd.DataFrame({
        "time_bkjd": time[~clean_mask],
        "flux": flux[~clean_mask],
        "flux_error": flux_error[~clean_mask],
        "quality": quality[~clean_mask]
    })

    median_flux = np.nanmedian(cleaned["flux"])
    cleaned["normalized_flux"] = cleaned["flux"] / median_flux

    return cleaned, flagged, raw_count, flux_source


def sinusoid(time, offset, amplitude, period, phase):
    return offset + amplitude * np.sin((2 * np.pi * time / period) + phase)


def run_lomb_scargle(time, normalized_flux, min_period=0.5, max_period=40):
    centered_flux = normalized_flux - np.nanmean(normalized_flux)

    frequency = np.linspace(1 / max_period, 1 / min_period, 25000)
    ls = LombScargle(time, centered_flux)
    power = ls.power(frequency)

    periods = 1 / frequency
    best_index = np.argmax(power)

    best_period = periods[best_index]
    best_power = power[best_index]

    false_alarm_level = float(np.asarray(ls.false_alarm_level(0.01)).mean())

    return periods, power, best_period, best_power, false_alarm_level


def estimate_period_uncertainty(periods, power, best_period):
    best_index = np.argmax(power)
    peak_power = power[best_index]

    baseline = np.nanmedian(power)
    half_height = baseline + 0.5 * (peak_power - baseline)

    left = best_index
    right = best_index

    while left > 0 and power[left] > half_height:
        left -= 1

    while right < len(power) - 1 and power[right] > half_height:
        right += 1

    lower_period = periods[min(left, right)]
    upper_period = periods[max(left, right)]

    uncertainty = abs(upper_period - lower_period) / 2

    if uncertainty == 0 or not np.isfinite(uncertainty):
        uncertainty = np.nan

    return uncertainty, lower_period, upper_period


def fit_sinusoid(time, normalized_flux, period_guess):
    initial = [
        1.0,
        0.5 * (np.nanmax(normalized_flux) - np.nanmin(normalized_flux)),
        period_guess,
        0.0
    ]

    bounds = (
        [0.5, -1.0, period_guess * 0.75, -2 * np.pi],
        [1.5, 1.0, period_guess * 1.25, 2 * np.pi]
    )

    params, covariance = curve_fit(
        sinusoid,
        time,
        normalized_flux,
        p0=initial,
        bounds=bounds,
        maxfev=20000
    )

    model_flux = sinusoid(time, *params)
    residuals = normalized_flux - model_flux

    fitted_period = params[2]

    if covariance is not None and np.isfinite(covariance[2, 2]):
        fitted_period_uncertainty = np.sqrt(covariance[2, 2])
    else:
        fitted_period_uncertainty = np.nan

    return params, model_flux, residuals, fitted_period, fitted_period_uncertainty


def phase_fold(time, flux, period):
    phase = (time % period) / period
    order = np.argsort(phase)

    return phase[order], flux[order]


def save_cleaned_light_curve_plot(data):
    plt.figure(figsize=(10, 5))
    plt.scatter(data["time_bkjd"], data["normalized_flux"], s=5, alpha=0.55)
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Normalized flux")
    plt.title("Cleaned Kepler Light Curve")
    plt.tight_layout()
    plt.savefig("outputs/cleaned_light_curve.png", dpi=300)
    plt.close()


def save_periodogram_plot(periods, power, best_period, false_alarm_level):
    plt.figure(figsize=(10, 5))
    plt.plot(periods, power, linewidth=1)
    plt.axvline(best_period, linestyle="--", label=f"Best period = {best_period:.4f} days")
    plt.axhline(false_alarm_level, linestyle=":", label="1% false-alarm level")
    plt.xlabel("Period (days)")
    plt.ylabel("Lomb-Scargle power")
    plt.title("Lomb-Scargle Periodogram")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/lomb_scargle_periodogram.png", dpi=300)
    plt.close()


def save_phase_folded_plot(time, flux, model_flux, period):
    phase, folded_flux = phase_fold(time, flux, period)
    model_phase, folded_model = phase_fold(time, model_flux, period)

    plt.figure(figsize=(10, 5))
    plt.scatter(phase, folded_flux, s=5, alpha=0.45, label="Observed flux")
    plt.plot(model_phase, folded_model, linewidth=2, label="Sinusoidal fit")
    plt.xlabel("Orbital phase")
    plt.ylabel("Normalized flux")
    plt.title("Phase-Folded Kepler Light Curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/phase_folded_light_curve.png", dpi=300)
    plt.close()


def save_model_plot(time, flux, model_flux):
    plt.figure(figsize=(10, 5))
    plt.scatter(time, flux, s=5, alpha=0.45, label="Observed flux")
    plt.plot(time, model_flux, linewidth=1.5, label="Sinusoidal model")
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Normalized flux")
    plt.title("Kepler Light Curve with Sinusoidal Model")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/model_fit.png", dpi=300)
    plt.close()


def save_residual_plot(time, residuals):
    plt.figure(figsize=(10, 5))
    plt.axhline(0, linestyle="--")
    plt.scatter(time, residuals, s=5, alpha=0.5)
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Residual flux")
    plt.title("Residuals from Sinusoidal Model")
    plt.tight_layout()
    plt.savefig("outputs/residuals.png", dpi=300)
    plt.close()


def main():
    os.makedirs("outputs", exist_ok=True)

    fits_path = find_fits_file()
    data, flagged, raw_count, flux_source = load_kepler_light_curve(fits_path)

    time = data["time_bkjd"].to_numpy()
    normalized_flux = data["normalized_flux"].to_numpy()

    periods, power, ls_period, ls_power, false_alarm_level = run_lomb_scargle(
        time,
        normalized_flux
    )

    ls_uncertainty, peak_lower, peak_upper = estimate_period_uncertainty(
        periods,
        power,
        ls_period
    )

    params, model_flux, residuals, fitted_period, fitted_period_uncertainty = fit_sinusoid(
        time,
        normalized_flux,
        ls_period
    )

    data["model_flux"] = model_flux
    data["residual_flux"] = residuals

    data.to_csv("outputs/cleaned_light_curve.csv", index=False)
    flagged.to_csv("outputs/flagged_removed_points.csv", index=False)

    residual_rms = np.sqrt(np.mean(residuals**2))
    residual_std = np.std(residuals, ddof=1)

    summary = pd.DataFrame({
        "metric": [
            "fits_file",
            "flux_source",
            "raw_point_count",
            "clean_point_count",
            "removed_point_count",
            "lomb_scargle_period_days",
            "lomb_scargle_period_uncertainty_days",
            "lomb_scargle_peak_power",
            "lomb_scargle_false_alarm_level",
            "sinusoid_fitted_period_days",
            "sinusoid_period_uncertainty_days",
            "period_peak_lower_bound_days",
            "period_peak_upper_bound_days",
            "residual_rms",
            "residual_std"
        ],
        "value": [
            fits_path,
            flux_source,
            raw_count,
            len(data),
            len(flagged),
            ls_period,
            ls_uncertainty,
            ls_power,
            false_alarm_level,
            fitted_period,
            fitted_period_uncertainty,
            peak_lower,
            peak_upper,
            residual_rms,
            residual_std
        ]
    })

    summary.to_csv("outputs/period_analysis_summary.csv", index=False)

    save_cleaned_light_curve_plot(data)
    save_periodogram_plot(periods, power, ls_period, false_alarm_level)
    save_model_plot(time, normalized_flux, model_flux)
    save_phase_folded_plot(time, normalized_flux, model_flux, fitted_period)
    save_residual_plot(time, residuals)

    print("Loaded FITS file:", fits_path)
    print("Flux source:", flux_source)
    print("Raw points:", raw_count)
    print("Clean points:", len(data))
    print("Removed or flagged points:", len(flagged))
    print("Lomb-Scargle period:", round(ls_period, 5), "days")
    print("Lomb-Scargle period uncertainty:", round(ls_uncertainty, 5), "days")
    print("Sinusoid fitted period:", round(fitted_period, 5), "days")
    print("Sinusoid period uncertainty:", round(fitted_period_uncertainty, 5), "days")
    print("Residual RMS:", round(residual_rms, 6))
    print("\nSaved outputs/period_analysis_summary.csv")


if __name__ == "__main__":
    main()
