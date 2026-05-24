import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit


MIN_PERIOD = 0.5
MAX_PERIOD = 40.0
N_FREQ = 25000


def find_fits(path=None):
    if path:
        if not os.path.exists(path):
            raise FileNotFoundError(f"FITS file not found: {path}")
        return path

    patterns = [
        "data/*.fits",
        "data/*.fit",
        "*.fits",
        "*.fit",
    ]

    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    files = sorted(files)

    if not files:
        raise FileNotFoundError(
            "No FITS file found. Put the Kepler light curve FITS file in data/ or the repo root."
        )

    return files[0]


def load_curve(path):
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

    finite = np.isfinite(time) & np.isfinite(flux)
    good_quality = quality == 0
    positive_flux = flux > 0
    keep = finite & good_quality & positive_flux

    data = pd.DataFrame(
        {
            "time_bkjd": time[keep],
            "flux": flux[keep],
            "flux_error": flux_error[keep],
            "quality": quality[keep],
        }
    )

    flagged = pd.DataFrame(
        {
            "time_bkjd": time[~keep],
            "flux": flux[~keep],
            "flux_error": flux_error[~keep],
            "quality": quality[~keep],
        }
    )

    median_flux = np.nanmedian(data["flux"])
    data["normalized_flux"] = data["flux"] / median_flux

    return data, flagged, raw_count, flux_source


def sine_model(time, offset, amplitude, period, phase):
    return offset + amplitude * np.sin((2 * np.pi * time / period) + phase)


def lomb_scan(time, flux, min_period=MIN_PERIOD, max_period=MAX_PERIOD, n_freq=N_FREQ):
    centered_flux = flux - np.nanmean(flux)

    frequency = np.linspace(1 / max_period, 1 / min_period, n_freq)
    ls = LombScargle(time, centered_flux)
    power = ls.power(frequency)

    periods = 1 / frequency
    best_index = np.argmax(power)

    best_period = periods[best_index]
    best_power = power[best_index]
    false_alarm_level = float(np.asarray(ls.false_alarm_level(0.01)).mean())

    return periods, power, best_period, best_power, false_alarm_level


def peak_width(periods, power):
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

    edge1 = periods[left]
    edge2 = periods[right]

    lower_period = min(edge1, edge2)
    upper_period = max(edge1, edge2)
    uncertainty = 0.5 * (upper_period - lower_period)

    if uncertainty == 0 or not np.isfinite(uncertainty):
        uncertainty = np.nan

    return uncertainty, lower_period, upper_period


def fit_wave(time, flux, period_guess):
    initial = [
        1.0,
        0.5 * (np.nanmax(flux) - np.nanmin(flux)),
        period_guess,
        0.0,
    ]

    bounds = (
        [0.5, -1.0, period_guess * 0.75, -2 * np.pi],
        [1.5, 1.0, period_guess * 1.25, 2 * np.pi],
    )

    params, covariance = curve_fit(
        sine_model,
        time,
        flux,
        p0=initial,
        bounds=bounds,
        maxfev=20000,
    )

    model_flux = sine_model(time, *params)
    residuals = flux - model_flux
    fitted_period = params[2]

    if covariance is not None and np.isfinite(covariance[2, 2]):
        formal_period_uncertainty = np.sqrt(covariance[2, 2])
    else:
        formal_period_uncertainty = np.nan

    return params, model_flux, residuals, fitted_period, formal_period_uncertainty


def fold(time, flux, period):
    phase = (time % period) / period
    order = np.argsort(phase)

    return phase[order], flux[order]


def plot_curve(data):
    plt.figure(figsize=(10, 5))
    plt.scatter(data["time_bkjd"], data["normalized_flux"], s=5, alpha=0.55)
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Normalized flux")
    plt.title("Cleaned light curve")
    plt.tight_layout()
    plt.savefig("outputs/cleaned_light_curve.png", dpi=300)
    plt.close()


def plot_periodogram(periods, power, best_period, false_alarm_level):
    plt.figure(figsize=(10, 5))
    plt.plot(periods, power, linewidth=1)
    plt.axvline(best_period, linestyle="--", label=f"Best period = {best_period:.4f} days")
    plt.axhline(false_alarm_level, linestyle=":", label="1% false-alarm level")
    plt.xlabel("Period (days)")
    plt.ylabel("Lomb-Scargle power")
    plt.title("Period search")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/lomb_scargle_periodogram.png", dpi=300)
    plt.close()


def plot_folded(time, flux, model_flux, period):
    phase, folded_flux = fold(time, flux, period)
    model_phase, folded_model = fold(time, model_flux, period)

    plt.figure(figsize=(10, 5))
    plt.scatter(phase, folded_flux, s=5, alpha=0.45, label="Observed flux")
    plt.plot(model_phase, folded_model, linewidth=2, label="Sine fit")
    plt.xlabel("Folded phase")
    plt.ylabel("Normalized flux")
    plt.title("Folded light curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/phase_folded_light_curve.png", dpi=300)
    plt.close()


def plot_fit(time, flux, model_flux):
    plt.figure(figsize=(10, 5))
    plt.scatter(time, flux, s=5, alpha=0.45, label="Observed flux")
    plt.plot(time, model_flux, linewidth=1.5, label="Sine model")
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Normalized flux")
    plt.title("Model fit")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/model_fit.png", dpi=300)
    plt.close()


def plot_residuals(time, residuals):
    plt.figure(figsize=(10, 5))
    plt.axhline(0, linestyle="--")
    plt.scatter(time, residuals, s=5, alpha=0.5)
    plt.xlabel("Time (BKJD)")
    plt.ylabel("Residual flux")
    plt.title("Fit residuals")
    plt.tight_layout()
    plt.savefig("outputs/residuals.png", dpi=300)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze a Kepler light curve FITS file."
    )
    parser.add_argument(
        "fits_file",
        nargs="?",
        help="Optional path to a Kepler FITS file. If omitted, the script searches data/ and the repo root.",
    )
    return parser.parse_args()


def main():
    os.makedirs("outputs", exist_ok=True)

    args = parse_args()
    fits_path = find_fits(args.fits_file)

    data, flagged, raw_count, flux_source = load_curve(fits_path)

    time = data["time_bkjd"].to_numpy()
    flux = data["normalized_flux"].to_numpy()

    periods, power, ls_period, ls_power, false_alarm_level = lomb_scan(time, flux)
    ls_uncertainty, peak_lower, peak_upper = peak_width(periods, power)

    params, model_flux, residuals, fitted_period, formal_fit_uncertainty = fit_wave(
        time,
        flux,
        ls_period,
    )

    data["model_flux"] = model_flux
    data["residual_flux"] = residuals

    data.to_csv("outputs/cleaned_light_curve.csv", index=False)
    flagged.to_csv("outputs/flagged_removed_points.csv", index=False)

    residual_rms = np.sqrt(np.mean(residuals**2))
    residual_std = np.std(residuals, ddof=1)

    summary = pd.DataFrame(
        {
            "metric": [
                "fits_file",
                "flux_source",
                "period_search_min_days",
                "period_search_max_days",
                "period_search_frequency_samples",
                "raw_point_count",
                "clean_point_count",
                "removed_point_count",
                "lomb_scargle_period_days",
                "lomb_scargle_peak_width_uncertainty_days",
                "lomb_scargle_peak_power",
                "lomb_scargle_false_alarm_level",
                "sine_fit_period_days",
                "formal_fit_period_uncertainty_days",
                "period_peak_lower_bound_days",
                "period_peak_upper_bound_days",
                "residual_rms",
                "residual_std",
            ],
            "value": [
                fits_path,
                flux_source,
                MIN_PERIOD,
                MAX_PERIOD,
                N_FREQ,
                raw_count,
                len(data),
                len(flagged),
                ls_period,
                ls_uncertainty,
                ls_power,
                false_alarm_level,
                fitted_period,
                formal_fit_uncertainty,
                peak_lower,
                peak_upper,
                residual_rms,
                residual_std,
            ],
        }
    )

    summary.to_csv("outputs/period_analysis_summary.csv", index=False)

    plot_curve(data)
    plot_periodogram(periods, power, ls_period, false_alarm_level)
    plot_fit(time, flux, model_flux)
    plot_folded(time, flux, model_flux, fitted_period)
    plot_residuals(time, residuals)

    print("Loaded FITS file:", fits_path)
    print("Flux source:", flux_source)
    print("Raw points:", raw_count)
    print("Clean points:", len(data))
    print("Removed or flagged points:", len(flagged))
    print("Period search window:", f"{MIN_PERIOD} to {MAX_PERIOD} days")
    print("Lomb-Scargle period:", round(ls_period, 5), "days")
    print("Lomb-Scargle peak-width uncertainty:", round(ls_uncertainty, 5), "days")
    print("Sine-fit period:", round(fitted_period, 5), "days")
    print("Formal fit period uncertainty:", round(formal_fit_uncertainty, 5), "days")
    print("Residual RMS:", round(residual_rms, 6))
    print("\nSaved outputs/period_analysis_summary.csv")


if __name__ == "__main__":
    main()
