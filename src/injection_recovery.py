import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from period_search import lomb_scan, bls_scan, top_peaks, fold


MIN_PERIOD = 1.0
MAX_PERIOD = 10.0
N_FREQ = 10000
OUTPUT_DIR = "outputs/injection"


def make_curve(days=30, points=3000, period=3.72, depth=0.018, duration=0.18, noise=0.004):
    np.random.seed(7)

    time = np.linspace(0, days, points)
    flux = np.ones_like(time)

    phase = ((time - 0.4 + 0.5 * period) % period) - 0.5 * period
    half_duration = 0.5 * duration
    in_transit = np.abs(phase) < half_duration

    flux[in_transit] -= depth

    trend = 0.002 * np.sin(2 * np.pi * time / 14)
    scatter = np.random.normal(0, noise, size=len(time))
    observed = flux + trend + scatter

    return time, observed, flux


def snr(observed, model):
    residual = observed - model
    scatter = np.std(residual)
    depth = 1 - np.min(model)

    return depth / scatter


def plot_curve(time, observed, model):
    plt.figure(figsize=(10, 5))
    plt.scatter(time, observed, s=5, alpha=0.45, label="Observed flux")
    plt.plot(time, model, linewidth=2, label="Injected transit")
    plt.xlabel("Time (days)")
    plt.ylabel("Relative flux")
    plt.title("Synthetic Transit Light Curve")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/transit_light_curve.png", dpi=300)
    plt.close()


def plot_lomb(periods, power, false_alarm):
    plt.figure(figsize=(10, 5))
    plt.plot(periods, power)
    plt.axhline(false_alarm, linestyle="--", label="1% false-alarm level")
    plt.xlabel("Period (days)")
    plt.ylabel("Lomb-Scargle power")
    plt.title("Lomb-Scargle Period Search")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/lomb_scargle_periodogram.png", dpi=300)
    plt.close()


def plot_bls(periods, power):
    plt.figure(figsize=(10, 5))
    plt.plot(periods, power)
    plt.xlabel("Period (days)")
    plt.ylabel("BLS power")
    plt.title("Box Least Squares Transit Search")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/box_least_squares_periodogram.png", dpi=300)
    plt.close()


def plot_folded(time, flux, period):
    phase, folded_flux = fold(time, flux, period)

    plt.figure(figsize=(10, 5))
    plt.scatter(phase, folded_flux, s=5, alpha=0.45)
    plt.xlabel("Orbital phase")
    plt.ylabel("Relative flux")
    plt.title("Phase-Folded Transit Light Curve")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/phase_folded_transit.png", dpi=300)
    plt.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    injected_period = 3.72

    time, observed, model = make_curve(period=injected_period)

    ls_periods, ls_power, ls_period, ls_peak, false_alarm = lomb_scan(
        time,
        observed,
        MIN_PERIOD,
        MAX_PERIOD,
        N_FREQ,
    )
    bls_periods, bls_power, bls_period, bls_duration, bls_depth, bls_peak = bls_scan(
        time,
        observed,
        MIN_PERIOD,
        MAX_PERIOD,
    )

    recovered = top_peaks(bls_periods, bls_power)
    transit_snr = snr(observed, model)

    plot_curve(time, observed, model)
    plot_lomb(ls_periods, ls_power, false_alarm)
    plot_bls(bls_periods, bls_power)
    plot_folded(time, observed, bls_period)

    table = pd.DataFrame(recovered, columns=["recovered_period_days", "bls_power"])
    table["injected_period_days"] = injected_period
    table["absolute_error_days"] = abs(table["recovered_period_days"] - injected_period)
    table.to_csv(f"{OUTPUT_DIR}/recovered_periods.csv", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "injected_period_days",
                "lomb_scargle_best_period_days",
                "box_least_squares_best_period_days",
                "box_least_squares_duration_days",
                "box_least_squares_depth",
                "transit_signal_to_noise",
                "lomb_scargle_false_alarm_level",
            ],
            "value": [
                injected_period,
                ls_period,
                bls_period,
                bls_duration,
                bls_depth,
                transit_snr,
                false_alarm,
            ],
        }
    )
    summary.to_csv(f"{OUTPUT_DIR}/detection_summary.csv", index=False)

    print("Injected period:", injected_period)
    print("Lomb-Scargle period:", round(ls_period, 4))
    print("BLS period:", round(bls_period, 4))
    print("BLS duration:", round(bls_duration, 4))
    print("BLS depth:", round(bls_depth, 4))
    print("Transit SNR:", round(transit_snr, 4))
    print("Lomb-Scargle 1% false-alarm level:", round(false_alarm, 4))

    print("\nTop BLS periods:")
    print(table)


if __name__ == "__main__":
    main()
