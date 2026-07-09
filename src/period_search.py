import numpy as np
from astropy.timeseries import LombScargle, BoxLeastSquares


def lomb_scan(time, flux, min_period, max_period, n_freq=10000):
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


def bls_scan(time, flux, min_period, max_period, n_periods=10000):
    normalized = flux / np.median(flux)

    periods = np.linspace(min_period, max_period, n_periods)
    durations = np.linspace(0.05, 0.25, 25)

    model = BoxLeastSquares(time, normalized)
    result = model.power(periods, durations)

    best = np.argmax(result.power)

    return (
        result.period,
        result.power,
        result.period[best],
        result.duration[best],
        result.depth[best],
        result.power[best],
    )


def top_peaks(periods, power, count=5, spacing=0.05):
    order = np.argsort(power)[::-1]
    picked = []
    used = []

    for index in order:
        period = periods[index]

        if all(abs(period - old) > spacing for old in used):
            picked.append((period, power[index]))
            used.append(period)

        if len(picked) == count:
            break

    return picked


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


def fold(time, flux, period):
    phase = (time % period) / period
    order = np.argsort(phase)

    return phase[order], flux[order]
