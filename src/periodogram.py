"""Module 4b: Period detection and variability classification.

Uses Lomb-Scargle periodogram to find dominant periodic signals in
light curve data, computes false alarm probability, and classifies
stellar variability type.
"""

import logging
import numpy as np

logger = logging.getLogger(__name__)


def detect_period(time, flux, flux_err=None, min_period=0.02, max_period=100.0,
                  oversample_factor=5):
    """Detect the dominant period in a light curve via Lomb-Scargle.

    Parameters
    ----------
    time : array-like
        Time values (days, e.g. BJD - offset).
    flux : array-like
        Normalized flux values.
    flux_err : array-like, optional
        Flux uncertainties. If None, equal weights assumed.
    min_period : float
        Minimum period to search (days).
    max_period : float
        Maximum period to search (days). Capped at half the time baseline.
    oversample_factor : int
        Frequency grid oversampling factor. Higher = finer resolution.

    Returns
    -------
    dict
        Period detection results with keys: period_days, period_power,
        period_fap, amplitude_ppt, frequency (array), power (array),
        detection_method.
    """
    from astropy.timeseries import LombScargle

    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)

    # Cap max_period at half the time baseline (Nyquist-like limit for periods)
    baseline = time[-1] - time[0]
    if baseline <= 0:
        logger.error("Time baseline is zero or negative")
        return _empty_result()

    max_period = min(max_period, baseline / 2.0)
    if max_period <= min_period:
        logger.warning("max_period (%.2f d) <= min_period (%.2f d) after capping to half baseline",
                       max_period, min_period)
        return _empty_result()

    logger.info("Running Lomb-Scargle: period range [%.4f, %.1f] days, %d points",
                min_period, max_period, len(time))

    # Build the Lomb-Scargle model
    if flux_err is not None:
        flux_err = np.asarray(flux_err, dtype=float)
        ls = LombScargle(time, flux, flux_err)
    else:
        ls = LombScargle(time, flux)

    # Compute frequency grid
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    frequency, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        samples_per_peak=oversample_factor,
    )

    if len(power) == 0:
        logger.warning("Periodogram produced no power values")
        return _empty_result()

    # Find dominant peak
    best_idx = np.argmax(power)
    best_freq = frequency[best_idx]
    best_power = float(power[best_idx])
    best_period = float(1.0 / best_freq)

    # False alarm probability
    fap = float(ls.false_alarm_probability(best_power))

    # Estimate peak-to-peak amplitude from phase-folded data
    amplitude_ppt = _estimate_amplitude(time, flux, best_period)

    logger.info("Best period: %.4f days (power=%.4f, FAP=%.2e, amplitude=%.2f ppt)",
                best_period, best_power, fap, amplitude_ppt)

    return {
        "period_days": best_period,
        "period_power": best_power,
        "period_fap": fap,
        "amplitude_ppt": amplitude_ppt,
        "detection_method": "lomb_scargle",
        "frequency": frequency,
        "power": power,
    }


def classify_variability(period_result, n_points=None, time_baseline_days=None):
    """Classify stellar variability based on period detection results.

    Simple rule-based classification for v2.0.

    Parameters
    ----------
    period_result : dict
        Output from detect_period().
    n_points : int, optional
        Number of data points (for quality flags).
    time_baseline_days : float, optional
        Total observation span (for quality flags).

    Returns
    -------
    dict
        Classification with keys: variability_class, variability_flag,
        period_days, period_fap, amplitude_ppt.
    """
    fap = period_result.get("period_fap")
    period = period_result.get("period_days")
    amplitude = period_result.get("amplitude_ppt", 0)
    power = period_result.get("period_power", 0)

    # Quality flags
    flag = "ok"
    if n_points is not None and n_points < 100:
        flag = "few_points"
    if time_baseline_days is not None and time_baseline_days < 1.0:
        flag = "short_baseline"

    # No detection
    if period is None or fap is None:
        return {
            "variability_class": "undetermined",
            "variability_flag": flag if flag != "ok" else "no_data",
            "period_days": None,
            "period_fap": None,
            "amplitude_ppt": None,
        }

    # Classification rules
    if fap < 0.001:
        var_class = "periodic"
    elif fap < 0.01:
        var_class = "possible_periodic"
    else:
        var_class = "non_variable"

    return {
        "variability_class": var_class,
        "variability_flag": flag,
        "period_days": period,
        "period_fap": fap,
        "amplitude_ppt": amplitude,
    }


def analyze_lightcurve(lc_data, min_period=0.02, max_period=100.0):
    """Full period analysis pipeline: detect period + classify.

    Parameters
    ----------
    lc_data : dict
        Output from lightcurve.retrieve_lightcurve() or similar.
        Must contain 'flux_flat' (or 'flux') and 'time' arrays.
    min_period : float
        Minimum period to search (days).
    max_period : float
        Maximum period to search (days).

    Returns
    -------
    dict
        Combined results: period detection + variability classification +
        light curve metadata. Arrays (frequency, power) are excluded from
        the summary dict for clean JSON serialization.
    """
    time = lc_data["time"]
    flux = lc_data.get("flux_flat", lc_data["flux"])
    flux_err = lc_data.get("flux_err")

    # Detect period
    period_result = detect_period(time, flux, flux_err,
                                  min_period=min_period, max_period=max_period)

    # Classify
    classification = classify_variability(
        period_result,
        n_points=lc_data.get("n_points_clean", len(time)),
        time_baseline_days=lc_data.get("time_baseline_days"),
    )

    # Build output dict (exclude large arrays for JSON output)
    result = {
        "lightcurve_available": True,
        "lc_mission": lc_data.get("mission"),
        "lc_author": lc_data.get("author"),
        "lc_n_sectors": lc_data.get("n_sectors"),
        "lc_n_points": lc_data.get("n_points_clean", len(time)),
        "lc_time_baseline_days": lc_data.get("time_baseline_days"),
        "lc_cadence_s": lc_data.get("cadence_s"),
    }
    result.update(classification)

    return result


def _estimate_amplitude(time, flux, period):
    """Estimate peak-to-peak amplitude in parts-per-thousand from phase-folded data.

    Bins the phase-folded light curve and takes the difference between
    the brightest and faintest bins to reduce noise influence.
    """
    phase = (time % period) / period
    n_bins = 20
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_medians = []
    for i in range(n_bins):
        in_bin = (phase >= bin_edges[i]) & (phase < bin_edges[i + 1])
        if np.any(in_bin):
            bin_medians.append(np.median(flux[in_bin]))

    if len(bin_medians) < 2:
        return 0.0

    amplitude = (max(bin_medians) - min(bin_medians)) * 1000.0  # parts per thousand
    return round(float(amplitude), 2)


def _empty_result():
    """Return empty period detection result."""
    return {
        "period_days": None,
        "period_power": None,
        "period_fap": None,
        "amplitude_ppt": None,
        "detection_method": "lomb_scargle",
        "frequency": np.array([]),
        "power": np.array([]),
    }
