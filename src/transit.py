"""Module 5: Transit detection and exoplanet characterization.

5a: BLS (Box Least Squares) transit detection on flattened light curves.
5b: Planet property derivation from transit observables + stellar properties.

References:
    - Kovacs, Zucker & Mazeh (2002), A&A 391, 369 (BLS algorithm)
    - Kopparapu et al. (2013), ApJ 765, 131 (habitable zone boundaries)
    - Fulton et al. (2017), AJ 154, 109 (planet radius gap / size classes)
    - Winn (2010), Exoplanet Transits and Occultations
"""

import logging
import math

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
R_SUN_REARTH = 109.076      # R_sun in Earth radii
R_SUN_AU = 0.00465047       # R_sun in AU
R_EARTH_RJUP = 0.08921      # R_earth in Jupiter radii

# ---------------------------------------------------------------------------
# 5a: BLS Transit Detection
# ---------------------------------------------------------------------------

def detect_transit(time, flux, flux_err=None, min_period=0.5, max_period=100.0,
                   duration_range=(0.01, 0.5), sde_threshold=6.0):
    """Detect transit signals using Box Least Squares periodogram.

    Parameters
    ----------
    time : array-like
        Time values (days).
    flux : array-like
        Normalized flux values (flattened light curve).
    flux_err : array-like, optional
        Flux uncertainties. If None, equal weights assumed.
    min_period : float
        Minimum orbital period to search (days).
    max_period : float
        Maximum orbital period to search (days). Capped at half the baseline.
    duration_range : tuple of float
        (min_duration, max_duration) in days for transit duration grid.
    sde_threshold : float
        Minimum SDE for a detection to be considered significant.

    Returns
    -------
    dict
        Transit detection results. If no significant transit is found,
        transit_detected is False.
    """
    from astropy.timeseries import BoxLeastSquares

    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)

    # Validate inputs
    baseline = time[-1] - time[0]
    if baseline <= 0 or len(time) < 50:
        logger.warning("Insufficient data for BLS: baseline=%.1f d, n_points=%d",
                       baseline, len(time))
        return _no_transit("insufficient_data")

    # Cap max_period at half the baseline
    max_period = min(max_period, baseline / 2.0)
    if max_period <= min_period:
        logger.warning("max_period (%.2f d) <= min_period (%.2f d) after capping",
                       max_period, min_period)
        return _no_transit("period_range_too_narrow")

    logger.info("Running BLS: period range [%.2f, %.1f] days, %d points, baseline %.1f days",
                min_period, max_period, len(time), baseline)

    # Build astropy BLS model directly (avoids lightkurve's period grid bloat)
    if flux_err is not None:
        flux_err = np.asarray(flux_err, dtype=float)
        bls = BoxLeastSquares(time, flux, dy=flux_err)
    else:
        bls = BoxLeastSquares(time, flux)

    # Build explicit period grid: uniform in frequency (finer at short periods).
    # 200k points is sufficient for transit detection without the 5M+ default.
    n_freqs = 200000
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    freqs = np.linspace(min_freq, max_freq, n_freqs)
    period_grid = (1.0 / freqs)[::-1]  # ascending period order

    # Duration grid (max duration must be shorter than min period)
    min_dur, max_dur = duration_range
    max_dur = min(max_dur, min_period * 0.9)
    durations = np.linspace(min_dur, max_dur, 20)

    logger.info("BLS period grid: %d points (%.2f to %.1f days), %d durations",
                len(period_grid), period_grid[0], period_grid[-1], len(durations))

    # Run BLS
    results = bls.power(period_grid, durations)

    # Extract best-fit parameters
    best_idx = np.argmax(results.power)
    best_period = float(period_grid[best_idx])
    best_power = float(results.power[best_idx])
    transit_depth = float(results.depth[best_idx])
    transit_duration = float(results.duration[best_idx])

    # Get transit epoch (time of first transit) via model
    transit_time = float(results.transit_time[best_idx])

    # Compute SDE (Signal Detection Efficiency)
    power = results.power
    mean_power = np.mean(power)
    std_power = np.std(power)
    if std_power > 0:
        sde = float((np.max(power) - mean_power) / std_power)
    else:
        sde = 0.0

    logger.info("BLS result: P=%.4f d, depth=%.6f, duration=%.4f d, SDE=%.1f",
                best_period, transit_depth, transit_duration, sde)

    # Check significance
    if sde < sde_threshold:
        logger.info("SDE %.1f below threshold %.1f -- no significant transit", sde, sde_threshold)
        return _no_transit("low_sde", period=best_period, sde=sde)

    # Number of observed transits
    n_transits = max(1, int(baseline / best_period))

    # Transit depth in ppm
    depth_ppm = transit_depth * 1e6

    # Transit duration in hours
    duration_hours = transit_duration * 24.0

    # Sanity checks -- flag suspicious detections
    warnings = []
    duration_fraction = transit_duration / best_period
    if duration_fraction > 0.15:
        # Typical transits last 1-5% of the orbital period.
        # >15% is physically implausible -- likely a spurious BLS fit.
        warnings.append("duration_exceeds_15pct_of_period")
        logger.warning("Transit duration (%.2f h) is %.0f%% of period (%.4f d) "
                       "-- likely spurious detection",
                       duration_hours, duration_fraction * 100, best_period)
    if n_transits < 3:
        warnings.append("fewer_than_3_transits")
        logger.warning("Only %d transit(s) observed -- period poorly constrained",
                       n_transits)
    if transit_depth > 0.05:
        # 5% depth = Rp/R* > 0.22 -- plausible for hot Jupiters but unusual
        warnings.append("very_deep_transit")
    if transit_depth < 0:
        warnings.append("negative_depth")

    transit_flag = "ok" if not warnings else "; ".join(warnings)

    return {
        "transit_detected": True,
        "transit_period_days": best_period,
        "transit_depth": transit_depth,
        "transit_depth_ppm": round(depth_ppm, 1),
        "transit_duration_hours": round(duration_hours, 2),
        "transit_epoch": transit_time,
        "transit_sde": round(sde, 1),
        "n_transits_observed": n_transits,
        "transit_flag": transit_flag,
        "detection_method": "bls",
    }


def _no_transit(reason, period=None, sde=None):
    """Return a no-detection result dict."""
    return {
        "transit_detected": False,
        "transit_flag": reason,
        "transit_period_days": period,
        "transit_sde": sde,
        "detection_method": "bls",
    }


# ---------------------------------------------------------------------------
# 5b: Planet Property Derivation
# ---------------------------------------------------------------------------

def classify_planet_size(radius_Rearth):
    """Classify planet by radius (Fulton et al. 2017 categories).

    Parameters
    ----------
    radius_Rearth : float
        Planet radius in Earth radii.

    Returns
    -------
    str
        Size class label.
    """
    if radius_Rearth < 0.8:
        return "sub_earth"
    elif radius_Rearth < 1.25:
        return "earth_sized"
    elif radius_Rearth < 2.0:
        return "super_earth"
    elif radius_Rearth < 4.0:
        return "sub_neptune"
    elif radius_Rearth < 6.0:
        return "neptune_sized"
    elif radius_Rearth < 10.0:
        return "sub_jupiter"
    elif radius_Rearth < 15.0:
        return "jupiter_sized"
    else:
        return "super_jupiter"


def compute_habitable_zone(luminosity_Lsun):
    """Compute conservative habitable zone boundaries (Kopparapu et al. 2013).

    Uses moist greenhouse (inner) and maximum greenhouse (outer) limits.

    Parameters
    ----------
    luminosity_Lsun : float
        Stellar luminosity in solar luminosities.

    Returns
    -------
    dict
        HZ boundaries in AU: hz_conservative_inner, hz_conservative_outer,
        hz_optimistic_inner, hz_optimistic_outer.
    """
    L = luminosity_Lsun
    sqrt_L = math.sqrt(L)

    return {
        "hz_conservative_inner_AU": round(sqrt_L / math.sqrt(1.107), 4),
        "hz_conservative_outer_AU": round(sqrt_L / math.sqrt(0.356), 4),
        "hz_optimistic_inner_AU": round(sqrt_L / math.sqrt(1.776), 4),
        "hz_optimistic_outer_AU": round(sqrt_L / math.sqrt(0.320), 4),
    }


def compute_equilibrium_temp(teff_K, radius_Rsun, a_AU, albedo=0.3):
    """Compute planetary equilibrium temperature.

    T_eq = T_star * sqrt(R_star / (2 * a)) * (1 - A_bond)^(1/4)

    Parameters
    ----------
    teff_K : float
        Stellar effective temperature (K).
    radius_Rsun : float
        Stellar radius in solar radii.
    a_AU : float
        Orbital semi-major axis in AU.
    albedo : float
        Bond albedo (default 0.3, Earth-like).

    Returns
    -------
    float
        Equilibrium temperature in Kelvin.
    """
    # Convert stellar radius to AU
    R_star_AU = radius_Rsun * R_SUN_AU
    T_eq = teff_K * math.sqrt(R_star_AU / (2.0 * a_AU)) * (1.0 - albedo) ** 0.25
    return round(T_eq, 1)


def compute_planet_properties(transit_result, stellar_props):
    """Derive planet physical properties from transit + stellar parameters.

    Parameters
    ----------
    transit_result : dict
        Output from detect_transit(). Must have transit_detected=True.
    stellar_props : dict
        Stellar properties with keys: radius_Rsun, mass_Msun, teff_K,
        luminosity_Lsun. Any can be None.

    Returns
    -------
    dict
        Planet properties: radius, orbit, temperature, habitable zone status.
    """
    if not transit_result.get("transit_detected"):
        logger.info("No transit detected; skipping planet property derivation")
        return {"planet_flag": "no_transit"}

    depth = transit_result["transit_depth"]
    period_days = transit_result["transit_period_days"]

    R_star = stellar_props.get("radius_Rsun")
    M_star = stellar_props.get("mass_Msun")
    T_star = stellar_props.get("teff_K")
    L_star = stellar_props.get("luminosity_Lsun")

    result = {
        "orbital_period_days": period_days,
    }

    # Radius ratio (always computable from transit depth alone)
    rp_rstar = math.sqrt(abs(depth))
    result["planet_radius_ratio"] = round(rp_rstar, 6)

    # Check if stellar properties are available for full characterization
    missing = []
    if R_star is None:
        missing.append("radius_Rsun")
    if M_star is None:
        missing.append("mass_Msun")
    if T_star is None:
        missing.append("teff_K")
    if L_star is None:
        missing.append("luminosity_Lsun")

    if R_star is None or M_star is None:
        logger.warning("Stellar properties incomplete (%s); "
                       "planet characterization limited to radius ratio",
                       ", ".join(missing))
        result["planet_flag"] = "stellar_params_missing"
        result["planet_radius_Rearth"] = None
        result["planet_radius_Rjup"] = None
        result["planet_size_class"] = None
        result["orbital_semi_major_axis_AU"] = None
        result["equilibrium_temp_K"] = None
        result["insolation_Searth"] = None
        result["in_habitable_zone"] = None
        return result

    # Planet radius
    Rp_Rsun = R_star * rp_rstar
    Rp_Rearth = Rp_Rsun * R_SUN_REARTH
    Rp_Rjup = Rp_Rearth * R_EARTH_RJUP

    result["planet_radius_Rearth"] = round(Rp_Rearth, 2)
    result["planet_radius_Rjup"] = round(Rp_Rjup, 3)
    result["planet_size_class"] = classify_planet_size(Rp_Rearth)

    logger.info("Planet radius: %.2f R_earth (%.3f R_jup) -- %s",
                Rp_Rearth, Rp_Rjup, result["planet_size_class"])

    # Orbital semi-major axis (Kepler's 3rd law)
    # a_AU = (M_star / M_sun)^(1/3) * (P_days / 365.25)^(2/3)
    P_years = period_days / 365.25
    a_AU = (M_star ** (1.0 / 3.0)) * (P_years ** (2.0 / 3.0))
    result["orbital_semi_major_axis_AU"] = round(a_AU, 5)

    logger.info("Orbital distance: %.5f AU (P=%.4f d, M*=%.3f Msun)", a_AU, period_days, M_star)

    # Insolation flux (relative to Earth)
    if L_star is not None and a_AU > 0:
        insolation = L_star / (a_AU ** 2)
        result["insolation_Searth"] = round(insolation, 1)
    else:
        result["insolation_Searth"] = None

    # Equilibrium temperature
    if T_star is not None and R_star is not None and a_AU > 0:
        T_eq = compute_equilibrium_temp(T_star, R_star, a_AU)
        result["equilibrium_temp_K"] = T_eq
    else:
        result["equilibrium_temp_K"] = None

    # Habitable zone check
    if L_star is not None:
        hz = compute_habitable_zone(L_star)
        result.update(hz)
        inner = hz["hz_conservative_inner_AU"]
        outer = hz["hz_conservative_outer_AU"]
        result["in_habitable_zone"] = bool(inner <= a_AU <= outer)
        logger.info("HZ: [%.4f, %.4f] AU, planet at %.5f AU -> in_HZ=%s",
                     inner, outer, a_AU, result["in_habitable_zone"])
    else:
        result["in_habitable_zone"] = None

    # Planet-level sanity checks
    warnings = []
    R_star_AU = R_star * R_SUN_AU
    if a_AU < R_star_AU:
        warnings.append("orbit_inside_star")
        logger.warning("Orbital distance %.5f AU is inside the star (R*=%.5f AU)",
                       a_AU, R_star_AU)
    if Rp_Rearth > R_star * R_SUN_REARTH * 0.5:
        warnings.append("planet_larger_than_half_star")
        logger.warning("Planet radius %.2f R_earth > 50%% of star radius -- "
                       "likely eclipsing binary or false positive", Rp_Rearth)
    if result.get("equilibrium_temp_K") is not None and result["equilibrium_temp_K"] > 4000:
        warnings.append("extreme_temperature")

    # Inherit transit-level warnings
    transit_flag = transit_result.get("transit_flag", "")
    if transit_flag and transit_flag != "ok":
        warnings.append(transit_flag)

    result["planet_flag"] = "ok" if not warnings else "; ".join(warnings)
    return result


def remove_variability(time, flux, period, n_bins=50):
    """Remove stellar variability by subtracting a phase-folded template.

    Phase-folds the light curve at the given period, computes binned medians
    to build a model-free template, then subtracts it. This removes starspot
    modulation, pulsations, or other periodic variability so BLS can find
    transit dips without confusion.

    Parameters
    ----------
    time : array-like
        Time values (days).
    flux : array-like
        Flux values.
    period : float
        Variability period to remove (days).
    n_bins : int
        Number of phase bins for the template.

    Returns
    -------
    ndarray
        Residual flux with variability removed (centered on 1.0).
    """
    phase = (time % period) / period
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_medians = np.ones(n_bins)

    for i in range(n_bins):
        in_bin = (phase >= bin_edges[i]) & (phase < bin_edges[i + 1])
        if np.any(in_bin):
            bin_medians[i] = np.median(flux[in_bin])

    # Interpolate template at each data point's phase
    # Wrap-around: append first bin at end for smooth interpolation
    centers_ext = np.concatenate([bin_centers - 1, bin_centers, bin_centers + 1])
    medians_ext = np.concatenate([bin_medians, bin_medians, bin_medians])
    template = np.interp(phase, centers_ext, medians_ext)

    # Subtract template and re-center on 1.0
    residual = flux - template + np.median(flux)

    logger.info("Removed variability at P=%.4f d (%d-bin template, "
                "template amplitude=%.2f ppt)",
                period, n_bins, (np.max(bin_medians) - np.min(bin_medians)) * 1000)

    return residual


def analyze_transit(lc_data, stellar_props, variability_period=None,
                    min_period=0.5, max_period=100.0, sde_threshold=6.0):
    """Full transit analysis pipeline: BLS detection + planet characterization.

    This is the main entry point called from the pipeline.

    Parameters
    ----------
    lc_data : dict
        Light curve data with 'time', 'flux_flat' (or 'flux'), 'flux_err'.
    stellar_props : dict
        Stellar properties: radius_Rsun, mass_Msun, teff_K, luminosity_Lsun.
    variability_period : float, optional
        If set, remove this periodic variability signal before running BLS.
        Should come from Module 4 Lomb-Scargle detection.
    min_period : float
        Minimum orbital period to search (days).
    max_period : float
        Maximum orbital period to search (days).
    sde_threshold : float
        Minimum SDE for a detection to be considered significant.

    Returns
    -------
    dict
        Combined transit detection + planet property results.
    """
    time = lc_data["time"]
    flux = lc_data.get("flux_flat", lc_data["flux"])
    flux_err = lc_data.get("flux_err")

    # Pre-whiten: remove stellar variability before BLS
    variability_removed = False
    if variability_period is not None and variability_period > 0:
        logger.info("Pre-whitening: removing variability at P=%.4f d before BLS",
                    variability_period)
        flux = remove_variability(time, flux, variability_period)
        variability_removed = True

    # Run BLS
    transit_result = detect_transit(
        time, flux, flux_err,
        min_period=min_period,
        max_period=max_period,
        sde_threshold=sde_threshold,
    )

    if variability_removed:
        transit_result["variability_removed_period_days"] = variability_period

    # Derive planet properties if transit detected
    planet_result = compute_planet_properties(transit_result, stellar_props)

    # Merge results
    result = {}
    result.update(transit_result)
    result.update(planet_result)
    return result
