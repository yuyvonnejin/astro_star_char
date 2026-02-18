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
                   duration_range=(0.01, 0.5), sde_threshold=6.0,
                   n_candidates=5):
    """Detect transit signals using Box Least Squares periodogram.

    Extracts multiple candidate peaks from the BLS power spectrum so that
    the true transit signal can be identified even when short-period noise
    or systematics produce a higher peak.

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
    n_candidates : int
        Number of top BLS peaks to extract and report.

    Returns
    -------
    dict
        Transit detection results for the best candidate. Includes a
        'transit_candidates' list with all significant peaks.
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

    # Log-uniform period grid: equal density per decade of period.
    # A frequency-uniform grid packs most points at short periods (e.g., 100k
    # of 200k points below 1 day), biasing BLS toward spurious short-period
    # detections. Log-uniform gives balanced sensitivity across all periods.
    n_periods = 50000
    period_grid = np.geomspace(min_period, max_period, n_periods)

    # Duration grid (max duration must be shorter than min period)
    min_dur, max_dur = duration_range
    max_dur = min(max_dur, min_period * 0.9)
    durations = np.linspace(min_dur, max_dur, 20)

    logger.info("BLS period grid: %d points (%.2f to %.1f days), %d durations",
                len(period_grid), period_grid[0], period_grid[-1], len(durations))

    # Run BLS
    results = bls.power(period_grid, durations)

    # Compute SDE for each period (vectorized)
    power = np.asarray(results.power, dtype=float)
    mean_power = float(np.mean(power))
    std_power = float(np.std(power))

    if std_power > 0:
        sde_array = (power - mean_power) / std_power
    else:
        sde_array = np.zeros_like(power)

    # Extract candidates: one per period decade for diversity
    candidates = _extract_bls_candidates_stratified(
        period_grid, power, sde_array, results, baseline=baseline,
    )

    # Convert all candidates to result dicts with sanity checks
    all_results = []
    for cand in candidates:
        cand_result = _candidate_to_result(cand, baseline)
        all_results.append(cand_result)

    # Filter to those above SDE threshold
    above_thresh = [r for r in all_results if r["transit_sde"] >= sde_threshold]

    if not above_thresh:
        best_period = float(period_grid[np.argmax(power)])
        best_sde = float(np.max(sde_array)) if std_power > 0 else 0.0
        logger.info("SDE %.1f below threshold %.1f -- no significant transit",
                     best_sde, sde_threshold)
        return _no_transit("low_sde", period=best_period, sde=round(best_sde, 1))

    # Select best candidate: prefer "clean" (no sanity warnings) over flagged.
    # A physically plausible detection at lower SDE is more trustworthy than
    # a high-SDE detection that fails basic sanity checks.
    clean = [r for r in above_thresh if r["transit_flag"] == "ok"]
    if clean:
        best = max(clean, key=lambda r: r["transit_sde"])
        logger.info("Selected clean candidate P=%.4f d (SDE=%.1f) over %d flagged candidates",
                    best["transit_period_days"], best["transit_sde"],
                    len(above_thresh) - len(clean))
    else:
        best = max(above_thresh, key=lambda r: r["transit_sde"])

    # Rank: best first, then others by SDE descending
    others = [r for r in above_thresh if r is not best]
    others.sort(key=lambda r: r["transit_sde"], reverse=True)
    ranked = [best] + others
    for i, r in enumerate(ranked):
        r["rank"] = i + 1

    result = dict(best)
    result["transit_candidates"] = ranked

    return result


def _extract_bls_candidates_stratified(period_grid, power, sde_array,
                                        bls_results, baseline):
    """Extract the best BLS peak from each period decade using per-bin SDE.

    Divides the period range into logarithmic bins and picks the highest-
    power peak from each bin. Uses per-bin (local) SDE so that a dominant
    short-period systematic cannot suppress detection at longer periods.

    Standard (global) SDE compares every peak to the mean/std of the full
    power spectrum. When a massive noise peak at P~0.5d inflates the global
    std, real transit signals at P~20d appear insignificant even though they
    are strong relative to their local noise floor.

    Parameters
    ----------
    period_grid : ndarray
        Period grid values.
    power : ndarray
        BLS power values.
    sde_array : ndarray
        Global SDE values (kept for reference).
    bls_results : BLS results object
        Contains depth, duration, transit_time arrays.
    baseline : float
        Time baseline of the light curve (days).

    Returns
    -------
    list of dict
        Candidate peaks (one per period bin), sorted by descending local SDE.
    """
    p_min = period_grid[0]
    p_max = period_grid[-1]

    # Create ~5 log-spaced bins across the period range.
    # Aim for roughly 3 bins per decade but at least 3 bins total.
    n_decades = np.log10(p_max / p_min)
    n_bins = max(3, int(round(n_decades * 3)))
    bin_edges = np.geomspace(p_min, p_max * 1.001, n_bins + 1)

    logger.info("Stratified BLS extraction: %d bins across [%.2f, %.1f] days",
                n_bins, p_min, p_max)

    candidates = []
    for i in range(n_bins):
        mask = (period_grid >= bin_edges[i]) & (period_grid < bin_edges[i + 1])
        if not np.any(mask):
            continue

        bin_indices = np.where(mask)[0]
        bin_power = power[bin_indices]

        # Find best peak in this bin
        local_best_pos = np.argmax(bin_power)
        global_idx = bin_indices[local_best_pos]

        # Per-bin (local) SDE: peak significance relative to local noise.
        # This prevents a dominant short-period peak from suppressing
        # SDE at longer periods where the real transit may live.
        local_mean = float(np.mean(bin_power))
        local_std = float(np.std(bin_power))

        if local_std > 0:
            local_sde = float((bin_power[local_best_pos] - local_mean) / local_std)
        else:
            local_sde = 0.0

        global_sde = float(sde_array[global_idx])

        # Skip bins where best peak is very weak even locally
        if local_sde < 3.0:
            logger.info("  Bin [%.2f, %.2f] d: best P=%.4f d, local_SDE=%.1f (skip)",
                        bin_edges[i], bin_edges[i + 1],
                        float(period_grid[global_idx]), local_sde)
            continue

        logger.info("  Bin [%.2f, %.2f] d: best P=%.4f d, local_SDE=%.1f, global_SDE=%.1f",
                    bin_edges[i], bin_edges[i + 1],
                    float(period_grid[global_idx]), local_sde, global_sde)

        candidates.append({
            "period": float(period_grid[global_idx]),
            "power": float(power[global_idx]),
            "sde": local_sde,
            "sde_global": global_sde,
            "depth": float(bls_results.depth[global_idx]),
            "duration": float(bls_results.duration[global_idx]),
            "transit_time": float(bls_results.transit_time[global_idx]),
        })

    # Sort by local SDE descending
    candidates.sort(key=lambda c: c["sde"], reverse=True)
    return candidates


def _candidate_to_result(candidate, baseline):
    """Convert a candidate dict to a transit result dict with sanity checks."""
    period = candidate["period"]
    depth = candidate["depth"]
    duration = candidate["duration"]
    sde = candidate["sde"]

    n_transits = max(1, int(baseline / period))
    depth_ppm = depth * 1e6
    duration_hours = duration * 24.0

    # Sanity checks
    warnings = []
    duration_fraction = duration / period
    if duration_fraction > 0.15:
        warnings.append("duration_exceeds_15pct_of_period")
    if n_transits < 3:
        warnings.append("fewer_than_3_transits")
    if depth > 0.05:
        warnings.append("very_deep_transit")
    if depth < 0:
        warnings.append("negative_depth")

    transit_flag = "ok" if not warnings else "; ".join(warnings)

    return {
        "transit_detected": True,
        "transit_period_days": period,
        "transit_depth": depth,
        "transit_depth_ppm": round(depth_ppm, 1),
        "transit_duration_hours": round(duration_hours, 2),
        "transit_epoch": candidate["transit_time"],
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


def compute_hz_period_range(mass_Msun, luminosity_Lsun, broadening_factor=2.0):
    """Compute the orbital period range corresponding to the habitable zone.

    Converts optimistic HZ boundaries (AU) to orbital periods via Kepler's
    3rd law, then broadens by a configurable factor to account for stellar
    property uncertainties.

    Parameters
    ----------
    mass_Msun : float or None
        Stellar mass in solar masses.
    luminosity_Lsun : float or None
        Stellar luminosity in solar luminosities.
    broadening_factor : float
        Factor to broaden the period range (default 2.0).
        Inner period is divided by this factor, outer is multiplied.

    Returns
    -------
    tuple of (float, float) or None
        (min_period_days, max_period_days), or None if stellar props missing.
    """
    if mass_Msun is None or luminosity_Lsun is None:
        logger.info("HZ period range: missing stellar properties (M=%s, L=%s)",
                    mass_Msun, luminosity_Lsun)
        return None
    if mass_Msun <= 0 or luminosity_Lsun <= 0:
        logger.info("HZ period range: invalid stellar properties (M=%s, L=%s)",
                    mass_Msun, luminosity_Lsun)
        return None

    hz = compute_habitable_zone(luminosity_Lsun)
    inner_AU = hz["hz_optimistic_inner_AU"]
    outer_AU = hz["hz_optimistic_outer_AU"]

    # Kepler's 3rd law: P_years = (a_AU^3 / M_Msun)^0.5
    inner_period_yr = math.sqrt(inner_AU ** 3 / mass_Msun)
    outer_period_yr = math.sqrt(outer_AU ** 3 / mass_Msun)

    inner_period_days = inner_period_yr * 365.25
    outer_period_days = outer_period_yr * 365.25

    # Broaden range
    min_period = inner_period_days / broadening_factor
    max_period = outer_period_days * broadening_factor

    # Floor at 0.5 days
    min_period = max(0.5, min_period)

    logger.info("HZ period range: [%.1f, %.1f] days (HZ: %.2f-%.2f AU, "
                "broadened %.1fx from [%.1f, %.1f] days)",
                min_period, max_period, inner_AU, outer_AU,
                broadening_factor, inner_period_days, outer_period_days)

    return (min_period, max_period)


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


def refine_transit_depth(time, flux, best, others):
    """Re-measure transit depth after masking transits from other candidates.

    In multi-planet systems, BLS depth for each candidate is contaminated by
    transits from other planets. This function masks the other candidates'
    transits, then phase-folds at the best candidate's period to get a
    cleaner depth measurement.

    Parameters
    ----------
    time : ndarray
        Time values (days).
    flux : ndarray
        Flux values.
    best : dict
        Best candidate result dict (needs transit_period_days, transit_epoch,
        transit_duration_hours).
    others : list of dict
        Other candidate result dicts to mask out.

    Returns
    -------
    float or None
        Refined transit depth, or None if measurement failed.
    """
    flux_clean = flux.copy()
    median_flux = np.median(flux_clean)

    # Mask transits from other candidates (only clean ones -- flagged candidates
    # are likely spurious and masking them corrupts the data, especially
    # short-period false positives that can mask >30% of all data points).
    n_masked = 0
    n_skipped = 0
    for cand in others:
        # Skip flagged (spurious) candidates
        if cand.get("transit_flag", "ok") != "ok":
            n_skipped += 1
            continue

        period = cand.get("transit_period_days")
        epoch = cand.get("transit_epoch")
        duration_d = cand.get("transit_duration_hours", 0) / 24.0
        if period is None or epoch is None or duration_d <= 0:
            continue

        # Phase relative to transit center, in days
        phase = (time - epoch) % period
        # Mask 1.5x the transit duration (total width), so half-width = 0.75 * duration
        half_dur = duration_d * 0.75
        in_transit = (phase < half_dur) | (phase > period - half_dur)
        flux_clean[in_transit] = median_flux
        n_masked += int(np.sum(in_transit))

    logger.info("Depth refinement: masked %d points from %d clean candidates "
                "(%d flagged skipped)",
                n_masked, len(others) - n_skipped, n_skipped)

    # Phase-fold at best candidate's period and measure depth
    period = best["transit_period_days"]
    epoch = best["transit_epoch"]
    duration_d = best.get("transit_duration_hours", 0) / 24.0
    if duration_d <= 0:
        logger.info("Depth refinement: best candidate has no duration, skipping")
        return None

    phase = (time - epoch) % period
    half_dur = duration_d / 2.0
    in_transit = (phase < half_dur) | (phase > period - half_dur)
    out_transit = ~in_transit

    n_in = int(np.sum(in_transit))
    n_out = int(np.sum(out_transit))

    if n_in < 5 or n_out < 50:
        logger.info("Depth refinement: too few points (in=%d, out=%d)", n_in, n_out)
        return None

    median_in = float(np.median(flux_clean[in_transit]))
    median_out = float(np.median(flux_clean[out_transit]))
    depth = median_out - median_in

    logger.info("Depth refinement: in-transit median=%.6f (%d pts), "
                "out-transit median=%.6f (%d pts), depth=%.1f ppm (was %.1f ppm)",
                median_in, n_in, median_out, n_out,
                depth * 1e6, best.get("transit_depth", 0) * 1e6)

    if depth <= 0:
        logger.info("Depth refinement: refined depth <= 0, keeping original")
        return None

    return depth


def validate_even_odd(time, flux, period, epoch, duration_days):
    """Validate transit by comparing even and odd transit depths.

    Real planet transits produce equal depths at every epoch. Eclipsing
    binaries often show different primary/secondary eclipse depths,
    producing a depth ratio far from 1.0.

    Parameters
    ----------
    time : array-like
        Time values (days).
    flux : array-like
        Flux values.
    period : float
        Transit period (days).
    epoch : float
        Transit epoch (time of first transit center).
    duration_days : float
        Transit duration (days).

    Returns
    -------
    dict
        Validation results with keys: depth_even_ppm, depth_odd_ppm,
        depth_ratio_even_odd, even_odd_validation_pass, n_even, n_odd,
        even_odd_flag.
    """
    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)

    if period is None or period <= 0 or duration_days is None or duration_days <= 0:
        return {
            "depth_even_ppm": None,
            "depth_odd_ppm": None,
            "depth_ratio_even_odd": None,
            "even_odd_validation_pass": None,
            "n_even": 0,
            "n_odd": 0,
            "even_odd_flag": "invalid_parameters",
        }

    # Assign each data point an epoch number
    epoch_num = np.round((time - epoch) / period).astype(int)

    # Identify in-transit points: within half the duration of transit center
    phase = (time - epoch) - epoch_num * period
    half_dur = duration_days / 2.0
    in_transit = np.abs(phase) < half_dur

    # Compute out-of-transit baseline
    out_transit = ~in_transit
    if np.sum(out_transit) < 10:
        return {
            "depth_even_ppm": None,
            "depth_odd_ppm": None,
            "depth_ratio_even_odd": None,
            "even_odd_validation_pass": None,
            "n_even": 0,
            "n_odd": 0,
            "even_odd_flag": "insufficient_data",
        }
    baseline = float(np.median(flux[out_transit]))

    # Split in-transit points by even/odd epoch
    is_even = (epoch_num % 2) == 0
    even_mask = in_transit & is_even
    odd_mask = in_transit & ~is_even

    n_even = int(np.sum(even_mask))
    n_odd = int(np.sum(odd_mask))

    if n_even < 5 or n_odd < 5:
        return {
            "depth_even_ppm": None,
            "depth_odd_ppm": None,
            "depth_ratio_even_odd": None,
            "even_odd_validation_pass": None,
            "n_even": n_even,
            "n_odd": n_odd,
            "even_odd_flag": "too_few_points",
        }

    depth_even = baseline - float(np.median(flux[even_mask]))
    depth_odd = baseline - float(np.median(flux[odd_mask]))

    depth_even_ppm = round(depth_even * 1e6, 1)
    depth_odd_ppm = round(depth_odd * 1e6, 1)

    # Compute ratio (avoid division by zero)
    if depth_odd > 0 and depth_even > 0:
        depth_ratio = depth_even / depth_odd
    elif depth_even <= 0 and depth_odd <= 0:
        depth_ratio = 1.0  # both non-detections
    else:
        depth_ratio = float('inf')

    # Pass if ratio within [0.5, 2.0]
    validation_pass = bool(0.5 <= depth_ratio <= 2.0)
    flag = "ok" if validation_pass else "even_odd_depth_mismatch"

    if not validation_pass:
        logger.warning("Even/odd validation FAILED: depth_even=%.0f ppm, "
                       "depth_odd=%.0f ppm, ratio=%.2f (possible eclipsing binary)",
                       depth_even_ppm, depth_odd_ppm, depth_ratio)
    else:
        logger.info("Even/odd validation passed: depth_even=%.0f ppm, "
                    "depth_odd=%.0f ppm, ratio=%.2f",
                    depth_even_ppm, depth_odd_ppm, depth_ratio)

    return {
        "depth_even_ppm": depth_even_ppm,
        "depth_odd_ppm": depth_odd_ppm,
        "depth_ratio_even_odd": round(depth_ratio, 3) if depth_ratio != float('inf') else None,
        "even_odd_validation_pass": validation_pass,
        "n_even": n_even,
        "n_odd": n_odd,
        "even_odd_flag": flag,
    }


def classify_transit_shape(time, flux, period, epoch, duration_days,
                           n_phase_bins=100):
    """Classify transit shape as U-shaped (planet) or V-shaped (EB).

    Phase-folds the light curve at the transit period, bins finely within
    the transit window, and measures the flat-bottom fraction to distinguish
    box-like (U-shaped, planet) from triangular (V-shaped, grazing EB) transits.

    Parameters
    ----------
    time : array-like
        Time values (days).
    flux : array-like
        Flux values.
    period : float
        Transit period (days).
    epoch : float
        Transit epoch (days).
    duration_days : float
        Transit duration (days).
    n_phase_bins : int
        Number of phase bins within the transit window.

    Returns
    -------
    dict
        Shape classification results with keys: shape_class, flat_bottom_fraction,
        n_in_transit, shape_flag.
    """
    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)

    if (period is None or period <= 0 or duration_days is None
            or duration_days <= 0):
        return {
            "shape_class": None,
            "flat_bottom_fraction": None,
            "n_in_transit": 0,
            "shape_flag": "invalid_parameters",
        }

    # Phase-fold relative to transit center
    phase = ((time - epoch) % period) / period
    # Center on transit: shift so transit is at phase=0.5
    phase = (phase + 0.5) % 1.0

    # Select points within 1.5x transit duration of transit center (at phase=0.5)
    dur_fraction = (duration_days * 1.5) / period
    half_window = dur_fraction / 2.0
    in_window = np.abs(phase - 0.5) < half_window

    n_in_transit = int(np.sum(in_window))
    if n_in_transit < 20:
        return {
            "shape_class": None,
            "flat_bottom_fraction": None,
            "n_in_transit": n_in_transit,
            "shape_flag": "too_few_points",
        }

    window_phase = phase[in_window]
    window_flux = flux[in_window]

    # Bin the transit window
    bin_edges = np.linspace(0.5 - half_window, 0.5 + half_window, n_phase_bins + 1)
    bin_medians = []
    for i in range(n_phase_bins):
        in_bin = (window_phase >= bin_edges[i]) & (window_phase < bin_edges[i + 1])
        if np.sum(in_bin) >= 2:
            bin_medians.append(float(np.median(window_flux[in_bin])))

    if len(bin_medians) < 10:
        return {
            "shape_class": None,
            "flat_bottom_fraction": None,
            "n_in_transit": n_in_transit,
            "shape_flag": "too_few_bins",
        }

    bin_medians = np.array(bin_medians)
    min_flux = np.min(bin_medians)
    max_flux = np.max(bin_medians)
    flux_range = max_flux - min_flux

    if flux_range <= 0:
        return {
            "shape_class": "ambiguous",
            "flat_bottom_fraction": 1.0,
            "n_in_transit": n_in_transit,
            "shape_flag": "no_depth",
        }

    # Flat-bottom fraction: bins within 20% of the depth from minimum
    threshold = min_flux + 0.2 * flux_range
    n_flat = int(np.sum(bin_medians <= threshold))
    flat_fraction = n_flat / len(bin_medians)

    # Classify
    if flat_fraction > 0.3:
        shape_class = "U_shape"
    elif flat_fraction < 0.15:
        shape_class = "V_shape"
    else:
        shape_class = "ambiguous"

    flag = "ok"
    if shape_class == "V_shape":
        logger.warning("Transit shape is V-shaped (flat_fraction=%.2f) -- "
                       "possible grazing eclipsing binary", flat_fraction)
    else:
        logger.info("Transit shape: %s (flat_fraction=%.2f)", shape_class, flat_fraction)

    return {
        "shape_class": shape_class,
        "flat_bottom_fraction": round(flat_fraction, 3),
        "n_in_transit": n_in_transit,
        "shape_flag": flag,
    }


def analyze_transit(lc_data, stellar_props, variability_period=None,
                    min_period=0.5, max_period=100.0, sde_threshold=6.0,
                    hz_targeted=False, hz_broadening=2.0):
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
    hz_targeted : bool
        If True, narrow BLS period search to the habitable zone.
    hz_broadening : float
        Broadening factor for HZ period range (default 2.0).

    Returns
    -------
    dict
        Combined transit detection + planet property results.
    """
    time = lc_data["time"]
    flux = lc_data.get("flux_flat", lc_data["flux"])
    flux_err = lc_data.get("flux_err")

    # HZ-targeted mode: narrow period range to habitable zone
    if hz_targeted:
        M_star = stellar_props.get("mass_Msun")
        L_star = stellar_props.get("luminosity_Lsun")
        hz_range = compute_hz_period_range(M_star, L_star, broadening_factor=hz_broadening)
        if hz_range is not None:
            min_period = max(min_period, hz_range[0])
            max_period = min(max_period, hz_range[1])
            logger.info("HZ-targeted mode: narrowed period range to [%.1f, %.1f] days",
                        min_period, max_period)
        else:
            logger.info("HZ-targeted mode: stellar properties missing, using full period range")

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

    # Derive planet properties for best candidate
    planet_result = compute_planet_properties(transit_result, stellar_props)

    # Derive planet properties for all candidates
    candidates = transit_result.get("transit_candidates", [])
    if candidates:
        enriched_candidates = []
        for cand in candidates:
            cand_planet = compute_planet_properties(cand, stellar_props)
            merged = dict(cand)
            merged.update(cand_planet)
            enriched_candidates.append(merged)

        # Re-rank: among clean candidates, prefer HZ planets.
        # Priority: (1) clean over flagged, (2) HZ over non-HZ, (3) highest SDE.
        enriched_candidates = _rerank_candidates(enriched_candidates)
        transit_result["transit_candidates"] = enriched_candidates

        # Update top-level result if re-ranking changed the best candidate
        new_best = enriched_candidates[0]
        if new_best["transit_period_days"] != transit_result.get("transit_period_days"):
            logger.info("Re-ranked: promoting HZ candidate P=%.4f d over P=%.4f d",
                        new_best["transit_period_days"],
                        transit_result.get("transit_period_days", 0))
            # Recompute planet result for the new best
            planet_result = compute_planet_properties(new_best, stellar_props)
            # Update transit-level fields from the new best
            for key in ("transit_detected", "transit_period_days", "transit_depth",
                        "transit_depth_ppm", "transit_duration_hours", "transit_epoch",
                        "transit_sde", "n_transits_observed", "transit_flag"):
                if key in new_best:
                    transit_result[key] = new_best[key]

        # Refine best candidate's depth by masking other candidates' transits.
        # If refinement fails (depth=0, meaning the "transit" was an artifact
        # of overlapping signals), try the next candidate.
        if len(enriched_candidates) > 1:
            refined = False
            for try_idx in range(len(enriched_candidates)):
                best_cand = enriched_candidates[try_idx]
                other_cands = [c for i, c in enumerate(enriched_candidates) if i != try_idx]
                refined_depth = refine_transit_depth(time, flux, best_cand, other_cands)
                if refined_depth is not None:
                    old_depth = best_cand["transit_depth"]
                    best_cand["transit_depth_raw"] = old_depth
                    best_cand["transit_depth_raw_ppm"] = best_cand["transit_depth_ppm"]
                    best_cand["transit_depth"] = refined_depth
                    best_cand["transit_depth_ppm"] = round(refined_depth * 1e6, 1)
                    # Recompute planet properties with refined depth
                    planet_result = compute_planet_properties(best_cand, stellar_props)
                    # If this wasn't the original #1, promote it
                    if try_idx > 0:
                        logger.info("Depth refinement: candidate #%d (P=%.4f d) has no "
                                    "independent signal; promoting #%d (P=%.4f d)",
                                    enriched_candidates[0]["rank"],
                                    enriched_candidates[0]["transit_period_days"],
                                    best_cand["rank"], best_cand["transit_period_days"])
                        enriched_candidates.insert(0, enriched_candidates.pop(try_idx))
                        for i, c in enumerate(enriched_candidates):
                            c["rank"] = i + 1
                    # Update top-level transit result
                    transit_result["transit_depth"] = refined_depth
                    transit_result["transit_depth_ppm"] = round(refined_depth * 1e6, 1)
                    transit_result["transit_depth_raw_ppm"] = round(old_depth * 1e6, 1)
                    for key in ("transit_period_days", "transit_duration_hours",
                                "transit_epoch", "transit_sde", "n_transits_observed",
                                "transit_flag"):
                        if key in best_cand:
                            transit_result[key] = best_cand[key]
                    refined = True
                    break
                else:
                    logger.info("Depth refinement failed for P=%.4f d, trying next candidate",
                                best_cand["transit_period_days"])
            if not refined:
                logger.info("Depth refinement: no candidate has independent transit signal")
            transit_result["transit_candidates"] = enriched_candidates

        # Log all candidates for comparison
        logger.info("BLS found %d candidate(s):", len(enriched_candidates))
        for cand in enriched_candidates:
            a_au = cand.get("orbital_semi_major_axis_AU", "?")
            insol = cand.get("insolation_Searth", "?")
            hz = " HZ" if cand.get("in_habitable_zone") else ""
            logger.info("  #%d: P=%.4f d, SDE=%.1f, depth=%.0f ppm, a=%s AU, S=%s S_earth%s",
                        cand["rank"], cand["transit_period_days"],
                        cand["transit_sde"], cand["transit_depth_ppm"],
                        a_au, insol, hz)

    # Validation: even/odd transit depth and shape classification
    if transit_result.get("transit_detected"):
        best_period = transit_result.get("transit_period_days")
        best_epoch = transit_result.get("transit_epoch")
        best_dur_d = (transit_result.get("transit_duration_hours", 0) / 24.0)

        even_odd_result = validate_even_odd(time, flux, best_period,
                                            best_epoch, best_dur_d)
        transit_result.update(even_odd_result)

        shape_result = classify_transit_shape(time, flux, best_period,
                                              best_epoch, best_dur_d)
        transit_result.update(shape_result)

    # Merge results
    result = {}
    result.update(transit_result)
    result.update(planet_result)
    return result


def _rerank_candidates(candidates):
    """Re-rank candidates: clean+HZ first, then clean, then flagged, by SDE within each tier."""
    def sort_key(cand):
        is_clean = cand.get("transit_flag", "") == "ok" or cand.get("planet_flag", "") == "ok"
        is_hz = bool(cand.get("in_habitable_zone"))
        sde = cand.get("transit_sde", 0)
        # Sort descending: (clean, hz, sde) -- higher = better
        return (is_clean, is_hz, sde)

    ranked = sorted(candidates, key=sort_key, reverse=True)
    for i, cand in enumerate(ranked):
        cand["rank"] = i + 1
    return ranked
