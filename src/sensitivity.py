"""Phase 7b: Combined detection sensitivity maps.

Computes minimum detectable planet parameters from transit photometry,
radial velocities, and astrometry, then merges into a unified period-mass
sensitivity grid showing which parameter space each method constrains.

References:
    - Winn (2010), transit detection SNR formalism
    - Cumming (2004), RV detection limits
    - Kervella et al. (2019), PMa companion mass sensitivity
"""

import logging
import numpy as np

logger = logging.getLogger(__name__)

# Physical constants
R_SUN_M = 6.957e8       # Solar radius in meters
R_EARTH_M = 6.371e6     # Earth radius in meters
R_EARTH_RSUN = R_EARTH_M / R_SUN_M  # Earth radius in solar radii
M_SUN_KG = 1.989e30
M_EARTH_KG = 5.972e24
M_JUP_KG = 1.898e27
AU_M = 1.496e11


def transit_sensitivity(stellar_radius_rsun, cdpp_ppm, observation_baseline_days,
                         period_grid=None, snr_threshold=7.1,
                         geometric_prob=True):
    """Compute minimum detectable planet radius vs orbital period from transits.

    Uses the matched-filter SNR approach: transit detection requires
    SNR = (depth / CDPP) * sqrt(n_transits * n_points_in_transit) > threshold.

    Simplified to: SNR ~ (Rp/R*)^2 * 1e6 / CDPP * sqrt(n_transits).
    Solving for Rp: Rp_min = R* * sqrt(CDPP * snr_threshold / (1e6 * sqrt(n_transits))).

    Parameters
    ----------
    stellar_radius_rsun : float
        Stellar radius in solar radii.
    cdpp_ppm : float
        Combined Differential Photometric Precision in ppm. Typical TESS CDPP
        is 50-200 ppm for bright stars at 2-min cadence.
    observation_baseline_days : float
        Total observation baseline in days.
    period_grid : array-like, optional
        Periods (days) at which to compute sensitivity.
        Default: log-spaced 0.5 to 1000 days (100 points).
    snr_threshold : float
        Minimum SNR for transit detection. Default 7.1 (Kepler-like).
    geometric_prob : bool
        If True, include geometric transit probability (R*/a) as a weighting
        factor on the detection probability.

    Returns
    -------
    dict
        Sensitivity results with keys: period_grid (array),
        min_radius_rearth (array), min_radius_rjup (array),
        n_transits (array), geometric_probability (array),
        detection_probability (array -- geometric_prob * (detected? 1:0)).
    """
    if period_grid is None:
        period_grid = np.logspace(np.log10(0.5), 3.0, 100)
    period_grid = np.asarray(period_grid, dtype=float)

    n_transits = np.floor(observation_baseline_days / period_grid).astype(float)
    n_transits = np.maximum(n_transits, 1.0)

    # Minimum detectable depth (ppm) for given SNR and n_transits
    # SNR = depth_ppm / cdpp_ppm * sqrt(n_transits)
    # => depth_min = cdpp_ppm * snr_threshold / sqrt(n_transits)
    depth_min_ppm = cdpp_ppm * snr_threshold / np.sqrt(n_transits)

    # Convert depth to planet radius
    # depth = (Rp/R*)^2 => Rp = R* * sqrt(depth)
    depth_min_frac = depth_min_ppm / 1e6
    # Guard against negative or zero depth
    depth_min_frac = np.maximum(depth_min_frac, 0.0)
    rp_rsun = stellar_radius_rsun * np.sqrt(depth_min_frac)
    rp_rearth = rp_rsun / R_EARTH_RSUN
    rp_rjup = rp_rearth * 0.08921  # Rearth to Rjup

    # Mark periods where fewer than 2 transits observed as undetectable
    too_few = n_transits < 2
    rp_rearth[too_few] = np.inf
    rp_rjup[too_few] = np.inf

    # Geometric transit probability: P_geo = R* / a
    # a from Kepler's 3rd law: a(AU) = P(yr)^(2/3) (for solar mass)
    # We use a simplification for solar-like stars
    p_years = period_grid / 365.25
    a_au = p_years ** (2.0 / 3.0)  # approximate for ~1 Msun
    a_rsun = a_au * AU_M / R_SUN_M
    geo_prob = np.minimum(stellar_radius_rsun / a_rsun, 1.0)

    if not geometric_prob:
        geo_prob = np.ones_like(period_grid)

    logger.info("Transit sensitivity: R*=%.3f Rsun, CDPP=%.0f ppm, "
                "baseline=%.0f d, SNR_thresh=%.1f",
                stellar_radius_rsun, cdpp_ppm, observation_baseline_days,
                snr_threshold)

    return {
        "period_grid": period_grid,
        "min_radius_rearth": rp_rearth,
        "min_radius_rjup": rp_rjup,
        "n_transits": n_transits,
        "geometric_probability": geo_prob,
        "cdpp_ppm": cdpp_ppm,
        "snr_threshold": snr_threshold,
    }


def rv_sensitivity(time, rv_err, stellar_mass_msun, period_grid=None):
    """Compute minimum detectable planet mass vs period from RV data.

    Wraps rv_detection_limit() and rv_to_planet_mass() to produce a
    sensitivity curve in (period, mass) space.

    Parameters
    ----------
    time : array-like
        BJD timestamps of actual RV observations.
    rv_err : array-like
        RV uncertainties (m/s).
    stellar_mass_msun : float
        Stellar mass in solar masses.
    period_grid : array-like, optional
        Periods (days). Default: log-spaced 1 to 10000 days.

    Returns
    -------
    dict
        Sensitivity results with keys: period_grid, k_min_ms,
        mass_min_mearth, mass_min_mjup, sigma_eff_ms, n_obs, baseline_days.
    """
    from src.rv_data import rv_detection_limit, rv_to_planet_mass

    det_limit = rv_detection_limit(time, rv_err, period_grid=period_grid)

    periods = det_limit["periods"]
    k_min = det_limit["k_min_ms"]

    # Convert K_min to planet mass at each period
    mass_mearth = np.zeros_like(k_min)
    for i in range(len(periods)):
        if np.isfinite(k_min[i]):
            mass_mearth[i] = rv_to_planet_mass(k_min[i], periods[i],
                                                stellar_mass_msun=stellar_mass_msun)
        else:
            mass_mearth[i] = np.inf

    mass_mjup = mass_mearth / 317.8

    logger.info("RV sensitivity: M*=%.2f Msun, sigma_eff=%.2f m/s, "
                "%d obs, %.0f d baseline",
                stellar_mass_msun, det_limit["sigma_eff_ms"],
                det_limit["n_obs"], det_limit["baseline_days"])

    return {
        "period_grid": periods,
        "k_min_ms": k_min,
        "mass_min_mearth": mass_mearth,
        "mass_min_mjup": mass_mjup,
        "sigma_eff_ms": det_limit["sigma_eff_ms"],
        "n_obs": det_limit["n_obs"],
        "baseline_days": det_limit["baseline_days"],
    }


def astrometric_sensitivity(distance_pc, stellar_mass_msun,
                              pma_detection_limit_mas_yr=None,
                              ruwe=None, separation_grid_au=None):
    """Compute companion mass constraints from Gaia astrometry.

    Uses the PMa formalism: at each assumed separation, what companion mass
    would produce a detectable PMa signal? Also considers RUWE as an indicator.

    Parameters
    ----------
    distance_pc : float
        Distance to the star in parsecs.
    stellar_mass_msun : float
        Stellar mass in solar masses.
    pma_detection_limit_mas_yr : float, optional
        Minimum detectable PMa in mas/yr (3-sigma). Default: 0.5 mas/yr.
    ruwe : float, optional
        RUWE from Gaia DR3. If > 1.4, indicates possible companion.
    separation_grid_au : array-like, optional
        Separation grid in AU. Default: log-spaced 0.5 to 100 AU.

    Returns
    -------
    dict
        Sensitivity results with keys: separation_au (array),
        mass_min_mearth (array), mass_min_mjup (array),
        period_grid_days (array -- from Kepler's law), ruwe_hint.
    """
    from src.proper_motion import pma_companion_mass, EPOCH_DIFF_YR

    if pma_detection_limit_mas_yr is None:
        pma_detection_limit_mas_yr = 0.5  # typical for bright nearby stars

    if separation_grid_au is None:
        separation_grid_au = np.logspace(np.log10(0.5), 2.0, 50)
    separation_grid_au = np.asarray(separation_grid_au, dtype=float)

    # Compute minimum companion mass at each separation
    mass_est = pma_companion_mass(
        pma_detection_limit_mas_yr, distance_pc, stellar_mass_msun,
        separation_au=separation_grid_au,
    )

    # Convert separations to periods via Kepler's 3rd law
    # P^2 = a^3 / M_star  (in years and AU)
    period_years = np.sqrt(separation_grid_au**3 / stellar_mass_msun)
    period_days = period_years * 365.25

    ruwe_hint = None
    if ruwe is not None:
        ruwe_hint = "elevated" if ruwe > 1.4 else "normal"

    logger.info("Astrometric sensitivity: d=%.1f pc, M*=%.2f Msun, "
                "PMa_min=%.2f mas/yr, RUWE=%s",
                distance_pc, stellar_mass_msun, pma_detection_limit_mas_yr,
                f"{ruwe:.2f}" if ruwe else "N/A")

    return {
        "separation_au": separation_grid_au,
        "mass_min_mearth": mass_est["mass_mearth"],
        "mass_min_mjup": mass_est["mass_mjup"],
        "period_grid_days": period_days,
        "pma_detection_limit_mas_yr": pma_detection_limit_mas_yr,
        "ruwe": ruwe,
        "ruwe_hint": ruwe_hint,
    }


def combined_sensitivity(transit_sens=None, rv_sens=None, astro_sens=None,
                          period_grid=None):
    """Merge transit, RV, and astrometric sensitivity into a unified map.

    Interpolates each method's sensitivity onto a common period grid and
    determines which regions of (period, mass/radius) space are detectable
    by which method(s).

    Parameters
    ----------
    transit_sens : dict, optional
        Output from transit_sensitivity().
    rv_sens : dict, optional
        Output from rv_sensitivity().
    astro_sens : dict, optional
        Output from astrometric_sensitivity().
    period_grid : array-like, optional
        Common period grid (days). Default: log-spaced 0.5 to 10000 days.

    Returns
    -------
    dict
        Combined sensitivity with keys: period_grid (array),
        rv_mass_min_mearth (array or None), transit_radius_min_rearth (array or None),
        astro_mass_min_mearth (array or None), methods_available (list),
        best_method_per_period (list of str).
    """
    if period_grid is None:
        period_grid = np.logspace(np.log10(0.5), np.log10(10000), 200)
    period_grid = np.asarray(period_grid, dtype=float)

    result = {
        "period_grid": period_grid,
        "methods_available": [],
    }

    # Interpolate RV sensitivity
    if rv_sens is not None:
        rv_mass = np.interp(
            period_grid, rv_sens["period_grid"], rv_sens["mass_min_mearth"],
            left=np.inf, right=np.inf,
        )
        result["rv_mass_min_mearth"] = rv_mass
        result["methods_available"].append("rv")
        logger.info("Combined: RV sensitivity interpolated (%d -> %d periods)",
                     len(rv_sens["period_grid"]), len(period_grid))
    else:
        result["rv_mass_min_mearth"] = None

    # Interpolate transit sensitivity (radius, not mass)
    if transit_sens is not None:
        transit_radius = np.interp(
            period_grid, transit_sens["period_grid"],
            transit_sens["min_radius_rearth"],
            left=np.inf, right=np.inf,
        )
        result["transit_radius_min_rearth"] = transit_radius
        result["methods_available"].append("transit")
        logger.info("Combined: transit sensitivity interpolated")
    else:
        result["transit_radius_min_rearth"] = None

    # Interpolate astrometric sensitivity
    if astro_sens is not None:
        astro_mass = np.interp(
            period_grid, astro_sens["period_grid_days"],
            astro_sens["mass_min_mearth"],
            left=np.inf, right=np.inf,
        )
        result["astro_mass_min_mearth"] = astro_mass
        result["methods_available"].append("astrometry")
        logger.info("Combined: astrometric sensitivity interpolated")
    else:
        result["astro_mass_min_mearth"] = None

    # Determine best mass-sensitive method per period
    # (transit provides radius, not mass, so excluded from mass comparison)
    best_method = []
    for k in range(len(period_grid)):
        candidates = {}
        if rv_sens is not None and np.isfinite(result["rv_mass_min_mearth"][k]):
            candidates["rv"] = result["rv_mass_min_mearth"][k]
        if astro_sens is not None and np.isfinite(result["astro_mass_min_mearth"][k]):
            candidates["astrometry"] = result["astro_mass_min_mearth"][k]

        if candidates:
            best = min(candidates, key=candidates.get)
            best_method.append(best)
        else:
            best_method.append("none")

    result["best_mass_method_per_period"] = best_method

    logger.info("Combined sensitivity: %d methods, %d period points",
                len(result["methods_available"]), len(period_grid))

    return result
