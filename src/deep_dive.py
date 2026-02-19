"""Phase 7b: Deep-dive orchestrator for single-target comprehensive analysis.

Chains all available analysis methods on a single target star:
1. Stellar properties via Gaia DR3 (Modules 1-3)
2. TESS light curve retrieval + transit detection (Modules 4-5)
3. DACE RV time series + multi-planet residual analysis
4. Injection-recovery RV sensitivity
5. Known planet query (NASA Exoplanet Archive)
6. Gaia-Hipparcos proper motion anomaly
7. Combined sensitivity map (transit + RV + astrometry)

Produces a structured result dict and saves JSON + markdown reports.
"""

import json
import logging
from datetime import datetime
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output" / "target_reports"


def run_deep_dive(target_name, max_tess_sectors=26,
                   injection_n_trials=50, injection_n_periods=15,
                   injection_n_amplitudes=10,
                   save_report=True):
    """Run comprehensive deep-dive analysis on a single target star.

    Parameters
    ----------
    target_name : str
        Target name (e.g., '82 G. Eridani', 'HD 20794').
    max_tess_sectors : int
        Maximum TESS sectors to download. Set to 0 for all available.
    injection_n_trials : int
        Number of noise realizations per injection-recovery grid point.
    injection_n_periods : int
        Number of periods in injection-recovery grid.
    injection_n_amplitudes : int
        Number of K amplitudes in injection-recovery grid.
    save_report : bool
        If True, save JSON and markdown reports to output/target_reports/.

    Returns
    -------
    dict
        Comprehensive analysis result with sections: target, stellar_properties,
        tess_lightcurve, transit_search, known_planets, rv_data, rv_residual_analysis,
        injection_recovery, proper_motion_anomaly, sensitivity_map, summary.
    """
    from src.targets import get_target, resolve_target_ids

    logger.info("=" * 70)
    logger.info("DEEP DIVE: %s", target_name)
    logger.info("=" * 70)

    # Step 1: Target lookup and ID resolution
    target = get_target(target_name)
    if target is None:
        logger.error("Target '%s' not found in catalog", target_name)
        return {"error": f"Target '{target_name}' not found in catalog"}

    target = resolve_target_ids(target)
    logger.info("Target: %s (HD=%s, HIP=%d, Gaia=%s, TIC=%s)",
                target["name"], target["hd"], target["hip"],
                target.get("gaia_dr3_id"), target.get("tic"))

    result = {
        "target": _target_info(target),
        "generated_utc": datetime.utcnow().isoformat(),
        "analysis_version": "phase7b",
    }

    # Step 2: Stellar properties via pipeline (Modules 1-3)
    result["stellar_properties"] = _run_stellar_properties(target)

    # Extract key stellar params for later use
    sp = result["stellar_properties"]
    stellar_mass = sp.get("mass_Msun") or target.get("mass_msun", 0.7)
    stellar_radius = sp.get("radius_Rsun") or 0.8  # fallback for G6V
    stellar_teff = sp.get("teff_K") or target.get("teff_k")
    stellar_lum = sp.get("luminosity_Lsun")
    distance_pc = sp.get("distance_pc") or target.get("distance_pc")

    # Step 3: Known planets query
    result["known_planets"] = _run_known_planets(target)
    known_periods = _extract_known_periods(result["known_planets"])

    # Step 4: TESS light curve + transit search
    tess_result = _run_tess_analysis(target, stellar_mass, stellar_radius,
                                      stellar_teff, stellar_lum,
                                      max_sectors=max_tess_sectors)
    result["tess_lightcurve"] = tess_result.get("lightcurve")
    result["transit_search"] = tess_result.get("transit")

    # Step 5: RV data retrieval and multi-planet analysis
    rv_result = _run_rv_analysis(target, known_periods)
    result["rv_data"] = rv_result.get("rv_data")
    result["rv_residual_analysis"] = rv_result.get("residual_analysis")

    # Step 6: Injection-recovery sensitivity
    result["injection_recovery"] = _run_injection_recovery(
        rv_result.get("rv_raw"), stellar_mass,
        n_trials=injection_n_trials,
        n_periods=injection_n_periods,
        n_amplitudes=injection_n_amplitudes,
    )

    # Step 7: Proper motion anomaly
    result["proper_motion_anomaly"] = _run_pma(target)

    # Step 8: Combined sensitivity map
    result["sensitivity_map"] = _run_sensitivity_map(
        stellar_radius, stellar_mass, distance_pc,
        rv_result.get("rv_raw"), tess_result,
        result.get("proper_motion_anomaly"),
    )

    # Summary
    result["summary"] = _build_summary(result)

    # Save report
    if save_report:
        result["saved_files"] = _save_deep_dive_report(result)

    logger.info("=" * 70)
    logger.info("DEEP DIVE COMPLETE: %s", target_name)
    logger.info("=" * 70)

    return result


# --- Section runners ---

def _target_info(target):
    """Extract target info for the report."""
    return {
        "name": target["name"],
        "hd": target.get("hd"),
        "hip": target.get("hip"),
        "tic": target.get("tic"),
        "gaia_dr3_id": target.get("gaia_dr3_id"),
        "spectral_type": target.get("spectral_type"),
        "teff_k": target.get("teff_k"),
        "distance_pc": target.get("distance_pc"),
        "mass_msun": target.get("mass_msun"),
    }


def _run_stellar_properties(target):
    """Run Modules 1-3 stellar properties via Gaia DR3."""
    logger.info("--- Stellar Properties (Modules 1-3) ---")
    try:
        from src.data_access import resolve_simbad_name, query_stars_by_id
        from src.pipeline import process_star

        gaia_id = target.get("gaia_dr3_id")
        if gaia_id is None:
            gaia_id = resolve_simbad_name(target["name"])
            if gaia_id is None and target.get("hd"):
                gaia_id = resolve_simbad_name(target["hd"])

        if gaia_id is None:
            logger.warning("No Gaia ID resolved for %s", target["name"])
            return {"status": "no_gaia_id"}

        stars = query_stars_by_id([gaia_id])
        if not stars:
            return {"status": "gaia_query_failed"}

        star_dict = stars[0]
        result = process_star(star_dict, include_lightcurve=False,
                               include_transit=False)
        result["status"] = "ok"
        return result

    except Exception as e:
        logger.error("Stellar properties failed: %s", e)
        return {"status": "error", "error": str(e)}


def _run_known_planets(target):
    """Query NASA Exoplanet Archive for known planets.

    Falls back to hardcoded reference data for key targets if the
    archive is unreachable (timeout/DNS issues are common).
    """
    logger.info("--- Known Planets ---")
    try:
        from src.rv_data import query_known_planets
        # Use common name (not HD) since _generate_search_names adds HD
        # from aliases, giving better coverage with one call
        planets = query_known_planets(target["name"])
        if not planets:
            planets = query_known_planets(target["hd"])
        if planets:
            return {
                "n_planets": len(planets),
                "planets": planets,
                "status": "ok",
            }
    except Exception as e:
        logger.warning("Known planets query failed: %s", e)

    # Fallback: hardcoded reference data for key targets
    fallback = _known_planets_fallback(target)
    if fallback:
        logger.info("Using hardcoded reference data for %s (%d planets)",
                     target["name"], len(fallback))
        return {
            "n_planets": len(fallback),
            "planets": fallback,
            "status": "fallback",
        }

    return {"n_planets": 0, "planets": [], "status": "none_found"}


def _known_planets_fallback(target):
    """Return hardcoded planet data for key targets when NASA is unreachable.

    Reference data from peer-reviewed publications (NASA Exoplanet Archive
    composite values, retrieved 2026-02). Only includes targets from the
    Phase 7 catalog that have confirmed planets.
    """
    # Key: HD number -> list of planet dicts
    reference = {
        "HD 20794": [
            {"pl_name": "HD 20794 b", "hostname": "HD 20794",
             "period_days": 18.315, "mass_earth": 2.7, "mass_jupiter": 0.0085,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2011,
             "semi_major_axis_au": 0.1207, "eq_temp_k": None,
             "eccentricity": 0.0, "distance_pc": 6.04},
            {"pl_name": "HD 20794 d", "hostname": "HD 20794",
             "period_days": 89.694, "mass_earth": 3.53, "mass_jupiter": 0.0111,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2011,
             "semi_major_axis_au": 0.3499, "eq_temp_k": None,
             "eccentricity": 0.0, "distance_pc": 6.04},
            {"pl_name": "HD 20794 e", "hostname": "HD 20794",
             "period_days": 147.025, "mass_earth": 4.77, "mass_jupiter": 0.015,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2011,
             "semi_major_axis_au": 0.5095, "eq_temp_k": None,
             "eccentricity": 0.0, "distance_pc": 6.04},
            {"pl_name": "HD 20794 f", "hostname": "HD 20794",
             "period_days": 647.6, "mass_earth": 5.35, "mass_jupiter": 0.0168,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2017,
             "semi_major_axis_au": 1.37, "eq_temp_k": None,
             "eccentricity": 0.0, "distance_pc": 6.04},
        ],
        "HD 10700": [
            {"pl_name": "tau Cet e", "hostname": "HD 10700",
             "period_days": 162.87, "mass_earth": 3.93, "mass_jupiter": 0.0124,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2012,
             "semi_major_axis_au": 0.538, "eq_temp_k": None,
             "eccentricity": 0.18, "distance_pc": 3.65},
            {"pl_name": "tau Cet f", "hostname": "HD 10700",
             "period_days": 636.13, "mass_earth": 3.93, "mass_jupiter": 0.0124,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2012,
             "semi_major_axis_au": 1.334, "eq_temp_k": None,
             "eccentricity": 0.16, "distance_pc": 3.65},
        ],
        "HD 115617": [
            {"pl_name": "61 Vir b", "hostname": "HD 115617",
             "period_days": 4.215, "mass_earth": 5.1, "mass_jupiter": 0.016,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2009,
             "semi_major_axis_au": 0.050, "eq_temp_k": None,
             "eccentricity": 0.12, "distance_pc": 8.56},
            {"pl_name": "61 Vir c", "hostname": "HD 115617",
             "period_days": 38.021, "mass_earth": 18.2, "mass_jupiter": 0.057,
             "radius_earth": None, "radius_jupiter": None,
             "discovery_method": "Radial Velocity", "discovery_year": 2009,
             "semi_major_axis_au": 0.217, "eq_temp_k": None,
             "eccentricity": 0.14, "distance_pc": 8.56},
        ],
    }

    hd = target.get("hd")
    if hd and hd in reference:
        return reference[hd]
    return None


def _extract_known_periods(known_planets_result):
    """Extract orbital periods from known planets result."""
    periods = []
    if known_planets_result and known_planets_result.get("planets"):
        for planet in known_planets_result["planets"]:
            p = planet.get("period_days")
            if p is not None and p > 0:
                periods.append(float(p))
    if periods:
        logger.info("Known planet periods: %s days",
                     ", ".join(f"{p:.2f}" for p in periods))
    return periods


def _run_tess_analysis(target, stellar_mass, stellar_radius, stellar_teff,
                        stellar_lum, max_sectors=26):
    """Run TESS light curve retrieval and transit search."""
    logger.info("--- TESS Light Curve + Transit Search ---")
    try:
        from src.lightcurve import (
            search_lightcurve, download_and_stitch,
            clean_lightcurve, bin_lightcurve, flatten_lightcurve,
        )
        from src.periodogram import analyze_lightcurve
        from src.transit import analyze_transit
        import lightkurve as lk

        # Search names in priority order
        search_names = []
        if target.get("tic"):
            search_names.append(f"TIC {target['tic']}")
        search_names.append(target["hd"])
        search_names.append(target["name"])

        # Use author='SPOC' to get consistent 2-min cadence data only.
        # Mixing authors (SPOC, QLP, TESS-SPOC) during stitch corrupts the
        # light curve because they use different flux normalization.
        lc_data = None
        for search_name in search_names:
            logger.info("Trying TESS SPOC light curve for '%s'", search_name)
            search_info = search_lightcurve(search_name, mission="TESS")
            if search_info is not None:
                # Filter to SPOC author for consistent flux format
                search_res = search_info["search_result"]
                if hasattr(search_res, "author"):
                    spoc_mask = search_res.author == "SPOC"
                    if np.any(spoc_mask):
                        filtered = search_res[spoc_mask]
                        search_info = {
                            "search_result": filtered,
                            "mission": "TESS",
                            "author": "SPOC",
                            "n_available": len(filtered),
                        }
                        logger.info("Filtered to %d SPOC sectors (of %d total)",
                                     len(filtered), len(search_res))

                lc_data = download_and_stitch(search_info, max_sectors=max_sectors)
                if lc_data is not None:
                    lc_data = clean_lightcurve(lc_data)
                    lc_data = bin_lightcurve(lc_data, bin_cadence_s=120.0)
                    lc_data = flatten_lightcurve(lc_data)
                    logger.info("TESS data found via '%s'", search_name)
                    break

        if lc_data is None:
            logger.info("No TESS light curve data available for %s", target["name"])
            return {
                "lightcurve": {"available": False, "status": "no_data"},
                "transit": {"transit_detected": None, "status": "no_lightcurve"},
            }

        # Variability analysis (Module 4)
        lc_analysis = analyze_lightcurve(lc_data)

        # Compute actual baseline from binned time array
        actual_baseline = float(lc_data["time"][-1] - lc_data["time"][0])

        lc_summary = {
            "available": True,
            "n_sectors": lc_data.get("n_sectors"),
            "n_points": lc_data.get("n_points_clean", len(lc_data.get("time", []))),
            "baseline_days": round(actual_baseline, 2),
            "cadence_s": lc_data.get("cadence_s"),
            "variability_class": lc_analysis.get("variability_class"),
            "amplitude_ppt": lc_analysis.get("amplitude_ppt"),
            "period_days": lc_analysis.get("period_days"),
            "status": "ok",
        }

        # Transit detection (Module 5)
        stellar_props = {
            "radius_Rsun": stellar_radius,
            "mass_Msun": stellar_mass,
            "teff_K": stellar_teff,
            "luminosity_Lsun": stellar_lum,
        }

        variability_period = None
        if lc_analysis.get("variability_class") == "periodic":
            variability_period = lc_analysis.get("period_days")

        # Try HZ-targeted first; fall back to full-range if baseline too short
        transit_result = analyze_transit(
            lc_data, stellar_props,
            variability_period=variability_period,
            hz_targeted=True,
            hz_broadening=2.0,
        )

        # If HZ-targeted search failed due to period range, retry full-range
        if (transit_result.get("transit_flag") == "period_range_too_narrow"
                or transit_result.get("transit_detected") is None):
            logger.info("HZ-targeted transit search yielded no result; "
                        "retrying with full period range")
            transit_result = analyze_transit(
                lc_data, stellar_props,
                variability_period=variability_period,
                hz_targeted=False,
            )
            transit_result["search_mode"] = "full_range_fallback"
        else:
            transit_result["search_mode"] = "hz_targeted"

        transit_result["status"] = "ok"

        return {
            "lightcurve": lc_summary,
            "transit": transit_result,
        }

    except Exception as e:
        logger.error("TESS analysis failed: %s", e)
        return {
            "lightcurve": {"available": False, "status": "error", "error": str(e)},
            "transit": {"transit_detected": None, "status": "error", "error": str(e)},
        }


def _run_rv_analysis(target, known_periods):
    """Run RV data retrieval and multi-planet residual analysis."""
    logger.info("--- RV Data + Residual Analysis ---")
    try:
        from src.rv_data import (
            query_dace_rv, rv_periodogram, rv_residual_analysis,
        )

        rv_data = query_dace_rv(target["hd"])
        if rv_data is None:
            logger.info("No DACE RV data for %s", target["name"])
            return {
                "rv_data": {"n_measurements": 0, "status": "no_data"},
                "residual_analysis": None,
                "rv_raw": None,
            }

        rv_summary = {
            "n_measurements": rv_data["n_measurements"],
            "time_baseline_days": rv_data["time_baseline_days"],
            "instruments": rv_data["instruments"],
            "instrument_summary": rv_data.get("instrument_summary"),
            "status": "ok",
        }

        # Run basic periodogram
        pg = rv_periodogram(rv_data["time"], rv_data["rv"], rv_data["rv_err"])
        rv_summary["periodogram"] = {
            "best_period": pg.get("best_period"),
            "best_power": pg.get("best_power"),
            "fap": pg.get("fap"),
            "peaks": pg.get("peaks", []),
        }

        # Run residual analysis if we have known planet periods
        residual = None
        if known_periods:
            residual = rv_residual_analysis(
                rv_data["time"], rv_data["rv"], rv_data["rv_err"],
                known_periods=known_periods,
                instruments=rv_data.get("instrument"),
            )
            # Strip large arrays for storage (keep only summary)
            residual_summary = {
                "known_periods_used": residual["known_periods_used"],
                "offset_correction": residual.get("offset_correction"),
                "sinusoid_subtraction": residual.get("sinusoid_subtraction"),
                "original_best_period": residual["original_periodogram"].get("best_period"),
                "original_fap": residual["original_periodogram"].get("fap"),
                "residual_best_period": residual["residual_periodogram"].get("best_period"),
                "residual_fap": residual["residual_periodogram"].get("fap"),
                "residual_peaks": residual["residual_periodogram"].get("peaks", []),
            }
        else:
            residual_summary = None

        return {
            "rv_data": rv_summary,
            "residual_analysis": residual_summary,
            "rv_raw": rv_data,
        }

    except Exception as e:
        logger.error("RV analysis failed: %s", e)
        return {
            "rv_data": {"n_measurements": 0, "status": "error", "error": str(e)},
            "residual_analysis": None,
            "rv_raw": None,
        }


def _run_injection_recovery(rv_raw, stellar_mass, n_trials=50,
                              n_periods=15, n_amplitudes=10):
    """Run injection-recovery sensitivity test."""
    logger.info("--- RV Injection-Recovery ---")
    if rv_raw is None:
        return {"status": "no_rv_data"}

    try:
        from src.rv_data import rv_injection_recovery

        baseline = rv_raw["time_baseline_days"]
        max_period = min(baseline / 2.0, 3000.0)
        period_grid = np.logspace(np.log10(5.0), np.log10(max_period), n_periods)
        # Start K grid at 0.05 m/s to probe sub-Earth detection threshold
        # for high-precision datasets (ESPRESSO+HARPS with N>10000)
        k_grid = np.logspace(np.log10(0.05), np.log10(10.0), n_amplitudes)

        result = rv_injection_recovery(
            rv_raw["time"], rv_raw["rv_err"],
            period_grid=period_grid, k_grid=k_grid,
            n_trials=n_trials, stellar_mass_msun=stellar_mass,
        )

        # Convert arrays for JSON serialization
        return {
            "period_grid": result["period_grid"].tolist(),
            "k_grid": result["k_grid"].tolist(),
            "detection_probability": result["detection_probability"].tolist(),
            "mass_grid_mearth": result["mass_grid_mearth"].tolist(),
            "n_trials": result["n_trials"],
            "n_obs": result["n_obs"],
            "baseline_days": result["baseline_days"],
            "status": "ok",
        }

    except Exception as e:
        logger.error("Injection-recovery failed: %s", e)
        return {"status": "error", "error": str(e)}


def _run_pma(target):
    """Run proper motion anomaly analysis."""
    logger.info("--- Proper Motion Anomaly ---")
    try:
        from src.proper_motion import analyze_pma

        gaia_id = target.get("gaia_dr3_id")
        hip_id = target.get("hip")

        if gaia_id is None:
            from src.data_access import resolve_simbad_name
            gaia_id = resolve_simbad_name(target["name"])
            if gaia_id is None and target.get("hd"):
                gaia_id = resolve_simbad_name(target["hd"])

        if gaia_id is None or hip_id is None:
            return {"status": "missing_ids"}

        return analyze_pma(
            gaia_id, hip_id,
            distance_pc=target.get("distance_pc"),
            stellar_mass_msun=target.get("mass_msun"),
        )

    except Exception as e:
        logger.error("PMa analysis failed: %s", e)
        return {"status": "error", "error": str(e)}


def _run_sensitivity_map(stellar_radius, stellar_mass, distance_pc,
                          rv_raw, tess_result, pma_result):
    """Compute combined sensitivity map."""
    logger.info("--- Combined Sensitivity Map ---")
    try:
        from src.sensitivity import (
            transit_sensitivity, rv_sensitivity,
            astrometric_sensitivity, combined_sensitivity,
        )

        common_periods = np.logspace(np.log10(0.5), np.log10(10000), 200)

        # Transit sensitivity
        transit_sens = None
        lc_info = tess_result.get("lightcurve") if tess_result else None
        if lc_info and lc_info.get("available"):
            baseline = lc_info.get("baseline_days", 0)
            # Estimate CDPP from light curve scatter (or use typical value)
            cdpp = 100.0  # ppm, conservative for bright star
            if baseline > 0:
                transit_sens = transit_sensitivity(
                    stellar_radius, cdpp, baseline,
                    period_grid=common_periods,
                )

        # RV sensitivity
        rv_sens = None
        if rv_raw is not None:
            rv_sens = rv_sensitivity(
                rv_raw["time"], rv_raw["rv_err"], stellar_mass,
                period_grid=common_periods,
            )

        # Astrometric sensitivity
        astro_sens = None
        if distance_pc is not None and pma_result and pma_result.get("status") == "ok":
            ruwe = pma_result.get("ruwe")
            astro_sens = astrometric_sensitivity(
                distance_pc, stellar_mass,
                ruwe=ruwe,
            )

        # Combined
        combined = combined_sensitivity(
            transit_sens=transit_sens,
            rv_sens=rv_sens,
            astro_sens=astro_sens,
            period_grid=common_periods,
        )

        # Convert arrays for JSON
        json_safe = {}
        for key, val in combined.items():
            if isinstance(val, np.ndarray):
                json_safe[key] = val.tolist()
            else:
                json_safe[key] = val
        json_safe["status"] = "ok"
        return json_safe

    except Exception as e:
        logger.error("Sensitivity map failed: %s", e)
        return {"status": "error", "error": str(e)}


def _build_summary(result):
    """Build a high-level summary of the deep-dive results."""
    summary = {
        "target": result["target"]["name"],
        "findings": [],
        "open_questions": [],
    }

    # Stellar properties
    sp = result.get("stellar_properties", {})
    if sp.get("status") == "ok":
        teff = sp.get("teff_K")
        if teff:
            summary["findings"].append(
                f"Pipeline Teff: {teff:.0f} K (catalog: {result['target'].get('teff_k')} K)"
            )

    # Known planets
    kp = result.get("known_planets", {})
    n_planets = kp.get("n_planets", 0)
    if n_planets > 0:
        names = [p["pl_name"] for p in kp["planets"] if p.get("pl_name")]
        summary["findings"].append(
            f"{n_planets} confirmed planet(s): {', '.join(names)}"
        )

    # Transit search
    tr = result.get("transit_search", {})
    if tr.get("transit_detected"):
        summary["findings"].append(
            f"Transit candidate: P={tr.get('transit_period_days'):.2f} d, "
            f"depth={tr.get('transit_depth_ppm'):.0f} ppm"
        )
    elif tr.get("transit_detected") is False:
        summary["findings"].append("No transit detected in TESS data")

    # RV data
    rv = result.get("rv_data", {})
    if rv.get("n_measurements", 0) > 0:
        summary["findings"].append(
            f"RV data: {rv['n_measurements']} measurements, "
            f"{rv.get('time_baseline_days', 0):.0f} day baseline"
        )

    # RV residuals
    resid = result.get("rv_residual_analysis")
    if resid and resid.get("residual_best_period"):
        summary["findings"].append(
            f"Residual RV best period: {resid['residual_best_period']:.2f} d "
            f"(after subtracting {len(resid.get('known_periods_used', []))} known planets)"
        )
        summary["open_questions"].append(
            f"Is the {resid['residual_best_period']:.2f}-day residual signal "
            "stellar activity, an alias, or a new planet?"
        )

    # PMa
    pma = result.get("proper_motion_anomaly", {})
    pma_data = pma.get("pma", {})
    if pma_data.get("significant"):
        summary["findings"].append(
            f"Significant PMa: {pma_data['pma_total_mas_yr']:.2f} mas/yr "
            f"(SNR={pma_data['pma_snr']:.1f})"
        )
        summary["open_questions"].append(
            "What is the mass and separation of the outer companion "
            "implied by the PMa signal?"
        )

    ruwe = pma.get("ruwe")
    if ruwe and ruwe > 1.4:
        summary["findings"].append(f"Elevated RUWE: {ruwe:.2f}")

    return summary


def _save_deep_dive_report(result):
    """Save deep-dive report as JSON and markdown."""
    output_dir = OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    target_name = result["target"]["name"].replace(" ", "_").lower()
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")

    # JSON report
    json_path = output_dir / f"{target_name}_deep_dive_{timestamp}.json"
    json_str = json.dumps(result, indent=2, default=_json_serializer)
    json_path.write_text(json_str, encoding="utf-8")
    logger.info("Deep-dive JSON report saved: %s", json_path)

    # Markdown report
    md_path = output_dir / f"{target_name}_deep_dive_{timestamp}.md"
    md_str = _format_deep_dive_markdown(result)
    md_path.write_text(md_str, encoding="utf-8")
    logger.info("Deep-dive markdown report saved: %s", md_path)

    return {"json_path": str(json_path), "md_path": str(md_path)}


def _format_deep_dive_markdown(result):
    """Format deep-dive result as markdown report."""
    lines = []
    t = result["target"]
    lines.append(f"# Deep Dive: {t['name']}")
    lines.append(f"**HD**: {t.get('hd')} | **HIP**: {t.get('hip')} | "
                 f"**SpT**: {t.get('spectral_type')}")
    lines.append(f"**Generated**: {result.get('generated_utc')}")
    lines.append("")

    # Summary
    summary = result.get("summary", {})
    if summary.get("findings"):
        lines.append("## Key Findings")
        lines.append("")
        for f in summary["findings"]:
            lines.append(f"- {f}")
        lines.append("")

    if summary.get("open_questions"):
        lines.append("## Open Questions")
        lines.append("")
        for q in summary["open_questions"]:
            lines.append(f"- {q}")
        lines.append("")

    # Stellar properties
    sp = result.get("stellar_properties", {})
    if sp and sp.get("status") == "ok":
        lines.append("## Stellar Properties")
        lines.append("")
        lines.append("| Property | Value |")
        lines.append("|----------|-------|")
        for key in ["distance_pc", "teff_K", "radius_Rsun", "mass_Msun",
                     "luminosity_Lsun"]:
            val = sp.get(key)
            if val is not None:
                lines.append(f"| {key} | {val:.4f} |")
        lines.append("")

    # Known planets
    kp = result.get("known_planets", {})
    if kp.get("planets"):
        lines.append(f"## Known Planets ({kp['n_planets']})")
        lines.append("")
        lines.append("| Planet | Period (d) | Mass (Me) | Method |")
        lines.append("|--------|-----------|-----------|--------|")
        for p in kp["planets"]:
            lines.append(
                f"| {p.get('pl_name', 'N/A')} | "
                f"{p.get('period_days', 'N/A')} | "
                f"{p.get('mass_earth', 'N/A')} | "
                f"{p.get('discovery_method', 'N/A')} |"
            )
        lines.append("")

    # TESS / Transit
    lc = result.get("tess_lightcurve", {})
    if lc:
        lines.append("## TESS Light Curve")
        lines.append("")
        if lc.get("available"):
            lines.append(f"- Sectors: {lc.get('n_sectors')}")
            lines.append(f"- Points: {lc.get('n_points')}")
            lines.append(f"- Baseline: {lc.get('baseline_days')} days")
            lines.append(f"- Variability: {lc.get('variability_class')}")
        else:
            lines.append("No TESS data available.")
        lines.append("")

    tr = result.get("transit_search", {})
    if tr:
        lines.append("## Transit Search")
        lines.append("")
        if tr.get("transit_detected"):
            lines.append(f"- Period: {tr.get('transit_period_days')} days")
            lines.append(f"- Depth: {tr.get('transit_depth_ppm')} ppm")
            lines.append(f"- SDE: {tr.get('transit_sde')}")
        elif tr.get("transit_detected") is False:
            lines.append(f"- No transit detected (best SDE: {tr.get('transit_sde')})")
        else:
            lines.append(f"- Status: {tr.get('status')}")
        lines.append("")

    # RV
    rv = result.get("rv_data", {})
    if rv and rv.get("n_measurements", 0) > 0:
        lines.append("## Radial Velocity Data")
        lines.append("")
        lines.append(f"- Measurements: {rv['n_measurements']}")
        lines.append(f"- Baseline: {rv.get('time_baseline_days')} days")
        lines.append(f"- Instruments: {', '.join(rv.get('instruments', []))}")
        pg = rv.get("periodogram", {})
        if pg:
            lines.append(f"- Best period: {pg.get('best_period')} days "
                         f"(FAP={pg.get('fap')})")
        lines.append("")

    # Residual analysis
    resid = result.get("rv_residual_analysis")
    if resid:
        lines.append("## RV Residual Analysis")
        lines.append("")
        lines.append(f"- Known periods subtracted: {resid.get('known_periods_used')}")
        sub = resid.get("sinusoid_subtraction", {})
        if sub:
            lines.append(f"- RMS before: {sub.get('rms_before_ms')} m/s")
            lines.append(f"- RMS after: {sub.get('rms_after_ms')} m/s")
        lines.append(f"- Residual best period: {resid.get('residual_best_period')} days")
        lines.append(f"- Residual FAP: {resid.get('residual_fap')}")
        lines.append("")

    # PMa
    pma = result.get("proper_motion_anomaly", {})
    pma_data = pma.get("pma", {})
    if pma_data:
        lines.append("## Proper Motion Anomaly")
        lines.append("")
        lines.append(f"- PMa total: {pma_data.get('pma_total_mas_yr')} mas/yr")
        lines.append(f"- SNR: {pma_data.get('pma_snr')}")
        lines.append(f"- Significant: {pma_data.get('significant')}")
        if pma.get("ruwe"):
            lines.append(f"- RUWE: {pma['ruwe']}")
        lines.append("")

    # Sensitivity
    sens = result.get("sensitivity_map", {})
    if sens and sens.get("status") == "ok":
        lines.append("## Sensitivity Map")
        lines.append("")
        lines.append(f"- Methods: {', '.join(sens.get('methods_available', []))}")
        lines.append("")

    return "\n".join(lines)


def _json_serializer(obj):
    """JSON serializer for numpy and other non-standard types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
