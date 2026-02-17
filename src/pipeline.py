"""Pipeline orchestration: chains distance -> temperature/luminosity -> mass -> lightcurve."""

import argparse
import json
import logging
import sys
from pathlib import Path

from src.distance import compute_distance
from src.temperature import compute_temperature_luminosity
from src.mass import compute_mass

LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "pipeline.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# Fields with default fallback values
DEFAULTS = {
    "ag_gspphot": 0.0,
    "ebpminrp_gspphot": 0.0,
    "feh": 0.0,
    "logg": 4.0,
}


def _apply_defaults(star_dict):
    """Fill in missing fields with spec-defined defaults. Returns a new dict."""
    out = dict(star_dict)
    for key, default in DEFAULTS.items():
        if out.get(key) is None:
            out[key] = default
    return out


def process_star(star_dict, include_lightcurve=False, lc_target=None,
                 lc_mission=None, lc_max_sectors=20):
    """Run the full pipeline on a single star record.

    Parameters
    ----------
    star_dict : dict
        Input record matching the input schema from docs.md.
    include_lightcurve : bool
        If True, run Module 4 (light curve analysis via MAST).
    lc_target : str, optional
        Target name for MAST light curve search. If None, not searched.
    lc_mission : str, optional
        Force a specific mission for light curve search (TESS, Kepler, K2).
    lc_max_sectors : int
        Maximum sectors/quarters to download for light curve.

    Returns
    -------
    dict
        Combined output from all modules.
    """
    star = _apply_defaults(star_dict)
    source_id = star.get("source_id", "unknown")
    logger.info("Processing star %s", source_id)

    result = {"source_id": source_id}

    # Forward reference (official) values for comparison
    REF_KEYS = [
        "teff_gspphot", "lum_gspphot",
        "ref_radius_Rsun", "ref_mass_Msun", "ref_distance_pc",
    ]
    for key in REF_KEYS:
        val = star.get(key)
        if val is not None:
            result[key] = val

    # Module 1: Distance
    dist_result = compute_distance(star)
    result.update(dist_result)
    distance_pc = dist_result["distance_pc"]

    if distance_pc is None:
        logger.warning("Distance computation failed for %s; aborting pipeline", source_id)
        return result

    # Module 2: Temperature and luminosity
    temp_result = compute_temperature_luminosity(star, distance_pc)
    result.update(temp_result)

    teff_K = temp_result.get("teff_K")
    luminosity = temp_result.get("luminosity_Lsun")

    # Module 3: Mass
    if teff_K is not None:
        mass_result = compute_mass(star, teff_K, luminosity)
        result.update(mass_result)
    else:
        result["mass_Msun"] = None
        result["mass_flag"] = "no_teff"
        result["is_main_sequence"] = None

    # Module 4: Light curve analysis (optional)
    if include_lightcurve and lc_target is not None:
        lc_result = _run_lightcurve_analysis(lc_target, lc_mission, lc_max_sectors)
        result.update(lc_result)

        # Feed detected Cepheid period back into distance if applicable
        if (lc_result.get("variability_class") == "periodic"
                and lc_result.get("period_days") is not None
                and star.get("is_cepheid")):
            detected_period = lc_result["period_days"]
            logger.info("Detected period %.4f d on Cepheid %s; recomputing distance",
                        detected_period, source_id)
            star_with_period = dict(star)
            star_with_period["cepheid_period_days"] = detected_period
            dist_recomputed = compute_distance(star_with_period)
            result["distance_pc_lc"] = dist_recomputed["distance_pc"]
            result["distance_method_lc"] = dist_recomputed.get("distance_method",
                                                                "cepheid_leavitt")
    elif include_lightcurve:
        result["lightcurve_available"] = False
        logger.info("Light curve requested but no target name provided for %s", source_id)

    logger.info("Finished star %s", source_id)
    return result


def _run_lightcurve_analysis(target, mission, max_sectors):
    """Run Module 4 light curve analysis. Returns result dict."""
    from src.lightcurve import retrieve_lightcurve
    from src.periodogram import analyze_lightcurve

    logger.info("Module 4: searching MAST light curves for '%s'", target)
    lc_data = retrieve_lightcurve(target, mission=mission, max_sectors=max_sectors)

    if lc_data is None:
        logger.info("No light curve data found for '%s'", target)
        return {"lightcurve_available": False}

    result = analyze_lightcurve(lc_data)
    return result


def main():
    parser = argparse.ArgumentParser(description="Astronomy Object Property Pipeline")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--json-file", type=str, help="Path to JSON file (single star or list)")
    group.add_argument("--json-str", type=str, help="Inline JSON string")
    parser.add_argument("--output", type=str, default=None, help="Output file path (default: stdout)")
    args = parser.parse_args()

    if args.json_file:
        with open(args.json_file, "r") as f:
            data = json.load(f)
    else:
        data = json.loads(args.json_str)

    # Accept single star or list
    if isinstance(data, dict):
        data = [data]

    results = []
    for star in data:
        result = process_star(star)
        results.append(result)

    output_json = json.dumps(results, indent=2)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output_json)
        logger.info("Results written to %s", args.output)
    else:
        sys.stdout.write(output_json + "\n")


if __name__ == "__main__":
    main()
