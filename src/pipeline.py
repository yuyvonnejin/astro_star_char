"""Pipeline orchestration: chains distance -> temperature/luminosity -> mass."""

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


def process_star(star_dict):
    """Run the full pipeline on a single star record.

    Parameters
    ----------
    star_dict : dict
        Input record matching the input schema from docs.md.

    Returns
    -------
    dict
        Combined output from all three modules.
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

    logger.info("Finished star %s", source_id)
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
