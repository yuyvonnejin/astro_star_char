"""Run the stellar property pipeline on predefined or custom stars."""

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.pipeline import process_star
from src.data_access import resolve_simbad_name, query_stars_by_id

LOG_DIR = Path(__file__).resolve().parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "run_stars.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# -- Predefined test stars --------------------------------------------------

STARS = {
    "sun": {
        "source_id": "sun_synthetic",
        "parallax_mas": 100.0,
        "parallax_error_mas": 0.01,
        "phot_g_mean_mag": 4.83,
        "phot_bp_mean_mag": 5.24,
        "phot_rp_mean_mag": 4.42,
        "bp_rp": 0.82,
        "ag_gspphot": 0.0,
        "ebpminrp_gspphot": 0.0,
        "feh": 0.0,
        "logg": 4.44,
        "is_cepheid": False,
        "cepheid_period_days": None,
        "teff_gspphot": 5772.0,
        "lum_gspphot": 1.0,
    },
    "proxima_cen": {
        "source_id": "5853498713190525696",
        "parallax_mas": 768.07,
        "parallax_error_mas": 0.03,
        "phot_g_mean_mag": 11.13,
        "phot_bp_mean_mag": 12.95,
        "phot_rp_mean_mag": 9.45,
        "bp_rp": 3.50,
        "ag_gspphot": 0.01,
        "ebpminrp_gspphot": 0.005,
        "feh": 0.0,
        "logg": 4.5,
        "is_cepheid": False,
        "cepheid_period_days": None,
        "teff_gspphot": None,
        "lum_gspphot": None,
    },
    "sirius_a": {
        "source_id": "sirius_a_synthetic",
        "parallax_mas": 379.21,
        "parallax_error_mas": 1.58,
        "phot_g_mean_mag": -1.09,
        "phot_bp_mean_mag": -1.17,
        "phot_rp_mean_mag": -0.81,
        "bp_rp": -0.36,
        "ag_gspphot": 0.0,
        "ebpminrp_gspphot": 0.0,
        "feh": 0.0,
        "logg": 4.3,
        "is_cepheid": False,
        "cepheid_period_days": None,
        "teff_gspphot": 9940.0,
        "lum_gspphot": 25.4,
    },
    "delta_cep": {
        "source_id": "delta_cep_synthetic",
        "parallax_mas": 3.66,
        "parallax_error_mas": 0.15,
        "phot_g_mean_mag": 3.95,
        "phot_bp_mean_mag": 4.30,
        "phot_rp_mean_mag": 3.40,
        "bp_rp": 0.90,
        "ag_gspphot": 0.25,
        "ebpminrp_gspphot": 0.10,
        "feh": 0.0,
        "logg": 2.0,
        "is_cepheid": True,
        "cepheid_period_days": 5.37,
        "teff_gspphot": 5900.0,
        "lum_gspphot": 2000.0,
    },
    "alpha_cen_a": {
        "source_id": "alpha_cen_a_synthetic",
        "parallax_mas": 747.17,
        "parallax_error_mas": 1.33,
        "phot_g_mean_mag": -0.01,
        "phot_bp_mean_mag": 0.39,
        "phot_rp_mean_mag": -0.52,
        "bp_rp": 0.91,
        "ag_gspphot": 0.0,
        "ebpminrp_gspphot": 0.0,
        "feh": 0.20,
        "logg": 4.30,
        "is_cepheid": False,
        "cepheid_period_days": None,
        "teff_gspphot": 5790.0,
        "lum_gspphot": 1.52,
    },
    "barnards_star": {
        "source_id": "4472832130942575872",
        "parallax_mas": 546.98,
        "parallax_error_mas": 0.04,
        "phot_g_mean_mag": 9.51,
        "phot_bp_mean_mag": 11.24,
        "phot_rp_mean_mag": 8.19,
        "bp_rp": 3.05,
        "ag_gspphot": 0.0,
        "ebpminrp_gspphot": 0.0,
        "feh": -0.39,
        "logg": 5.05,
        "is_cepheid": False,
        "cepheid_period_days": None,
        "teff_gspphot": 3278.0,
        "lum_gspphot": 0.0035,
    },
}


def format_result(name, result):
    """Format a single star result for display."""
    lines = []
    lines.append(f"  {'Star':20s}: {name}")
    lines.append(f"  {'Source ID':20s}: {result.get('source_id', '?')}")

    d = result.get("distance_pc")
    method = result.get("distance_method", "?")
    if d is not None:
        lo = result.get("distance_lower_pc")
        hi = result.get("distance_upper_pc")
        interval = f" [{lo:.2f}, {hi:.2f}]" if lo and hi else ""
        lines.append(f"  {'Distance':20s}: {d:.4f} pc{interval}  ({method})")
    else:
        lines.append(f"  {'Distance':20s}: FAILED")

    teff = result.get("teff_K")
    if teff is not None:
        flag = result.get("teff_flag", "")
        unc = result.get("teff_uncertainty_K", "?")
        lines.append(f"  {'Temperature':20s}: {teff:.0f} +/- {unc} K  ({flag})")

    lum = result.get("luminosity_Lsun")
    if lum is not None:
        lines.append(f"  {'Luminosity':20s}: {lum:.5f} Lsun")
        ratio = result.get("luminosity_validation_ratio")
        if ratio is not None:
            lines.append(f"  {'  validation ratio':20s}: {ratio:.3f}")
    else:
        lines.append(f"  {'Luminosity':20s}: N/A")

    radius = result.get("radius_Rsun")
    if radius is not None:
        lines.append(f"  {'Radius':20s}: {radius:.4f} Rsun")
    else:
        lines.append(f"  {'Radius':20s}: N/A")

    mass = result.get("mass_Msun")
    if mass is not None:
        lines.append(f"  {'Mass':20s}: {mass:.3f} Msun")
    else:
        flag = result.get("mass_flag", "?")
        lines.append(f"  {'Mass':20s}: N/A ({flag})")

    ms = result.get("is_main_sequence")
    if ms is not None:
        lines.append(f"  {'Main sequence':20s}: {ms}")

    return "\n".join(lines)


def resolve_and_query(simbad_name):
    """Resolve a SIMBAD name to Gaia data and return (display_name, star_dict) or None."""
    source_id = resolve_simbad_name(simbad_name)
    if source_id is None:
        logger.warning("Could not resolve '%s'", simbad_name)
        return None

    stars = query_stars_by_id([source_id])
    if not stars:
        logger.warning("Gaia query returned no data for source_id %s", source_id)
        return None

    return simbad_name, stars[0]


def run(predefined_names=None, simbad_names=None):
    """Run pipeline on predefined and/or SIMBAD-resolved stars."""
    print("=" * 60)
    print("Stellar Property Pipeline")
    print("=" * 60)

    results = {}

    # Predefined stars
    if predefined_names is None and simbad_names is None:
        predefined_names = list(STARS.keys())

    for name in (predefined_names or []):
        star_data = STARS.get(name)
        if star_data is None:
            logger.warning("Unknown predefined star: %s (available: %s)", name, ", ".join(STARS.keys()))
            continue
        print(f"\n--- {name} ---")
        result = process_star(star_data)
        results[name] = result
        print(format_result(name, result))

    # SIMBAD-resolved stars
    for simbad_name in (simbad_names or []):
        print(f"\n--- {simbad_name} (SIMBAD lookup) ---")
        resolved = resolve_and_query(simbad_name)
        if resolved is None:
            print(f"  FAILED to resolve '{simbad_name}'")
            continue
        display_name, star_data = resolved
        result = process_star(star_data)
        results[display_name] = result
        print(format_result(display_name, result))

    # Save full JSON output
    output_path = Path(__file__).resolve().parent / "output" / "results.json"
    output_path.parent.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Full results saved to %s", output_path)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stellar Property Pipeline")
    parser.add_argument("predefined", nargs="*", help="Predefined star keys (e.g. sun, proxima_cen)")
    parser.add_argument("--name", nargs="+", dest="simbad_names",
                        help="SIMBAD names to resolve (e.g. Vega, 'Alp Lyr', 'HD 172167')")
    args = parser.parse_args()

    predefined = args.predefined if args.predefined else None
    run(predefined_names=predefined, simbad_names=args.simbad_names)
