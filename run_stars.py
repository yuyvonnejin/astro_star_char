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
        "ref_radius_Rsun": 1.0,
        "ref_mass_Msun": 1.0,
        "ref_distance_pc": 0.0,  # N/A for synthetic
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
        "teff_gspphot": 3042.0,
        "lum_gspphot": 0.00155,
        "ref_radius_Rsun": 0.1542,
        "ref_mass_Msun": 0.1221,
        "ref_distance_pc": 1.301,
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
        "ref_radius_Rsun": 1.711,
        "ref_mass_Msun": 2.063,
        "ref_distance_pc": 2.637,
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
        "ref_radius_Rsun": 44.5,
        "ref_mass_Msun": 4.5,
        "ref_distance_pc": 273.0,
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
        "ref_radius_Rsun": 1.2234,
        "ref_mass_Msun": 1.1055,
        "ref_distance_pc": 1.339,
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
        "ref_radius_Rsun": 0.187,
        "ref_mass_Msun": 0.144,
        "ref_distance_pc": 1.828,
    },
}


def _fmt_compare(label, computed, ref, unit, comp_fmt=".4f", ref_fmt=".4f"):
    """Format a line showing computed vs reference value with percent error."""
    if computed is not None and ref is not None and ref != 0:
        pct = abs(computed - ref) / abs(ref) * 100
        return f"  {label:20s}: {computed:{comp_fmt}} {unit}  (ref: {ref:{ref_fmt}} {unit}, err: {pct:.1f}%)"
    elif computed is not None:
        return f"  {label:20s}: {computed:{comp_fmt}} {unit}"
    else:
        return f"  {label:20s}: N/A"


def format_result(name, result):
    """Format a single star result for display with reference comparisons."""
    lines = []
    lines.append(f"  {'Star':20s}: {name}")
    lines.append(f"  {'Source ID':20s}: {result.get('source_id', '?')}")

    # Distance
    d = result.get("distance_pc")
    method = result.get("distance_method", "?")
    ref_d = result.get("ref_distance_pc")
    if d is not None:
        lo = result.get("distance_lower_pc")
        hi = result.get("distance_upper_pc")
        interval = f" [{lo:.2f}, {hi:.2f}]" if lo and hi else ""
        base = f"  {'Distance':20s}: {d:.4f} pc{interval}  ({method})"
        if ref_d is not None and ref_d > 0:
            pct = abs(d - ref_d) / abs(ref_d) * 100
            base += f"\n  {'  (reference)':20s}: {ref_d:.4f} pc  (err: {pct:.1f}%)"
        lines.append(base)
    else:
        lines.append(f"  {'Distance':20s}: FAILED")

    # Temperature
    teff = result.get("teff_K")
    ref_teff = result.get("teff_gspphot")
    if teff is not None:
        flag = result.get("teff_flag", "")
        unc = result.get("teff_uncertainty_K", "?")
        base = f"  {'Temperature':20s}: {teff:.0f} +/- {unc} K  ({flag})"
        if ref_teff is not None and ref_teff > 0:
            pct = abs(teff - ref_teff) / ref_teff * 100
            base += f"\n  {'  (reference)':20s}: {ref_teff:.0f} K  (err: {pct:.1f}%)"
        lines.append(base)

    # Luminosity
    lum = result.get("luminosity_Lsun")
    ref_lum = result.get("lum_gspphot")
    if lum is not None:
        base = f"  {'Luminosity':20s}: {lum:.5f} Lsun"
        if ref_lum is not None and ref_lum > 0:
            pct = abs(lum - ref_lum) / ref_lum * 100
            base += f"\n  {'  (reference)':20s}: {ref_lum:.5f} Lsun  (err: {pct:.1f}%)"
        lines.append(base)
    else:
        lines.append(f"  {'Luminosity':20s}: N/A")

    # Radius
    radius = result.get("radius_Rsun")
    ref_radius = result.get("ref_radius_Rsun")
    lines.append(_fmt_compare("Radius", radius, ref_radius, "Rsun"))

    # Mass
    mass = result.get("mass_Msun")
    ref_mass = result.get("ref_mass_Msun")
    if mass is not None:
        lines.append(_fmt_compare("Mass", mass, ref_mass, "Msun", ".3f", ".3f"))
    else:
        flag = result.get("mass_flag", "?")
        ref_str = f"  (ref: {ref_mass:.3f} Msun)" if ref_mass is not None else ""
        lines.append(f"  {'Mass':20s}: N/A ({flag}){ref_str}")

    ms = result.get("is_main_sequence")
    if ms is not None:
        lines.append(f"  {'Main sequence':20s}: {ms}")

    # Module 4: Light curve results
    if result.get("lightcurve_available") is True:
        lines.append(f"  {'--- Light Curve ---':20s}")
        lc_mission = result.get("lc_mission", "?")
        lc_sectors = result.get("lc_n_sectors", "?")
        lc_points = result.get("lc_n_points", "?")
        lc_baseline = result.get("lc_time_baseline_days", "?")
        lines.append(f"  {'LC Mission':20s}: {lc_mission} ({lc_sectors} sectors, {lc_points} points, {lc_baseline} d)")

        var_class = result.get("variability_class", "?")
        var_flag = result.get("variability_flag", "")
        period = result.get("period_days")
        fap = result.get("period_fap")
        amplitude = result.get("amplitude_ppt")

        if period is not None:
            lines.append(f"  {'Period':20s}: {period:.4f} days (FAP: {fap:.2e})")
            lines.append(f"  {'Amplitude':20s}: {amplitude:.2f} ppt")
        lines.append(f"  {'Variability':20s}: {var_class} ({var_flag})")

        # Cepheid distance recomputation
        d_lc = result.get("distance_pc_lc")
        if d_lc is not None:
            method_lc = result.get("distance_method_lc", "?")
            lines.append(f"  {'Distance (LC)':20s}: {d_lc:.4f} pc  ({method_lc})")
    elif result.get("lightcurve_available") is False:
        lines.append(f"  {'Light Curve':20s}: not available")

    # Module 5: Transit detection and planet properties
    if result.get("transit_detected") is True:
        lines.append(f"  {'--- Transit ---':20s}")
        t_flag = result.get("transit_flag", "ok")
        if t_flag != "ok":
            lines.append(f"  {'WARNING':20s}: {t_flag}")
        t_period = result.get("transit_period_days")
        t_depth_ppm = result.get("transit_depth_ppm")
        t_dur = result.get("transit_duration_hours")
        t_sde = result.get("transit_sde")
        n_tr = result.get("n_transits_observed")
        lines.append(f"  {'Transit Period':20s}: {t_period:.4f} days (SDE: {t_sde}, {n_tr} transits)")
        lines.append(f"  {'Transit Depth':20s}: {t_depth_ppm:.0f} ppm")
        lines.append(f"  {'Transit Duration':20s}: {t_dur:.2f} hours")

        # Planet properties
        rp = result.get("planet_radius_Rearth")
        if rp is not None:
            lines.append(f"  {'--- Planet ---':20s}")
            rp_jup = result.get("planet_radius_Rjup")
            size_class = result.get("planet_size_class", "?")
            lines.append(f"  {'Planet Radius':20s}: {rp:.2f} R_earth ({rp_jup:.3f} R_jup) [{size_class}]")

            a_AU = result.get("orbital_semi_major_axis_AU")
            if a_AU is not None:
                lines.append(f"  {'Orbital Distance':20s}: {a_AU:.5f} AU")

            T_eq = result.get("equilibrium_temp_K")
            if T_eq is not None:
                lines.append(f"  {'Equilibrium Temp':20s}: {T_eq:.0f} K")

            insol = result.get("insolation_Searth")
            if insol is not None:
                lines.append(f"  {'Insolation':20s}: {insol:.1f} S_earth")

            in_hz = result.get("in_habitable_zone")
            if in_hz is not None:
                hz_inner = result.get("hz_conservative_inner_AU", "?")
                hz_outer = result.get("hz_conservative_outer_AU", "?")
                lines.append(f"  {'Habitable Zone':20s}: {'YES' if in_hz else 'no'} (HZ: {hz_inner}-{hz_outer} AU)")

            pflag = result.get("planet_flag", "ok")
            if pflag != "ok":
                lines.append(f"  {'Planet WARNING':20s}: {pflag}")
        else:
            pflag = result.get("planet_flag", "?")
            lines.append(f"  {'Planet Properties':20s}: incomplete ({pflag})")

        # Show alternative transit candidates
        candidates = result.get("transit_candidates", [])
        if len(candidates) > 1:
            lines.append(f"  {'--- Candidates ---':20s}")
            for cand in candidates:
                rank = cand.get("rank", "?")
                cp = cand.get("transit_period_days", 0)
                cs = cand.get("transit_sde", 0)
                cd = cand.get("transit_depth_ppm", 0)
                ca = cand.get("orbital_semi_major_axis_AU")
                ci = cand.get("insolation_Searth")
                ct = cand.get("equilibrium_temp_K")
                chz = cand.get("in_habitable_zone")
                cflag = cand.get("transit_flag", "ok")
                flag_str = "" if cflag == "ok" else f" [{cflag}]"
                a_str = f"{ca:.4f} AU" if ca is not None else "?"
                i_str = f"{ci:.1f} S_e" if ci is not None else "?"
                t_str = f"{ct:.0f}K" if ct is not None else "?"
                hz_str = " HZ" if chz else ""
                marker = " <-- best" if rank == 1 else ""
                lines.append(f"    #{rank}: P={cp:.4f}d  SDE={cs:.1f}  "
                              f"depth={cd:.0f}ppm  a={a_str}  T={t_str}  "
                              f"S={i_str}{hz_str}{flag_str}{marker}")
    elif result.get("transit_detected") is False:
        tflag = result.get("transit_flag", "?")
        lines.append(f"  {'Transit':20s}: not detected ({tflag})")

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


def run(predefined_names=None, simbad_names=None, include_lightcurve=False,
        include_transit=False):
    """Run pipeline on predefined and/or SIMBAD-resolved stars."""
    print("=" * 60)
    print("Stellar Property Pipeline")
    if include_lightcurve:
        print("  (with light curve analysis)")
    if include_transit:
        print("  (with transit detection)")
    print("=" * 60)

    results = {}

    # Predefined stars
    if predefined_names is None and simbad_names is None:
        predefined_names = list(STARS.keys())

    # Map predefined keys to MAST-searchable names for light curve lookup
    PREDEFINED_LC_NAMES = {
        "proxima_cen": "Proxima Cen",
        "sirius_a": "Sirius",
        "delta_cep": "Delta Cep",
        "alpha_cen_a": "Alpha Cen A",
        "barnards_star": "Barnard's Star",
        # sun has no MAST light curve
    }

    for name in (predefined_names or []):
        star_data = STARS.get(name)
        if star_data is None:
            logger.warning("Unknown predefined star: %s (available: %s)", name, ", ".join(STARS.keys()))
            continue
        print(f"\n--- {name} ---")
        lc_target = PREDEFINED_LC_NAMES.get(name) if (include_lightcurve or include_transit) else None
        result = process_star(star_data, include_lightcurve=include_lightcurve,
                              include_transit=include_transit, lc_target=lc_target)
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
        # Use the SIMBAD name directly as the MAST search target
        lc_target = simbad_name if (include_lightcurve or include_transit) else None
        result = process_star(star_data, include_lightcurve=include_lightcurve,
                              include_transit=include_transit, lc_target=lc_target)
        results[display_name] = result
        print(format_result(display_name, result))

    # Save full JSON output -- filename based on star names
    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(exist_ok=True)
    import re
    star_tag = "_".join(results.keys())
    star_tag = re.sub(r"[^\w]+", "_", star_tag).strip("_").lower()
    if not star_tag:
        star_tag = "results"
    output_path = output_dir / f"results_{star_tag}.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Full results saved to %s", output_path)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stellar Property Pipeline")
    parser.add_argument("predefined", nargs="*", help="Predefined star keys (e.g. sun, proxima_cen)")
    parser.add_argument("--name", nargs="+", dest="simbad_names",
                        help="SIMBAD names to resolve (e.g. Vega, 'Alp Lyr', 'HD 172167')")
    parser.add_argument("--lightcurve", action="store_true",
                        help="Enable Module 4: download and analyze MAST light curves")
    parser.add_argument("--transit", action="store_true",
                        help="Enable Module 5: BLS transit detection and planet characterization")
    args = parser.parse_args()

    predefined = args.predefined if args.predefined else None
    run(predefined_names=predefined, simbad_names=args.simbad_names,
        include_lightcurve=args.lightcurve, include_transit=args.transit)
