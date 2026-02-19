"""Batch process TESS confirmed planets and compare with catalog values.

Reads ref_planets.csv, runs the pipeline on each target in parallel via
subprocess calls to run_stars.py, then compares pipeline output against
the reference catalog (period, planet radius, stellar properties).

Usage:
    .\\venv\\Scripts\\python run_toi_batch.py
    .\\venv\\Scripts\\python run_toi_batch.py --workers 4
    .\\venv\\Scripts\\python run_toi_batch.py --no-reuse      # re-run all targets
    .\\venv\\Scripts\\python run_toi_batch.py --compare-only   # skip processing
"""

import argparse
import csv
import json
import logging
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

LOG_DIR = Path(__file__).resolve().parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

VENV_PYTHON = str(Path(__file__).resolve().parent / "venv" / "Scripts" / "python.exe")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "toi_batch.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _float_or_none(val):
    """Parse a float value, returning None for empty/invalid."""
    if val is None or val == "":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _hostname_to_slug(hostname):
    """Convert hostname to output file slug (matching run_stars.py logic)."""
    return re.sub(r"[^\w]+", "_", hostname).strip("_").lower()


def _pct_err(computed, reference):
    """Compute percent error, or None if either value is missing."""
    if computed is None or reference is None or reference == 0:
        return None
    return abs(computed - reference) / abs(reference) * 100


# ---------------------------------------------------------------------------
# CSV parsing
# ---------------------------------------------------------------------------

def load_reference_catalog(csv_path):
    """Parse ref_planets.csv into a list of reference dicts."""
    targets = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            target = {
                "pl_name": row["pl_name"],
                "hostname": row["hostname"],
                "ref_period_days": _float_or_none(row.get("pl_orbper")),
                "ref_planet_radius_Rearth": _float_or_none(row.get("pl_rade")),
                "ref_planet_radius_Rjup": _float_or_none(row.get("pl_radj")),
                "ref_semi_major_axis_AU": _float_or_none(row.get("pl_orbsmax")),
                "ref_equilibrium_temp_K": _float_or_none(row.get("pl_eqt")),
                "ref_insolation_Searth": _float_or_none(row.get("pl_insol")),
                "ref_st_teff_K": _float_or_none(row.get("st_teff")),
                "ref_st_radius_Rsun": _float_or_none(row.get("st_rad")),
                "ref_st_mass_Msun": _float_or_none(row.get("st_mass")),
                "ref_distance_pc": _float_or_none(row.get("sy_dist")),
            }
            # Estimate reference transit depth from (Rp/Rs)^2
            if target["ref_planet_radius_Rearth"] and target["ref_st_radius_Rsun"]:
                rp_rsun = target["ref_planet_radius_Rearth"] / 109.076  # R_earth -> R_sun
                ratio = rp_rsun / target["ref_st_radius_Rsun"]
                target["ref_depth_ppm"] = round(ratio ** 2 * 1e6, 1)
            else:
                target["ref_depth_ppm"] = None
            targets.append(target)
    return targets


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def run_single_target(hostname, timeout_s=600):
    """Run pipeline for a single target via subprocess.

    Returns (hostname, returncode, duration_s).
    """
    slug = _hostname_to_slug(hostname)
    log_path = LOG_DIR / f"toi_batch_{slug}.log"

    cmd = [VENV_PYTHON, "run_stars.py", "--name", hostname, "--transit"]

    logger.info("Starting: %s", hostname)
    t0 = time.time()

    try:
        with open(log_path, "w") as log_file:
            proc = subprocess.run(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                timeout=timeout_s,
                cwd=str(Path(__file__).resolve().parent),
            )
        duration = time.time() - t0
        logger.info("Finished: %s (rc=%d, %.1fs)", hostname, proc.returncode, duration)
        return hostname, proc.returncode, duration
    except subprocess.TimeoutExpired:
        duration = time.time() - t0
        logger.warning("Timeout: %s after %.0fs", hostname, duration)
        return hostname, -1, duration
    except Exception as e:
        duration = time.time() - t0
        logger.error("Error processing %s: %s", hostname, e)
        return hostname, -2, duration


def load_pipeline_result(hostname):
    """Load pipeline output JSON for a hostname. Returns dict or None."""
    slug = _hostname_to_slug(hostname)
    output_path = OUTPUT_DIR / f"results_{slug}.json"
    if not output_path.exists():
        return None
    with open(output_path, "r") as f:
        data = json.load(f)
    # run_stars.py stores as {display_name: result_dict}
    if isinstance(data, dict) and hostname in data:
        return data[hostname]
    # Fallback: single-entry dict with slightly different key
    if isinstance(data, dict) and len(data) == 1:
        return next(iter(data.values()))
    return data


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def find_best_period_match(candidates, ref_period, tolerance=0.05):
    """Find the candidate whose period best matches the reference.

    Checks direct period and common aliases (2x, 0.5x).
    Returns (matching_candidate, relative_error) or (None, best_error).
    """
    if not candidates or ref_period is None:
        return None, None

    best_match = None
    best_error = float("inf")

    for cand in candidates:
        p = cand.get("transit_period_days")
        if p is None:
            continue
        # Check direct period, double, and half (common BLS aliases)
        for test_p in [p, p * 2, p / 2]:
            rel_err = abs(test_p - ref_period) / ref_period
            if rel_err < best_error:
                best_error = rel_err
                best_match = cand

    if best_error <= tolerance:
        return best_match, best_error
    return None, best_error


def compare_results(targets, results):
    """Compare pipeline results against reference catalog. Returns list of comparison dicts."""
    comparisons = []

    for target in targets:
        hostname = target["hostname"]
        result = results.get(hostname)

        comp = {
            "hostname": hostname,
            "pl_name": target["pl_name"],
            "status": "no_result",
        }

        if result is None:
            comparisons.append(comp)
            continue

        comp["status"] = "ok"

        # -- Stellar properties --
        comp["distance_pc"] = result.get("distance_pc")
        comp["ref_distance_pc"] = target["ref_distance_pc"]
        comp["distance_err_pct"] = _pct_err(result.get("distance_pc"), target["ref_distance_pc"])

        comp["teff_K"] = result.get("teff_K")
        comp["ref_teff_K"] = target["ref_st_teff_K"]
        comp["teff_err_pct"] = _pct_err(result.get("teff_K"), target["ref_st_teff_K"])

        comp["radius_Rsun"] = result.get("radius_Rsun")
        comp["ref_radius_Rsun"] = target["ref_st_radius_Rsun"]
        comp["radius_err_pct"] = _pct_err(result.get("radius_Rsun"), target["ref_st_radius_Rsun"])

        comp["mass_Msun"] = result.get("mass_Msun")
        comp["ref_mass_Msun"] = target["ref_st_mass_Msun"]
        comp["mass_err_pct"] = _pct_err(result.get("mass_Msun"), target["ref_st_mass_Msun"])

        # -- Transit detection --
        comp["transit_detected"] = result.get("transit_detected")
        candidates = result.get("transit_candidates", [])
        comp["n_candidates"] = len(candidates)

        # Check quiet-star skip
        if result.get("planet_flag") == "too_variable":
            comp["status"] = "too_variable"
        elif result.get("transit_detected") is None and result.get("planet_flag") == "no_lightcurve":
            comp["status"] = "no_lightcurve"

        # Find candidate matching reference period
        ref_period = target["ref_period_days"]
        comp["ref_period_days"] = ref_period

        match, period_err = find_best_period_match(candidates, ref_period)
        if match is not None:
            comp["period_matched"] = True
            comp["matched_period_days"] = match.get("transit_period_days")
            comp["period_err_pct"] = period_err * 100 if period_err is not None else None
            comp["matched_rank"] = match.get("rank")
            comp["matched_sde"] = match.get("transit_sde")
            comp["matched_depth_ppm"] = match.get("transit_depth_ppm")
            comp["ref_depth_ppm"] = target["ref_depth_ppm"]
            comp["depth_err_pct"] = _pct_err(match.get("transit_depth_ppm"), target["ref_depth_ppm"])

            comp["matched_planet_radius_Rearth"] = match.get("planet_radius_Rearth")
            comp["ref_planet_radius_Rearth"] = target["ref_planet_radius_Rearth"]
            comp["planet_radius_err_pct"] = _pct_err(
                match.get("planet_radius_Rearth"), target["ref_planet_radius_Rearth"]
            )

            comp["matched_semi_major_axis_AU"] = match.get("orbital_semi_major_axis_AU")
            comp["ref_semi_major_axis_AU"] = target["ref_semi_major_axis_AU"]
            comp["semi_major_axis_err_pct"] = _pct_err(
                match.get("orbital_semi_major_axis_AU"), target["ref_semi_major_axis_AU"]
            )

            comp["matched_eq_temp_K"] = match.get("equilibrium_temp_K")
            comp["ref_eq_temp_K"] = target["ref_equilibrium_temp_K"]
            comp["eq_temp_err_pct"] = _pct_err(
                match.get("equilibrium_temp_K"), target["ref_equilibrium_temp_K"]
            )

            comp["matched_flag"] = match.get("transit_flag")
            comp["matched_planet_flag"] = match.get("planet_flag")
        else:
            comp["period_matched"] = False
            comp["closest_period_err_pct"] = period_err * 100 if period_err is not None else None
            comp["ref_planet_radius_Rearth"] = target["ref_planet_radius_Rearth"]
            comp["ref_depth_ppm"] = target["ref_depth_ppm"]
            if candidates:
                comp["best_period_days"] = candidates[0].get("transit_period_days")
                comp["best_sde"] = candidates[0].get("transit_sde")

        comparisons.append(comp)

    return comparisons


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_summary(comparisons):
    """Print formatted summary tables of the comparison results."""

    print("\n" + "=" * 130)
    print("TESS CONFIRMED PLANET PIPELINE COMPARISON")
    print("=" * 130)

    # -- Stellar properties table --
    print("\n--- Stellar Properties ---")
    hdr = (f"{'Target':<18} {'Dist(pc)':>9} {'Ref':>9} {'Err%':>6}  "
           f"{'Teff(K)':>7} {'Ref':>7} {'Err%':>6}  "
           f"{'R*(Rs)':>7} {'Ref':>7} {'Err%':>6}  "
           f"{'M*(Ms)':>7} {'Ref':>7} {'Err%':>6}")
    print(hdr)
    print("-" * 130)

    for c in comparisons:
        if c["status"] == "no_result":
            print(f"{c['hostname']:<18} --- no result ---")
            continue

        def _v(val, fmt):
            return f"{val:{fmt}}" if val is not None else "N/A"

        d = _v(c.get("distance_pc"), ".1f")
        d_r = _v(c.get("ref_distance_pc"), ".1f")
        d_e = _v(c.get("distance_err_pct"), ".1f")
        t = _v(c.get("teff_K"), ".0f")
        t_r = _v(c.get("ref_teff_K"), ".0f")
        t_e = _v(c.get("teff_err_pct"), ".1f")
        r = _v(c.get("radius_Rsun"), ".3f")
        r_r = _v(c.get("ref_radius_Rsun"), ".3f")
        r_e = _v(c.get("radius_err_pct"), ".1f")
        m = _v(c.get("mass_Msun"), ".3f")
        m_r = _v(c.get("ref_mass_Msun"), ".3f")
        m_e = _v(c.get("mass_err_pct"), ".1f")

        print(f"{c['hostname']:<18} {d:>9} {d_r:>9} {d_e:>6}  "
              f"{t:>7} {t_r:>7} {t_e:>6}  "
              f"{r:>7} {r_r:>7} {r_e:>6}  "
              f"{m:>7} {m_r:>7} {m_e:>6}")

    # -- Transit detection table --
    print("\n--- Transit Detection ---")
    hdr2 = (f"{'Target':<18} {'Det':>3} {'#C':>3} {'PMatch':>6} "
            f"{'P_det(d)':>10} {'P_ref(d)':>10} {'P_err%':>7} "
            f"{'Rk':>3} {'SDE':>5} "
            f"{'Rp(Re)':>7} {'Ref':>7} {'Rp_err%':>7} "
            f"{'a(AU)':>7} {'Ref':>7} {'a_err%':>7}")
    print(hdr2)
    print("-" * 130)

    n_ok = 0
    n_detected = 0
    n_matched = 0

    for c in comparisons:
        if c["status"] == "no_result":
            print(f"{c['hostname']:<18} --- no result ---")
            continue

        n_ok += 1
        det = "Y" if c.get("transit_detected") else "N"
        if c.get("status") == "too_variable":
            det = "skp"

        n_cand = c.get("n_candidates", 0)

        if c.get("transit_detected"):
            n_detected += 1

        if c.get("period_matched"):
            n_matched += 1
            p_det = f"{c['matched_period_days']:.4f}"
            p_ref = f"{c['ref_period_days']:.4f}"
            p_err = f"{c.get('period_err_pct', 0):.2f}"
            rank = f"{c['matched_rank']}"
            sde = f"{c['matched_sde']:.1f}"
            rp = f"{c.get('matched_planet_radius_Rearth', 0):.2f}" if c.get("matched_planet_radius_Rearth") else "N/A"
            rp_ref = f"{c.get('ref_planet_radius_Rearth', 0):.2f}" if c.get("ref_planet_radius_Rearth") else "N/A"
            rp_err = f"{c.get('planet_radius_err_pct', 0):.1f}" if c.get("planet_radius_err_pct") is not None else "N/A"
            a_det = f"{c.get('matched_semi_major_axis_AU', 0):.4f}" if c.get("matched_semi_major_axis_AU") else "N/A"
            a_ref = f"{c.get('ref_semi_major_axis_AU', 0):.4f}" if c.get("ref_semi_major_axis_AU") else "N/A"
            a_err = f"{c.get('semi_major_axis_err_pct', 0):.1f}" if c.get("semi_major_axis_err_pct") is not None else "N/A"
            pmatch = "YES"
        else:
            p_det = "N/A"
            if c.get("best_period_days"):
                p_det = f"{c['best_period_days']:.4f}"
            p_ref = f"{c.get('ref_period_days', 0):.4f}" if c.get("ref_period_days") else "N/A"
            p_err = f"{c.get('closest_period_err_pct', 0):.1f}" if c.get("closest_period_err_pct") is not None else "N/A"
            rank = "-"
            sde = "-"
            rp = "-"
            rp_ref = f"{c.get('ref_planet_radius_Rearth', 0):.2f}" if c.get("ref_planet_radius_Rearth") else "N/A"
            rp_err = "-"
            a_det = "-"
            a_ref = f"{c.get('ref_semi_major_axis_AU', 0):.4f}" if c.get("ref_semi_major_axis_AU") else "N/A"
            a_err = "-"
            pmatch = "no"

        print(f"{c['hostname']:<18} {det:>3} {n_cand:>3} {pmatch:>6} "
              f"{p_det:>10} {p_ref:>10} {p_err:>7} "
              f"{rank:>3} {sde:>5} "
              f"{rp:>7} {rp_ref:>7} {rp_err:>7} "
              f"{a_det:>7} {a_ref:>7} {a_err:>7}")

    # -- Summary statistics --
    print("\n--- Summary ---")
    n_total = len(comparisons)
    print(f"Total targets:       {n_total}")
    print(f"Successfully run:    {n_ok}")
    print(f"Transit detected:    {n_detected} / {n_ok}")
    print(f"Period matched:      {n_matched} / {n_ok} (within 5%% tolerance, incl. 2x/0.5x aliases)")

    # Average errors for targets with matched periods
    matched = [c for c in comparisons if c.get("period_matched")]
    if matched:
        avg_p = sum(c.get("period_err_pct", 0) for c in matched) / len(matched)
        rp_errs = [c["planet_radius_err_pct"] for c in matched if c.get("planet_radius_err_pct") is not None]
        a_errs = [c["semi_major_axis_err_pct"] for c in matched if c.get("semi_major_axis_err_pct") is not None]
        d_errs = [c["distance_err_pct"] for c in matched if c.get("distance_err_pct") is not None]
        t_errs = [c["teff_err_pct"] for c in matched if c.get("teff_err_pct") is not None]
        r_errs = [c["radius_err_pct"] for c in matched if c.get("radius_err_pct") is not None]
        m_errs = [c["mass_err_pct"] for c in matched if c.get("mass_err_pct") is not None]

        print(f"\nFor {len(matched)} period-matched targets (avg errors):")
        print(f"  Period:            {avg_p:.3f}%")
        if rp_errs:
            print(f"  Planet radius:     {sum(rp_errs)/len(rp_errs):.1f}%")
        if a_errs:
            print(f"  Semi-major axis:   {sum(a_errs)/len(a_errs):.1f}%")
        if d_errs:
            print(f"  Distance:          {sum(d_errs)/len(d_errs):.1f}%")
        if t_errs:
            print(f"  Teff:              {sum(t_errs)/len(t_errs):.1f}%")
        if r_errs:
            print(f"  Stellar radius:    {sum(r_errs)/len(r_errs):.1f}%")
        if m_errs:
            print(f"  Stellar mass:      {sum(m_errs)/len(m_errs):.1f}%")

    # All-target stellar property averages
    all_ok = [c for c in comparisons if c["status"] not in ("no_result",)]
    if all_ok:
        d_all = [c["distance_err_pct"] for c in all_ok if c.get("distance_err_pct") is not None]
        t_all = [c["teff_err_pct"] for c in all_ok if c.get("teff_err_pct") is not None]
        r_all = [c["radius_err_pct"] for c in all_ok if c.get("radius_err_pct") is not None]
        m_all = [c["mass_err_pct"] for c in all_ok if c.get("mass_err_pct") is not None]

        print(f"\nFor all {len(all_ok)} processed targets (stellar property avg errors):")
        if d_all:
            print(f"  Distance:          {sum(d_all)/len(d_all):.1f}%")
        if t_all:
            print(f"  Teff:              {sum(t_all)/len(t_all):.1f}%")
        if r_all:
            print(f"  Stellar radius:    {sum(r_all)/len(r_all):.1f}%")
        if m_all:
            print(f"  Stellar mass:      {sum(m_all)/len(m_all):.1f}%")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Batch process TESS confirmed planets")
    parser.add_argument("--csv", type=str, default="ref_planets.csv",
                        help="Path to reference catalog CSV (default: ref_planets.csv)")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    parser.add_argument("--timeout", type=int, default=600,
                        help="Timeout per target in seconds (default: 600)")
    parser.add_argument("--no-reuse", action="store_true",
                        help="Re-run all targets even if output already exists")
    parser.add_argument("--compare-only", action="store_true",
                        help="Skip processing, only compare existing results")
    args = parser.parse_args()

    csv_path = Path(__file__).resolve().parent / args.csv
    targets = load_reference_catalog(csv_path)
    logger.info("Loaded %d targets from %s", len(targets), csv_path)

    hostnames = [t["hostname"] for t in targets]

    if not args.compare_only:
        # Determine which targets need processing
        to_process = []
        for hostname in hostnames:
            slug = _hostname_to_slug(hostname)
            output_path = OUTPUT_DIR / f"results_{slug}.json"
            if output_path.exists() and not args.no_reuse:
                logger.info("Reusing existing result for %s (%s)", hostname, output_path.name)
            else:
                to_process.append(hostname)

        if to_process:
            logger.info("Processing %d targets with %d workers: %s",
                        len(to_process), args.workers, ", ".join(to_process))

            t0 = time.time()
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = {
                    executor.submit(run_single_target, h, args.timeout): h
                    for h in to_process
                }
                for future in as_completed(futures):
                    hostname_done, rc, duration = future.result()
                    status = "OK" if rc == 0 else f"FAILED (rc={rc})"
                    logger.info("  %s: %s (%.1fs)", hostname_done, status, duration)

            total_time = time.time() - t0
            logger.info("All targets processed in %.1fs", total_time)
        else:
            logger.info("All targets have existing results. Use --no-reuse to reprocess.")

    # Load all results
    results = {}
    for hostname in hostnames:
        result = load_pipeline_result(hostname)
        if result is not None:
            results[hostname] = result
    logger.info("Loaded results for %d / %d targets", len(results), len(hostnames))

    # Compare and report
    comparisons = compare_results(targets, results)
    print_summary(comparisons)

    # Save full comparison JSON
    comp_path = OUTPUT_DIR / "toi_batch_comparison.json"
    with open(comp_path, "w") as f:
        json.dump(comparisons, f, indent=2)
    logger.info("Comparison saved to %s", comp_path)


if __name__ == "__main__":
    main()
