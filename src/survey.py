"""Survey driver: run the validated RV chain over a target list.

Phase 8.3. Turns the single-target deep-dive chain (filter -> nightly
bin -> Keplerian fit -> activity decorrelation -> GP -> optional MCMC)
into a proximity-ordered survey:

- Iterates targets in distance order (nearest first).
- One scorecard per target: data quality, recovered planets vs
  literature reference, vetoed candidate list, analytic detection
  limits. A target failure never kills the survey.
- Targets with zero confirmed planets (e.g. Tau Ceti) get candidates
  and detection limits from the filtered/binned data directly; no
  Keplerian fit is forced.
- Aggregate outputs: output/survey/survey_results.json and
  survey_summary.md.

Usage:
    from src.survey import run_survey
    result = run_survey(["HD 20794", "HD 115617", "HD 10700"])

    python -m src.survey --targets "HD 20794" "HD 115617" --mcmc
"""

import argparse
import json
import logging
import os
from datetime import datetime

import numpy as np

logger = logging.getLogger(__name__)

# Reference periods for candidate vetoes (days)
SIDEREAL_DAY = 0.99727
ANNUAL_DAYS = 365.25

DEFAULT_SURVEY_TARGETS = ["HD 10700", "HD 20794", "HD 115617"]

OUTPUT_DIR = os.path.join("output", "survey")


# ============================================================
# Candidate vetting (basic set; Phase 8.5 extends this)
# ============================================================

def vet_candidate_period(period, known_periods=(), rotation_days=None,
                         activity_cycle_days=None, baseline_days=None,
                         rel_tol=0.05):
    """Check one periodogram peak against known non-planet explanations.

    Parameters
    ----------
    period : float
        Peak period in days.
    known_periods : iterable of float
        Confirmed planet periods (a match is not a new candidate).
    rotation_days : float, optional
        Stellar rotation period.
    activity_cycle_days : float, optional
        Magnetic activity cycle period.
    baseline_days : float, optional
        Observing baseline; periods without ~2 full cycles are
        unconstrained.
    rel_tol : float
        Relative tolerance for period matching (default 5%).

    Returns
    -------
    (bool, str)
        (True, "pass") if the peak survives all vetoes, else
        (False, reason).
    """
    period = float(period)
    if period <= 0:
        return False, "nonphysical period"

    def close(p_ref, tol=rel_tol):
        return abs(period - p_ref) <= tol * p_ref

    # Known planet: recovered signal, not a new candidate
    for p_known in known_periods:
        if close(float(p_known)):
            return False, f"known planet ({p_known:.2f} d)"

    # Daily sampling aliases of known planets: f_alias = |1/1d -/+ 1/P|
    # (e.g. 61 Vir b at 4.215 d aliases to 1.31 d)
    for p_known in known_periods:
        p_known = float(p_known)
        if p_known <= 1.0:
            continue
        for f_alias in (1.0 - 1.0 / p_known, 1.0 + 1.0 / p_known):
            if f_alias > 0 and close(1.0 / f_alias):
                return False, (f"daily alias of known planet "
                               f"({p_known:.2f} d)")

    # 1-day alias family (sampling: one observation per night)
    for p_alias in (SIDEREAL_DAY, 1.0, 0.5, 2.0):
        if close(p_alias):
            return False, f"1-day alias family ({p_alias:.3f} d)"

    # Annual sampling alias and first harmonic
    for p_annual in (ANNUAL_DAYS, ANNUAL_DAYS / 2.0):
        if abs(period - p_annual) <= 0.10 * p_annual:
            return False, f"annual alias ({p_annual:.1f} d)"

    # Stellar rotation and harmonics (P, P/2, 2P)
    if rotation_days:
        for p_rot in (rotation_days, rotation_days / 2.0,
                      rotation_days * 2.0):
            if abs(period - p_rot) <= 0.15 * p_rot:
                return False, f"rotation harmonic ({p_rot:.1f} d)"

    # Magnetic cycle, its annual beat (the HD 20794 414 d lesson),
    # and anything within 30% of the cycle itself
    if activity_cycle_days:
        cycle = float(activity_cycle_days)
        if abs(period - cycle) <= 0.30 * cycle:
            return False, f"magnetic cycle ({cycle:.0f} d)"
        if cycle > ANNUAL_DAYS:
            beat = 1.0 / (1.0 / ANNUAL_DAYS - 1.0 / cycle)
            if abs(period - beat) <= 0.10 * beat:
                return False, f"cycle-annual beat ({beat:.0f} d)"

    # Coverage: need at least ~2 full cycles to constrain a period
    if baseline_days and period > baseline_days / 2.0:
        return False, f"insufficient coverage (baseline {baseline_days:.0f} d)"

    return True, "pass"


def vet_candidates(peaks, known_periods=(), rotation_days=None,
                   activity_cycle_days=None, baseline_days=None):
    """Apply vet_candidate_period to a periodogram peak list.

    Parameters
    ----------
    peaks : list of dict
        Each with period_days and power (rv_periodogram peak format).

    Returns
    -------
    (list, list)
        (surviving candidates, vetoed peaks with reasons); both lists
        of dicts with period_days, power, and veto_reason on vetoed.
    """
    candidates, vetoed = [], []
    for peak in peaks or []:
        p = peak.get("period_days")
        if p is None:
            continue
        ok, reason = vet_candidate_period(
            p, known_periods=known_periods, rotation_days=rotation_days,
            activity_cycle_days=activity_cycle_days,
            baseline_days=baseline_days,
        )
        entry = {"period_days": float(p), "power": peak.get("power")}
        if ok:
            candidates.append(entry)
        else:
            entry["veto_reason"] = reason
            vetoed.append(entry)
    return candidates, vetoed


# ============================================================
# Scorecard construction
# ============================================================

def _match_recovered_planet(ref_period, fitted_planets, rel_tol=0.10):
    """Find the fitted planet whose period is nearest a reference period."""
    best = None
    for planet in fitted_planets or []:
        p_fit = planet.get("period")
        if p_fit is None:
            continue
        if abs(p_fit - ref_period) <= rel_tol * ref_period:
            if best is None or (abs(p_fit - ref_period)
                                < abs(best.get("period") - ref_period)):
                best = planet
    return best


def _cross_match_candidates(candidates, unconfirmed_planets, rel_tol=0.05):
    """Annotate surviving candidates that match literature candidates."""
    for cand in candidates:
        for planet in unconfirmed_planets:
            p_ref = planet.get("period_days")
            if p_ref and abs(cand["period_days"] - p_ref) <= rel_tol * p_ref:
                cand["matches_literature_candidate"] = planet.get("pl_name")
    return candidates


def _detection_limit_summary(time, residual_rms, baseline_days):
    """Analytic K detection limit from final residual scatter.

    Uses rv_detection_limit with the residual RMS as the effective
    per-point noise (jitter-dominated regime). Analytic placeholder;
    injection-recovery limits are a later refinement.
    """
    from src.rv_data import rv_detection_limit

    time = np.asarray(time, dtype=float)
    if len(time) < 10 or not residual_rms:
        return None
    limits = rv_detection_limit(
        time, np.full(len(time), float(residual_rms)))
    report_periods = [10.0, 30.0, 100.0, 300.0, 1000.0]
    k_at = {}
    for p in report_periods:
        if p > baseline_days / 1.5:
            k_at[f"{p:.0f}_d"] = None
            continue
        idx = int(np.argmin(np.abs(limits["periods"] - p)))
        k_min = limits["k_min_ms"][idx]
        k_at[f"{p:.0f}_d"] = round(float(k_min), 3) if np.isfinite(k_min) else None
    return {
        "method": "analytic (5-sigma, residual RMS as noise)",
        "sigma_eff_ms": round(float(residual_rms), 3),
        "n_obs": int(len(time)),
        "k_min_ms_at_period": k_at,
    }


def build_scorecard(target, planets_ref, config, rv_result,
                    use_gp=True, mcmc=False):
    """Build one survey scorecard from the RV chain result.

    Parameters
    ----------
    target : dict
        Target catalog entry (name, hd, distance_pc, ...).
    planets_ref : list of dict or None
        Literature reference planets (known_planets.json schema).
    config : dict
        Per-target config from known_planets.json.
    rv_result : dict
        Output of deep_dive._run_rv_analysis: rv_data summary,
        residual_analysis, rv_raw.

    Returns
    -------
    dict
        Scorecard with sections: target, data, planets, candidates,
        detection_limit, status.
    """
    planets_ref = planets_ref or []
    confirmed_ref = [p for p in planets_ref if p.get("confirmed")]
    unconfirmed_ref = [p for p in planets_ref
                       if p.get("confirmed") is False]

    card = {
        "target": {
            "name": target.get("name"),
            "hd": target.get("hd"),
            "distance_pc": target.get("distance_pc"),
            "spectral_type": target.get("spectral_type"),
        },
        "chain": {
            "use_gp": bool(use_gp),
            "mcmc": bool(mcmc),
        },
        "status": "ok",
    }

    # --- Data section ---
    rv_summary = rv_result.get("rv_data") or {}
    card["data"] = {
        "status": rv_summary.get("status"),
        "n_measurements_raw": rv_summary.get("n_measurements_raw"),
        "n_measurements_filtered": rv_summary.get("n_measurements_filtered"),
        "n_measurements_binned": rv_summary.get("n_measurements"),
        "time_baseline_days": rv_summary.get("time_baseline_days"),
        "instruments": rv_summary.get("instruments"),
        "excluded_instruments": rv_summary.get("excluded_instruments"),
    }
    if rv_summary.get("status") != "ok":
        card["status"] = "no_rv_data"
        card["planets"] = []
        card["candidates"] = {"surviving": [], "vetoed": []}
        card["detection_limit"] = None
        return card

    baseline = rv_summary.get("time_baseline_days") or 0.0
    rotation = config.get("rotation_days")
    cycle = config.get("activity_cycle_days")

    residual = rv_result.get("residual_analysis")
    rv_raw = rv_result.get("rv_raw") or {}

    # --- Planets section: recovered vs reference ---
    planets_card = []
    fitted = (residual or {}).get("keplerian_planets", [])
    known_periods = [p["period_days"] for p in confirmed_ref
                     if p.get("period_days")]
    for ref in confirmed_ref:
        p_ref = ref.get("period_days")
        entry = {
            "pl_name": ref.get("pl_name"),
            "reference": ref.get("reference"),
            "period_days_ref": p_ref,
            "k_ms_ref": ref.get("k_ms"),
            "recovered": False,
        }
        match = _match_recovered_planet(p_ref, fitted) if p_ref else None
        if match:
            k_fit = match.get("k")
            k_err = match.get("k_err")
            entry.update({
                "recovered": True,
                "period_days_fit": match.get("period"),
                "period_err": match.get("period_err"),
                "k_ms_fit": abs(k_fit) if k_fit is not None else None,
                "k_err": k_err,
                "e_fit": match.get("e"),
            })
            if k_fit is not None and ref.get("k_ms"):
                entry["k_ratio_fit_over_ref"] = round(
                    abs(k_fit) / ref["k_ms"], 3)
                if k_err:
                    entry["k_n_sigma_vs_ref"] = round(
                        abs(abs(k_fit) - ref["k_ms"]) / k_err, 2)
        planets_card.append(entry)
    for ref in unconfirmed_ref:
        planets_card.append({
            "pl_name": ref.get("pl_name"),
            "reference": ref.get("reference"),
            "period_days_ref": ref.get("period_days"),
            "k_ms_ref": ref.get("k_ms"),
            "recovered": None,
            "note": "unconfirmed in literature; not fitted",
        })
    card["planets"] = planets_card
    card["n_planets_confirmed_ref"] = len(confirmed_ref)
    card["n_planets_recovered"] = sum(
        1 for p in planets_card if p.get("recovered") is True)

    # --- Fit quality ---
    if residual:
        card["fit"] = {
            "method": residual.get("method"),
            "rms_before_ms": residual.get("keplerian_rms_before_ms"),
            "rms_after_ms": residual.get("keplerian_rms_after_ms"),
            "decorrelation": bool(residual.get("activity_decorrelation")),
            "gp": bool(residual.get("gp_activity_model")),
        }
    else:
        card["fit"] = None

    # --- Candidates: residual peaks (fit) or raw peaks (no fit) ---
    if residual:
        peaks = residual.get("residual_peaks", [])
        peak_source = "keplerian_residuals"
        residual_rms = residual.get("keplerian_rms_after_ms")
    else:
        peaks = (rv_summary.get("periodogram") or {}).get("peaks", [])
        peak_source = "binned_rv_no_fit"
        # No model subtracted: use centered scatter as the noise floor
        rv = np.asarray(rv_raw.get("rv", []), dtype=float)
        inst = np.asarray(rv_raw.get("instrument", []))
        if len(rv) and len(inst) == len(rv):
            centered = np.array([
                v - np.median(rv[inst == i]) for v, i in zip(rv, inst)])
            residual_rms = float(np.std(centered))
        elif len(rv):
            residual_rms = float(np.std(rv - np.median(rv)))
        else:
            residual_rms = None

    candidates, vetoed = vet_candidates(
        peaks, known_periods=known_periods, rotation_days=rotation,
        activity_cycle_days=cycle, baseline_days=baseline,
    )
    candidates = _cross_match_candidates(candidates, unconfirmed_ref)
    card["candidates"] = {
        "peak_source": peak_source,
        "surviving": candidates,
        "vetoed": vetoed,
        "note": ("candidates are flags for follow-up, not detections; "
                 "FAP per peak and injection tests pending (8.5)"),
    }

    # --- Detection limit ---
    time = rv_raw.get("time")
    card["detection_limit"] = (
        _detection_limit_summary(time, residual_rms, baseline)
        if time is not None else None)
    card["residual_rms_ms"] = (round(float(residual_rms), 3)
                               if residual_rms else None)

    return card


# ============================================================
# Per-target runner and survey loop
# ============================================================

def survey_target(target_name, use_gp=True, mcmc=False, bin_nightly=True,
                  decorrelate_activity=True):
    """Run the RV chain on one target and return its scorecard.

    Reference-first: literature tables in data/known_planets.json
    drive which planets are fitted; NASA archive is only queried for
    targets without a curated entry.
    """
    from src.targets import get_target
    from src.reference_data import load_known_planets, load_target_config
    from src.deep_dive import _run_rv_analysis, _run_known_planets

    target = get_target(target_name)
    if target is None:
        raise ValueError(f"Target '{target_name}' not found in catalog")

    logger.info("=" * 60)
    logger.info("SURVEY TARGET: %s (%s, %.2f pc)", target["name"],
                target.get("hd"), target.get("distance_pc") or -1)
    logger.info("=" * 60)

    hd = target.get("hd")
    planets_ref = load_known_planets(hd) if hd else None
    config = load_target_config(hd) if hd else {}
    reference_source = "known_planets.json"

    if planets_ref is None:
        # Uncurated target: fall back to NASA archive
        nasa = _run_known_planets(target)
        planets_ref = nasa.get("planets") or []
        reference_source = f"nasa_archive ({nasa.get('status')})"

    known_periods = [p["period_days"] for p in planets_ref
                     if p.get("period_days")
                     and p.get("confirmed") is not False]
    logger.info("Reference: %d planets (%d confirmed) from %s",
                len(planets_ref), len(known_periods), reference_source)

    rv_result = _run_rv_analysis(
        target, known_periods, bin_nightly=bin_nightly,
        decorrelate_activity=decorrelate_activity,
        use_gp=use_gp, mcmc_final=mcmc,
    )

    card = build_scorecard(target, planets_ref, config, rv_result,
                           use_gp=use_gp, mcmc=mcmc)
    card["reference_source"] = reference_source
    return card


def run_survey(target_names=None, use_gp=True, mcmc=False,
               bin_nightly=True, decorrelate_activity=True,
               out_dir=OUTPUT_DIR, save=True):
    """Run the RV survey over a target list, nearest star first.

    Parameters
    ----------
    target_names : list of str, optional
        Catalog names or HD identifiers. Default: the Phase 8.3
        anchor set (Tau Ceti, HD 20794, 61 Vir).
    use_gp : bool
        Final GP activity-model fit per target (default True).
    mcmc : bool
        MCMC on the final fit for K error bars. Slow (~35 min per
        target with GP); scorecard-grade runs should enable it.
    save : bool
        Write survey_results.json + survey_summary.md to out_dir.

    Returns
    -------
    dict
        {"targets": [scorecards], "n_ok", "n_failed", ...}
    """
    from src.targets import get_target

    if target_names is None:
        target_names = list(DEFAULT_SURVEY_TARGETS)

    # Order by distance, nearest first (survey principle)
    def dist(name):
        t = get_target(name)
        return t.get("distance_pc") or 1e9 if t else 1e9
    ordered = sorted(target_names, key=dist)

    cards = []
    for name in ordered:
        try:
            card = survey_target(
                name, use_gp=use_gp, mcmc=mcmc, bin_nightly=bin_nightly,
                decorrelate_activity=decorrelate_activity,
            )
        except Exception as e:
            logger.error("Survey target '%s' failed: %s", name, e)
            card = {
                "target": {"name": name},
                "status": "error",
                "error": str(e),
            }
        cards.append(card)

    result = {
        "generated_utc": datetime.utcnow().isoformat(),
        "survey_version": "phase8.3",
        "chain": {"use_gp": use_gp, "mcmc": mcmc,
                  "bin_nightly": bin_nightly,
                  "decorrelate_activity": decorrelate_activity},
        "targets": cards,
        "n_targets": len(cards),
        "n_ok": sum(1 for c in cards if c.get("status") == "ok"),
        "n_failed": sum(1 for c in cards if c.get("status") == "error"),
    }

    if save:
        result["saved_files"] = _save_survey_outputs(result, out_dir)

    logger.info("SURVEY COMPLETE: %d targets, %d ok, %d failed",
                result["n_targets"], result["n_ok"], result["n_failed"])
    return result


# ============================================================
# Output
# ============================================================

def _save_survey_outputs(result, out_dir):
    """Write survey_results.json and survey_summary.md."""
    os.makedirs(out_dir, exist_ok=True)

    json_path = os.path.join(out_dir, "survey_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=_json_default)

    md_path = os.path.join(out_dir, "survey_summary.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(_render_summary_md(result))

    logger.info("Survey outputs saved: %s, %s", json_path, md_path)
    return [json_path, md_path]


def _json_default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


def _render_summary_md(result):
    """Markdown summary table across targets."""
    lines = [
        "# RV Survey Summary",
        "",
        f"Generated: {result['generated_utc']} UTC",
        f"Chain: GP={result['chain']['use_gp']}, "
        f"MCMC={result['chain']['mcmc']}, "
        f"binned={result['chain']['bin_nightly']}, "
        f"decorrelated={result['chain']['decorrelate_activity']}",
        "",
        "| Target | Dist (pc) | N binned | Baseline (d) | "
        "RMS (m/s) | Planets (rec/ref) | Candidates | K_lim@100d (m/s) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for card in result["targets"]:
        t = card.get("target", {})
        if card.get("status") == "error":
            lines.append(f"| {t.get('name')} | - | - | - | - | - | - | "
                         f"ERROR: {card.get('error', '')[:60]} |")
            continue
        data = card.get("data", {})
        limit = card.get("detection_limit") or {}
        k100 = (limit.get("k_min_ms_at_period") or {}).get("100_d")
        n_cand = len((card.get("candidates") or {}).get("surviving", []))
        baseline = data.get("time_baseline_days")
        rms = card.get("residual_rms_ms")
        row = "| {name} ({hd}) | {dist} | {n} | {base} | {rms} | {rec}/{ref} | {cand} | {k100} |".format(
            name=t.get("name"), hd=t.get("hd"),
            dist=t.get("distance_pc"),
            n=data.get("n_measurements_binned"),
            base=f"{baseline:.0f}" if baseline is not None else "-",
            rms=f"{rms:.2f}" if rms is not None else "-",
            rec=card.get("n_planets_recovered", 0),
            ref=card.get("n_planets_confirmed_ref", 0),
            cand=n_cand,
            k100=k100 if k100 is not None else "-",
        )
        lines.append(row)

    lines += ["", "## Per-target notes", ""]
    for card in result["targets"]:
        t = card.get("target", {})
        lines.append(f"### {t.get('name')}")
        if card.get("status") == "error":
            lines.append(f"- ERROR: {card.get('error')}")
            lines.append("")
            continue
        for planet in card.get("planets", []):
            if planet.get("recovered") is True:
                k_fit = planet.get("k_ms_fit")
                k_err = planet.get("k_err")
                k_txt = f"{k_fit:.3f}" if k_fit is not None else "?"
                err_txt = f" +/- {k_err:.3f}" if k_err else ""
                lines.append(
                    f"- {planet['pl_name']}: recovered, "
                    f"K={k_txt}{err_txt} m/s "
                    f"(ref {planet.get('k_ms_ref')})")
            elif planet.get("recovered") is False:
                lines.append(f"- {planet['pl_name']}: NOT recovered")
            else:
                lines.append(f"- {planet['pl_name']}: unconfirmed, "
                             f"not fitted")
        cands = (card.get("candidates") or {}).get("surviving", [])
        for cand in cands:
            match = cand.get("matches_literature_candidate")
            match_txt = f" (matches {match})" if match else ""
            power = cand.get("power")
            power_txt = f"{power:.3f}" if power is not None else "?"
            lines.append(
                f"- Candidate peak: {cand['period_days']:.1f} d, "
                f"power {power_txt}{match_txt}")
        if not cands:
            lines.append("- No candidate peaks survive vetoes")
        lines.append("")
    return "\n".join(lines)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run the RV survey over catalog targets")
    parser.add_argument("--targets", nargs="+",
                        default=None, help="Target names or HD ids")
    parser.add_argument("--mcmc", action="store_true",
                        help="MCMC error bars on final fits (slow)")
    parser.add_argument("--no-gp", action="store_true",
                        help="Skip the GP activity-model stage")
    parser.add_argument("--out-dir", default=OUTPUT_DIR)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s %(message)s")

    run_survey(target_names=args.targets, use_gp=not args.no_gp,
               mcmc=args.mcmc, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
