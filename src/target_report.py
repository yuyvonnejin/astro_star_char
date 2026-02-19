"""Phase 7: Per-target report generator.

Orchestrates all available analysis methods for a single target star:
stellar properties (Modules 1-3), light curve + transit (Modules 4-5),
known planet query, RV data inventory, and proper motion anomaly.

Produces structured report dicts and optional markdown summaries.
"""

import json
import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output" / "target_reports"


def generate_report(target, include_stellar=True, include_lightcurve=True,
                    include_transit=True, include_rv=True, include_pma=True,
                    include_known_planets=True):
    """Generate a comprehensive report for a single target star.

    Parameters
    ----------
    target : dict
        Target dictionary from src.targets (must have name, hd, hip, etc.).
    include_stellar : bool
        Run Modules 1-3 (stellar properties via Gaia DR3).
    include_lightcurve : bool
        Run Module 4 (TESS light curve retrieval + variability analysis).
    include_transit : bool
        Run Module 5 (BLS transit detection).
    include_rv : bool
        Query RV archives and run RV periodogram.
    include_pma : bool
        Compute Gaia-Hipparcos proper motion anomaly.
    include_known_planets : bool
        Query NASA Exoplanet Archive for confirmed planets.

    Returns
    -------
    dict
        Structured report with sections: target, stellar_properties,
        lightcurve, transit, known_planets, rv_data, proper_motion_anomaly,
        data_summary.
    """
    report = {
        "target": _target_summary(target),
        "generated_utc": datetime.utcnow().isoformat(),
        "sections_included": [],
    }

    # 1. Stellar properties via existing pipeline
    if include_stellar:
        report["sections_included"].append("stellar_properties")
        report["stellar_properties"] = _run_stellar_properties(target)

    # 2. Known planets from NASA Exoplanet Archive
    if include_known_planets:
        report["sections_included"].append("known_planets")
        report["known_planets"] = _run_known_planets(target)

    # 3. Light curve analysis
    if include_lightcurve:
        report["sections_included"].append("lightcurve")
        report["lightcurve"] = _run_lightcurve(target)

    # 4. Transit detection (requires light curve)
    if include_transit and include_lightcurve:
        report["sections_included"].append("transit")
        report["transit"] = _run_transit(target, report.get("stellar_properties"),
                                          report.get("lightcurve"))

    # 5. RV data inventory and analysis
    if include_rv:
        report["sections_included"].append("rv_data")
        report["rv_data"] = _run_rv_analysis(target)

    # 6. Proper motion anomaly
    if include_pma:
        report["sections_included"].append("proper_motion_anomaly")
        report["proper_motion_anomaly"] = _run_pma(target)

    # 7. Data availability summary
    report["data_summary"] = _summarize_data(report)

    return report


def format_report_markdown(report):
    """Format a report dict as a markdown string.

    Parameters
    ----------
    report : dict
        Output from generate_report().

    Returns
    -------
    str
        Markdown-formatted report.
    """
    lines = []
    t = report["target"]
    lines.append(f"# Target Report: {t['name']}")
    lines.append(f"**HD**: {t.get('hd', 'N/A')} | **HIP**: {t.get('hip', 'N/A')} | "
                 f"**Spectral Type**: {t.get('spectral_type', 'N/A')}")
    lines.append(f"**Generated**: {report.get('generated_utc', 'N/A')}")
    lines.append("")

    # Stellar properties
    sp = report.get("stellar_properties", {})
    if sp:
        lines.append("## Stellar Properties")
        lines.append("")
        lines.append("| Property | Pipeline Value | Catalog Value |")
        lines.append("|----------|---------------|---------------|")
        for key in ["distance_pc", "teff_K", "radius_Rsun", "mass_Msun", "luminosity_Lsun"]:
            pipe_val = sp.get(key)
            cat_key = {"distance_pc": "distance_pc", "teff_K": "teff_k",
                       "radius_Rsun": None, "mass_Msun": "mass_msun",
                       "luminosity_Lsun": None}.get(key)
            cat_val = t.get(cat_key) if cat_key else "N/A"
            pipe_str = f"{pipe_val:.4f}" if isinstance(pipe_val, float) else str(pipe_val)
            cat_str = f"{cat_val:.4f}" if isinstance(cat_val, float) else str(cat_val)
            lines.append(f"| {key} | {pipe_str} | {cat_str} |")
        lines.append("")

    # Known planets
    kp = report.get("known_planets", {})
    if kp:
        planets = kp.get("planets", [])
        lines.append(f"## Known Planets ({len(planets)} found)")
        lines.append("")
        if planets:
            lines.append("| Planet | Period (d) | Mass (Mearth) | Radius (Rearth) | Method | Year |")
            lines.append("|--------|-----------|---------------|-----------------|--------|------|")
            for p in planets:
                lines.append(
                    f"| {p.get('pl_name', 'N/A')} | "
                    f"{_fmt(p.get('period_days'))} | "
                    f"{_fmt(p.get('mass_earth'))} | "
                    f"{_fmt(p.get('radius_earth'))} | "
                    f"{p.get('discovery_method', 'N/A')} | "
                    f"{p.get('discovery_year', 'N/A')} |"
                )
        else:
            lines.append("No confirmed planets in NASA Exoplanet Archive.")
        lines.append("")

    # Light curve
    lc = report.get("lightcurve", {})
    if lc:
        lines.append("## Light Curve Analysis")
        lines.append("")
        if lc.get("available"):
            lines.append(f"- **Mission**: {lc.get('mission', 'N/A')}")
            lines.append(f"- **Sectors**: {lc.get('n_sectors', 'N/A')}")
            lines.append(f"- **Variability**: {lc.get('variability_class', 'N/A')}")
            lines.append(f"- **Amplitude**: {_fmt(lc.get('amplitude_ppt'))} ppt")
        else:
            lines.append("No TESS light curve data available.")
        lines.append("")

    # Transit
    tr = report.get("transit", {})
    if tr:
        lines.append("## Transit Search")
        lines.append("")
        if tr.get("transit_detected"):
            n_cand = tr.get("n_candidates", 0)
            lines.append(f"- **Transit detected**: Yes ({n_cand} candidate(s))")
            candidates = tr.get("candidates", [])
            if candidates:
                lines.append("")
                lines.append("| # | Period (d) | Depth (ppm) | Rp (Rearth) | Shape |")
                lines.append("|---|-----------|-------------|-------------|-------|")
                for i, c in enumerate(candidates):
                    lines.append(
                        f"| {i+1} | {_fmt(c.get('period_days'))} | "
                        f"{_fmt(c.get('depth_ppm'))} | "
                        f"{_fmt(c.get('planet_radius_Rearth'))} | "
                        f"{c.get('shape', 'N/A')} |"
                    )
        else:
            lines.append("- **Transit detected**: No")
        lines.append("")

    # RV data
    rv = report.get("rv_data", {})
    if rv:
        lines.append("## Radial Velocity Data")
        lines.append("")
        if rv.get("n_measurements", 0) > 0:
            lines.append(f"- **Source**: {rv.get('source', 'N/A')}")
            lines.append(f"- **Measurements**: {rv['n_measurements']}")
            lines.append(f"- **Baseline**: {_fmt(rv.get('time_baseline_days'))} days")
            lines.append(f"- **Instruments**: {', '.join(rv.get('instruments', []))}")
            if rv.get("periodogram"):
                pg = rv["periodogram"]
                lines.append(f"- **Best RV period**: {_fmt(pg.get('best_period'))} days "
                             f"(FAP={_fmt(pg.get('fap'))})")
        elif rv.get("status") == "orbital_solutions_only":
            n_sol = rv.get("n_solutions", 0)
            lines.append(f"- **RV orbital solutions**: {n_sol} (from NASA Exoplanet Archive)")
            lines.append(f"- **Note**: {rv.get('note', 'Raw RV time series not available')}")
        else:
            note = rv.get("note", "No public RV data found.")
            lines.append(note)
        lines.append("")

    # PMa
    pma = report.get("proper_motion_anomaly", {})
    if pma:
        lines.append("## Proper Motion Anomaly")
        lines.append("")
        pma_data = pma.get("pma", {})
        if pma_data:
            lines.append(f"- **PMa total**: {_fmt(pma_data.get('pma_total_mas_yr'))} mas/yr")
            lines.append(f"- **SNR**: {_fmt(pma_data.get('pma_snr'))}")
            lines.append(f"- **Significant**: {pma_data.get('significant', 'N/A')}")
        ruwe = pma.get("ruwe")
        if ruwe is not None:
            lines.append(f"- **RUWE**: {ruwe:.3f}")
        lines.append("")

    # Data summary
    ds = report.get("data_summary", {})
    if ds:
        lines.append("## Data Availability Summary")
        lines.append("")
        methods = ds.get("methods_available", [])
        lines.append(f"- **Methods available**: {', '.join(methods) if methods else 'None'}")
        notes = ds.get("notes", [])
        for note in notes:
            lines.append(f"- {note}")
        lines.append("")

    return "\n".join(lines)


def save_report(report, output_dir=None):
    """Save report as JSON and markdown files.

    Parameters
    ----------
    report : dict
        Output from generate_report().
    output_dir : str or Path, optional
        Output directory. Default: output/target_reports/.

    Returns
    -------
    dict
        Paths to saved files: json_path, md_path.
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    target_name = report["target"]["name"].replace(" ", "_").lower()
    timestamp = datetime.utcnow().strftime("%Y%m%d")

    # Save JSON (convert numpy arrays to lists)
    json_path = output_dir / f"{target_name}_{timestamp}.json"
    json_str = json.dumps(report, indent=2, default=_json_serializer)
    json_path.write_text(json_str, encoding="utf-8")
    logger.info("Report saved: %s", json_path)

    # Save markdown
    md_path = output_dir / f"{target_name}_{timestamp}.md"
    md_str = format_report_markdown(report)
    md_path.write_text(md_str, encoding="utf-8")
    logger.info("Report saved: %s", md_path)

    return {"json_path": str(json_path), "md_path": str(md_path)}


# --- Internal section runners ---

def _target_summary(target):
    """Extract target summary from catalog entry."""
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
        "feh": target.get("feh"),
        "tier": target.get("tier"),
        "known_planets": target.get("known_planets"),
        "planet_status": target.get("planet_status"),
        "notes": target.get("notes"),
    }


def _run_stellar_properties(target):
    """Run existing pipeline Modules 1-3 on the target."""
    from src.data_access import resolve_simbad_name, query_stars_by_id
    from src.pipeline import process_star

    gaia_id = target.get("gaia_dr3_id")
    if gaia_id is None:
        # Try common name first, then HD number
        gaia_id = resolve_simbad_name(target["name"])
        if gaia_id is None and target.get("hd"):
            gaia_id = resolve_simbad_name(target["hd"])

    if gaia_id is None:
        logger.warning("Cannot run stellar properties for %s: no Gaia ID", target["name"])
        return {"status": "no_gaia_id"}

    stars = query_stars_by_id([gaia_id])
    if not stars:
        return {"status": "gaia_query_failed"}

    star_dict = stars[0]
    result = process_star(star_dict, include_lightcurve=False, include_transit=False)
    result["status"] = "ok"
    return result


def _run_known_planets(target):
    """Query NASA Exoplanet Archive for known planets."""
    from src.rv_data import query_known_planets

    # Try HD name first
    planets = query_known_planets(target["hd"])

    if not planets:
        # Try common name
        planets = query_known_planets(target["name"])

    return {
        "n_planets": len(planets),
        "planets": planets,
        "status": "ok" if planets else "none_found",
    }


def _run_lightcurve(target):
    """Run TESS light curve retrieval and variability analysis."""
    from src.lightcurve import search_lightcurve

    # Check availability first
    search_names = []
    if target.get("tic"):
        search_names.append(f"TIC {target['tic']}")
    search_names.append(target["hd"])
    search_names.append(target["name"])

    for search_name in search_names:
        try:
            info = search_lightcurve(search_name, mission="TESS")
            if info is not None:
                # Just report availability; full download deferred to transit analysis
                return {
                    "available": True,
                    "n_sectors": info["n_available"],
                    "mission": info.get("mission"),
                    "author": info.get("author"),
                    "search_name": search_name,
                    "status": "ok",
                    # Variability analysis would need full download
                    "variability_class": None,
                    "amplitude_ppt": None,
                }
        except Exception as e:
            logger.warning("Light curve search failed for '%s': %s", search_name, e)

    return {"available": False, "status": "no_data"}


def _run_transit(target, stellar_props, lc_info):
    """Run transit detection if light curve is available."""
    if not lc_info or not lc_info.get("available"):
        return {"transit_detected": None, "status": "no_lightcurve"}

    # Transit detection requires downloading the actual light curve data.
    # For the report, we just note availability; full transit analysis is
    # done via the main pipeline when running deep-dive.
    return {
        "transit_detected": None,
        "status": "deferred_to_deep_dive",
        "note": "Full transit analysis requires light curve download; use pipeline.process_star()",
    }


def _run_rv_analysis(target):
    """Query RV archives and run periodogram if data available."""
    from src.rv_data import query_dace_rv, query_nasa_rv_data, rv_periodogram

    # Try DACE first for raw RV time series
    rv_data = query_dace_rv(target["hd"])

    if rv_data is not None:
        result = {
            "n_measurements": rv_data["n_measurements"],
            "time_baseline_days": rv_data["time_baseline_days"],
            "instruments": rv_data["instruments"],
            "source": "DACE",
            "status": "ok",
        }

        # Run periodogram if enough data
        if rv_data["n_measurements"] >= 20:
            pg = rv_periodogram(rv_data["time"], rv_data["rv"], rv_data["rv_err"])
            result["periodogram"] = {
                "best_period": pg.get("best_period"),
                "best_power": pg.get("best_power"),
                "fap": pg.get("fap"),
                "peaks": pg.get("peaks", []),
            }
        else:
            result["periodogram"] = None
        return result

    # Fallback: get orbital solution parameters from NASA (not raw time series)
    nasa_rv = query_nasa_rv_data(target["hd"])
    if nasa_rv is not None:
        return {
            "n_measurements": 0,
            "instruments": [],
            "source": "NASA_orbital_solutions",
            "n_solutions": nasa_rv.get("n_solutions", 0),
            "solutions": nasa_rv.get("solutions", []),
            "status": "orbital_solutions_only",
            "note": "Raw RV time series unavailable; DACE unreachable. "
                    "Orbital solutions from NASA Exoplanet Archive shown instead.",
        }

    return {
        "n_measurements": 0,
        "instruments": [],
        "status": "no_data",
        "note": "DACE unreachable and no NASA RV data. "
                "Raw RV time series requires direct ESO archive access.",
    }


def _run_pma(target):
    """Run proper motion anomaly analysis."""
    from src.proper_motion import analyze_pma

    gaia_id = target.get("gaia_dr3_id")
    hip_id = target.get("hip")

    if gaia_id is None:
        # Try to resolve via common name, then HD number
        from src.data_access import resolve_simbad_name
        gaia_id = resolve_simbad_name(target["name"])
        if gaia_id is None and target.get("hd"):
            gaia_id = resolve_simbad_name(target["hd"])

    if gaia_id is None or hip_id is None:
        return {"status": "missing_ids", "gaia_dr3_id": gaia_id, "hip_id": hip_id}

    return analyze_pma(
        gaia_id, hip_id,
        distance_pc=target.get("distance_pc"),
        stellar_mass_msun=target.get("mass_msun"),
    )


def _summarize_data(report):
    """Create a data availability summary from the report sections."""
    methods = []
    notes = []

    # Stellar properties
    sp = report.get("stellar_properties", {})
    if sp and sp.get("status") == "ok":
        methods.append("Gaia DR3 stellar properties")

    # Known planets
    kp = report.get("known_planets", {})
    if kp:
        n = kp.get("n_planets", 0)
        if n > 0:
            notes.append(f"{n} confirmed planet(s) in NASA Exoplanet Archive")
        else:
            notes.append("No confirmed planets in literature")

    # Light curve
    lc = report.get("lightcurve", {})
    if lc and lc.get("available"):
        n_sectors = lc.get("n_sectors", 0)
        methods.append(f"TESS photometry ({n_sectors} sector(s))")
    else:
        notes.append("No TESS light curve available")

    # RV
    rv = report.get("rv_data", {})
    if rv and rv.get("n_measurements", 0) > 0:
        instruments = rv.get("instruments", [])
        methods.append(f"RV time series ({rv['n_measurements']} measurements, "
                       f"{', '.join(instruments)})")
    else:
        notes.append("No public RV time series found")

    # PMa
    pma = report.get("proper_motion_anomaly", {})
    if pma and pma.get("status") == "ok":
        methods.append("Gaia-Hipparcos proper motion anomaly")
        pma_data = pma.get("pma", {})
        if pma_data.get("significant"):
            notes.append(f"SIGNIFICANT PMa detected (SNR={pma_data.get('pma_snr', 'N/A')})")

    return {
        "methods_available": methods,
        "n_methods": len(methods),
        "notes": notes,
    }


# --- Formatting helpers ---

def _fmt(value, decimals=4):
    """Format a value for display, handling None."""
    if value is None:
        return "N/A"
    if isinstance(value, float):
        return f"{value:.{decimals}f}"
    return str(value)


def _json_serializer(obj):
    """JSON serializer for numpy types and other non-standard types."""
    import numpy as np
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
