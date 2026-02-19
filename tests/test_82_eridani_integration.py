"""Integration test: Generate full Phase 7 report for 82 G. Eridani (HD 20794).

This test requires network access to query Gaia, SIMBAD, NASA Exoplanet Archive,
DACE, and Hipparcos catalogs.

Run with: .\\venv\\Scripts\\python.exe tests\\test_82_eridani_integration.py
"""

import json
import logging
import sys
from pathlib import Path

# Setup logging
LOG_DIR = Path(__file__).resolve().parent.parent / "logs"
LOG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "test_82_eridani.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.targets import get_target, resolve_target_ids
from src.rv_data import query_known_planets
from src.proper_motion import analyze_pma
from src.target_report import generate_report, format_report_markdown, save_report


def test_resolve_ids():
    """Step 1: Resolve Gaia DR3 and TIC IDs via SIMBAD."""
    logger.info("=" * 60)
    logger.info("Step 1: Resolving identifiers for 82 G. Eridani")
    logger.info("=" * 60)

    target = get_target("82 G. Eridani")
    assert target is not None, "82 G. Eridani not in catalog"

    resolved = resolve_target_ids(target)
    logger.info("Resolved Gaia DR3 ID: %s", resolved.get("gaia_dr3_id"))
    logger.info("Resolved TIC: %s", resolved.get("tic"))

    return resolved


def test_known_planets():
    """Step 2: Query NASA Exoplanet Archive for known planets."""
    logger.info("=" * 60)
    logger.info("Step 2: Querying known planets for HD 20794")
    logger.info("=" * 60)

    planets = query_known_planets("HD 20794")
    logger.info("Found %d planet(s)", len(planets))
    for p in planets:
        logger.info("  %s: P=%.2f d, M=%.2f Mearth, method=%s",
                     p.get("pl_name"), p.get("period_days", 0),
                     p.get("mass_earth", 0), p.get("discovery_method"))

    return planets


def test_pma(target):
    """Step 3: Compute proper motion anomaly."""
    logger.info("=" * 60)
    logger.info("Step 3: Proper motion anomaly analysis")
    logger.info("=" * 60)

    gaia_id = target.get("gaia_dr3_id")
    hip_id = target.get("hip")

    if gaia_id is None:
        logger.warning("No Gaia ID resolved, skipping PMa")
        return None

    result = analyze_pma(
        gaia_id, hip_id,
        distance_pc=target.get("distance_pc"),
        stellar_mass_msun=target.get("mass_msun"),
    )

    pma = result.get("pma", {})
    logger.info("PMa total: %.4f mas/yr, SNR: %.1f, significant: %s",
                pma.get("pma_total_mas_yr", 0),
                pma.get("pma_snr", 0),
                pma.get("significant", "N/A"))
    logger.info("RUWE: %s", result.get("ruwe"))

    return result


def test_full_report(target):
    """Step 4: Generate full report (excluding light curve download for speed)."""
    logger.info("=" * 60)
    logger.info("Step 4: Generating full target report")
    logger.info("=" * 60)

    report = generate_report(
        target,
        include_stellar=True,
        include_lightcurve=True,  # Just checks availability, no download
        include_transit=False,     # Skip transit to avoid long download
        include_rv=True,
        include_pma=True,
        include_known_planets=True,
    )

    # Print summary
    ds = report.get("data_summary", {})
    logger.info("Methods available: %s", ds.get("methods_available", []))
    logger.info("Notes: %s", ds.get("notes", []))

    # Save report
    paths = save_report(report)
    logger.info("Report saved: %s", paths)

    # Also print markdown
    md = format_report_markdown(report)
    logger.info("\n%s", md)

    return report


def main():
    logger.info("Integration test: 82 G. Eridani (HD 20794)")
    logger.info("This test requires network access.")

    # Step 1: Resolve IDs
    target = test_resolve_ids()

    # Step 2: Known planets
    planets = test_known_planets()

    # Step 3: PMa
    pma_result = test_pma(target)

    # Step 4: Full report
    report = test_full_report(target)

    # Validation checks
    logger.info("=" * 60)
    logger.info("Validation summary")
    logger.info("=" * 60)

    checks = []
    checks.append(("Gaia ID resolved", target.get("gaia_dr3_id") is not None))
    checks.append(("Known planets found", len(planets) >= 1))
    checks.append(("PMa computed", pma_result is not None and pma_result.get("status") == "ok"))
    checks.append(("Report generated", report is not None))
    checks.append(("Stellar properties", report.get("stellar_properties", {}).get("status") == "ok"))

    all_passed = True
    for name, passed in checks:
        status = "PASS" if passed else "FAIL"
        logger.info("  [%s] %s", status, name)
        if not passed:
            all_passed = False

    if all_passed:
        logger.info("All validation checks passed.")
    else:
        logger.warning("Some validation checks failed.")

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
