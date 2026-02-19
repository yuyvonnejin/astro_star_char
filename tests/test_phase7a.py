"""Tests for Phase 7a modules: targets, rv_data, proper_motion, target_report.

Unit tests use synthetic/mock data where possible to avoid network dependency.
Integration tests (marked with comments) require network access.
"""

import numpy as np
import pytest

from src.targets import get_targets, get_target, TARGET_CATALOG
from src.rv_data import (
    rv_periodogram, rv_detection_limit, rv_to_planet_mass,
    _generate_search_names, _find_peaks,
)
from src.proper_motion import compute_pma, pma_companion_mass, assess_ruwe


# ============================================================
# Target catalog tests
# ============================================================

class TestTargetCatalog:
    """Tests for src.targets module."""

    def test_catalog_has_10_targets(self):
        """All 10 Phase 7 targets are present."""
        targets = get_targets()
        assert len(targets) == 10

    def test_all_required_fields_present(self):
        """Each target has all required fields."""
        required = ["name", "hd", "hip", "spectral_type", "teff_k",
                     "distance_pc", "mass_msun", "tier"]
        for target in TARGET_CATALOG:
            for field in required:
                assert field in target, f"{target['name']} missing field '{field}'"
                assert target[field] is not None, f"{target['name']} has None '{field}'"

    def test_target_lookup_by_name(self):
        """Lookup by common name works."""
        target = get_target("Tau Ceti")
        assert target is not None
        assert target["hd"] == "HD 10700"

    def test_target_lookup_by_hd(self):
        """Lookup by HD number works."""
        target = get_target("HD 20794")
        assert target is not None
        assert target["name"] == "82 G. Eridani"

    def test_target_lookup_case_insensitive(self):
        """Lookup is case insensitive."""
        target = get_target("tau ceti")
        assert target is not None
        assert target["name"] == "Tau Ceti"

    def test_target_lookup_not_found(self):
        """Lookup returns None for unknown target."""
        target = get_target("Betelgeuse")
        assert target is None

    def test_teff_in_range(self):
        """All targets have Teff within the G/early-F selection range."""
        for target in TARGET_CATALOG:
            assert 5200 <= target["teff_k"] <= 6100, (
                f"{target['name']} Teff={target['teff_k']} out of selection range"
            )

    def test_distance_within_25pc(self):
        """All targets are within 25 pc (maximum for bonus targets)."""
        for target in TARGET_CATALOG:
            assert target["distance_pc"] <= 25.0, (
                f"{target['name']} at {target['distance_pc']} pc exceeds 25 pc limit"
            )

    def test_priority_ordering(self):
        """Targets are roughly ordered by priority (distance)."""
        targets = get_targets()
        # First 8 should be within 10 pc (the core selection)
        for t in targets[:8]:
            assert t["distance_pc"] <= 10.0, (
                f"{t['name']} at {t['distance_pc']} pc should be within 10 pc"
            )

    def test_82_eridani_is_validation_target(self):
        """82 G. Eridani (HD 20794) is in the catalog with confirmed planets."""
        target = get_target("82 G. Eridani")
        assert target is not None
        assert target["known_planets"] >= 3
        assert target["hd"] == "HD 20794"

    def test_search_names_generation(self):
        """_generate_search_names produces expected aliases."""
        names = _generate_search_names("82 G. Eridani")
        assert "HD 20794" in names
        assert "82 G. Eridani" in names


# ============================================================
# RV analysis tests (using synthetic data)
# ============================================================

class TestRVAnalysis:
    """Tests for src.rv_data analysis functions."""

    def test_rv_periodogram_synthetic(self):
        """RV periodogram detects injected sinusoidal signal."""
        np.random.seed(42)
        n_obs = 200
        true_period = 50.0  # days
        k_amplitude = 3.0  # m/s

        time = np.sort(np.random.uniform(0, 500, n_obs))
        rv = k_amplitude * np.sin(2 * np.pi * time / true_period)
        rv += np.random.normal(0, 0.5, n_obs)  # noise
        rv_err = np.full(n_obs, 0.5)

        result = rv_periodogram(time, rv, rv_err, min_period=5.0, max_period=200.0)
        assert result["best_period"] is not None
        # Check that detected period is within 10% of true period
        assert abs(result["best_period"] - true_period) / true_period < 0.10

    def test_rv_periodogram_no_signal(self):
        """RV periodogram returns result even with pure noise."""
        np.random.seed(123)
        time = np.sort(np.random.uniform(0, 500, 100))
        rv = np.random.normal(0, 1.0, 100)
        rv_err = np.ones(100)

        result = rv_periodogram(time, rv, rv_err)
        assert result["best_period"] is not None
        # FAP should be high for noise
        assert result["fap"] > 0.01

    def test_rv_periodogram_empty(self):
        """RV periodogram handles edge cases."""
        result = rv_periodogram(np.array([1.0, 1.0]), np.array([0.0, 0.0]))
        assert result["best_period"] is None

    def test_rv_detection_limit(self):
        """Detection limit produces reasonable values."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        rv_err = np.full(300, 1.0)  # 1 m/s precision

        result = rv_detection_limit(time, rv_err)
        assert len(result["k_min_ms"]) > 0
        assert result["sigma_eff_ms"] == pytest.approx(1.0, rel=0.01)
        # Detection limit should be finite for well-sampled periods
        finite_mask = np.isfinite(result["k_min_ms"])
        assert np.any(finite_mask)
        # At well-sampled periods, K_min should be order of sigma/sqrt(N)
        k_min_short = result["k_min_ms"][finite_mask][0]
        assert k_min_short < 5.0  # should be much less than 5 m/s with 300 obs

    def test_rv_to_planet_mass(self):
        """Planet mass conversion is physically reasonable."""
        # Earth around the Sun: K ~ 0.089 m/s, P ~ 365.25 days
        mass = rv_to_planet_mass(0.089, 365.25, stellar_mass_msun=1.0)
        # Should be close to 1 Earth mass (within factor of 2 for this approximation)
        assert 0.3 < mass < 3.0, f"Expected ~1 Mearth, got {mass:.2f}"

    def test_rv_to_planet_mass_jupiter(self):
        """Jupiter-mass conversion is reasonable."""
        # Jupiter around the Sun: K ~ 12.5 m/s, P ~ 4332 days
        mass = rv_to_planet_mass(12.5, 4332.0, stellar_mass_msun=1.0)
        mass_mjup = mass / 317.8  # convert Mearth to Mjup
        assert 0.3 < mass_mjup < 3.0, f"Expected ~1 Mjup, got {mass_mjup:.2f}"


# ============================================================
# Proper motion anomaly tests (using synthetic data)
# ============================================================

class TestProperMotionAnomaly:
    """Tests for src.proper_motion analysis functions."""

    def test_pma_zero_difference(self):
        """PMa is zero when Gaia and Hipparcos agree."""
        gaia = {"pmra": 100.0, "pmdec": -200.0, "pmra_error": 0.01, "pmdec_error": 0.01}
        hip = {"pmra": 100.0, "pmdec": -200.0, "pmra_error": 0.5, "pmdec_error": 0.5}

        result = compute_pma(gaia, hip)
        assert result["pma_total_mas_yr"] == pytest.approx(0.0, abs=0.001)
        assert not result["significant"]

    def test_pma_significant_detection(self):
        """PMa is significant with large difference."""
        gaia = {"pmra": 100.0, "pmdec": -200.0, "pmra_error": 0.01, "pmdec_error": 0.01}
        hip = {"pmra": 110.0, "pmdec": -200.0, "pmra_error": 0.5, "pmdec_error": 0.5}

        result = compute_pma(gaia, hip)
        assert result["pma_total_mas_yr"] == pytest.approx(10.0, abs=0.1)
        assert result["significant"]
        assert result["pma_snr"] > 3.0

    def test_pma_companion_mass_physically_reasonable(self):
        """Companion mass estimate is order-of-magnitude correct."""
        # 1 mas/yr PMa at 5 pc for a solar-mass star
        result = pma_companion_mass(1.0, 5.0, 1.0, separation_au=np.array([5.0]))
        # Should give a substellar to planetary mass at 5 AU
        assert result["mass_mjup"][0] > 0
        assert result["mass_mjup"][0] < 1000  # less than 1 solar mass in Mjup

    def test_pma_companion_mass_closer_stronger(self):
        """Closer companion needs higher mass to produce same PMa."""
        result = pma_companion_mass(1.0, 5.0, 1.0,
                                     separation_au=np.array([2.0, 5.0, 10.0]))
        # At smaller separation, same PMa implies *smaller* companion mass
        # (because closer companions produce more reflex motion per unit mass)
        assert result["mass_mjup"][0] < result["mass_mjup"][1]
        assert result["mass_mjup"][1] < result["mass_mjup"][2]

    def test_ruwe_normal(self):
        """Normal RUWE is correctly classified."""
        result = assess_ruwe(1.05)
        assert not result["companion_hint"]
        assert "normal" in result["interpretation"]

    def test_ruwe_elevated(self):
        """Elevated RUWE triggers companion hint."""
        result = assess_ruwe(2.3)
        assert result["companion_hint"]

    def test_ruwe_none(self):
        """None RUWE handled gracefully."""
        result = assess_ruwe(None)
        assert not result["companion_hint"]


# ============================================================
# Report structure tests
# ============================================================

class TestReportStructure:
    """Tests for report generation structure (no network calls)."""

    def test_target_summary(self):
        """_target_summary extracts correct fields."""
        from src.target_report import _target_summary
        target = get_target("Tau Ceti")
        summary = _target_summary(target)
        assert summary["name"] == "Tau Ceti"
        assert summary["hd"] == "HD 10700"
        assert summary["spectral_type"] == "G8.5V"

    def test_format_report_markdown(self):
        """format_report_markdown produces valid markdown."""
        from src.target_report import format_report_markdown

        # Minimal report
        report = {
            "target": {
                "name": "Test Star",
                "hd": "HD 99999",
                "hip": 99999,
                "spectral_type": "G2V",
            },
            "generated_utc": "2026-02-19T00:00:00",
            "sections_included": [],
            "data_summary": {
                "methods_available": ["Gaia DR3"],
                "n_methods": 1,
                "notes": ["Test note"],
            },
        }

        md = format_report_markdown(report)
        assert "# Target Report: Test Star" in md
        assert "HD 99999" in md
        assert "Test note" in md

    def test_json_serializer(self):
        """JSON serializer handles numpy types."""
        from src.target_report import _json_serializer
        import json

        data = {
            "array": np.array([1.0, 2.0, 3.0]),
            "int": np.int64(42),
            "float": np.float64(3.14),
        }
        # Should not raise
        json_str = json.dumps(data, default=_json_serializer)
        assert "42" in json_str
        assert "3.14" in json_str
