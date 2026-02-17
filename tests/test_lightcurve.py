"""Unit tests for Module 4a: light curve retrieval.

Tests marked with @pytest.mark.network require MAST network access.
Run with: pytest -m network tests/test_lightcurve.py
Skip with: pytest -m "not network"
"""

import numpy as np
import pytest
from src.lightcurve import clean_lightcurve, flatten_lightcurve


# -- Offline tests (no network) ---------------------------------------------

class TestCleanLightcurve:

    def _make_lc_data(self, n=1000):
        """Helper: create a simple light curve dict."""
        rng = np.random.default_rng(42)
        time = np.linspace(0, 30, n)
        flux = 1.0 + rng.normal(0, 0.001, n)
        flux_err = np.full(n, 0.001)
        return {
            "time": time,
            "flux": flux,
            "flux_err": flux_err,
            "mission": "Test",
            "author": "Test",
            "n_sectors": 1,
            "n_points_raw": n,
            "cadence_s": 120.0,
            "time_baseline_days": 30.0,
        }

    def test_removes_nans(self):
        """Should remove NaN flux values."""
        lc = self._make_lc_data(100)
        lc["flux"][10] = np.nan
        lc["flux"][50] = np.nan

        cleaned = clean_lightcurve(lc)
        assert cleaned["n_points_clean"] == 98
        assert not np.any(np.isnan(cleaned["flux"]))

    def test_removes_outliers(self):
        """Should remove flux outliers beyond sigma threshold."""
        lc = self._make_lc_data(1000)
        # Inject a huge outlier
        lc["flux"][500] = 2.0  # way above 5-sigma

        cleaned = clean_lightcurve(lc, sigma_clip=5.0)
        assert cleaned["n_points_clean"] < 1000
        assert np.max(cleaned["flux"]) < 2.0

    def test_preserves_clean_data(self):
        """Clean data should pass through mostly unchanged."""
        lc = self._make_lc_data(1000)
        cleaned = clean_lightcurve(lc, sigma_clip=5.0)
        # Should keep nearly all points (only statistical outliers removed)
        assert cleaned["n_points_clean"] >= 990

    def test_arrays_have_consistent_length(self):
        """time, flux, flux_err arrays should all have same length after cleaning."""
        lc = self._make_lc_data(500)
        lc["flux"][0] = np.nan
        lc["flux"][100] = 99.0  # outlier

        cleaned = clean_lightcurve(lc)
        assert len(cleaned["time"]) == len(cleaned["flux"])
        assert len(cleaned["time"]) == len(cleaned["flux_err"])
        assert len(cleaned["time"]) == cleaned["n_points_clean"]


class TestFlattenLightcurve:

    def _make_trended_lc(self, n=2000):
        """Create light curve with a linear trend + periodic signal."""
        time = np.linspace(0, 30, n)
        # Linear trend (spacecraft drift) + sinusoidal signal
        trend = 0.01 * (time - 15)  # slope of 0.01/day
        signal = 0.005 * np.sin(2 * np.pi * time / 3.0)
        flux = 1.0 + trend + signal
        flux_err = np.full(n, 0.001)
        return {
            "time": time,
            "flux": flux,
            "flux_err": flux_err,
            "mission": "Test",
            "author": "Test",
            "n_sectors": 1,
            "n_points_raw": n,
            "n_points_clean": n,
            "cadence_s": 120.0,
            "time_baseline_days": 30.0,
        }

    def test_removes_trend(self):
        """Flattening should remove the linear trend."""
        lc = self._make_trended_lc()
        original_range = np.ptp(lc["flux"])

        flat = flatten_lightcurve(lc, window_length=501)
        flat_range = np.ptp(flat["flux_flat"])

        # Flattened range should be much smaller than original
        assert flat_range < original_range * 0.5

    def test_preserves_periodic_signal(self):
        """Flattening should preserve the periodic signal."""
        lc = self._make_trended_lc(n=5000)
        flat = flatten_lightcurve(lc, window_length=1001)

        # Check that the periodic signal survives by looking at RMS
        flux_flat = flat["flux_flat"]
        rms = np.std(flux_flat)
        # Signal amplitude is 0.005, so RMS should be > 0.002
        assert rms > 0.002

    def test_output_has_flux_flat_key(self):
        """Output should contain flux_flat and window_length."""
        lc = self._make_trended_lc(n=500)
        flat = flatten_lightcurve(lc, window_length=101)
        assert "flux_flat" in flat
        assert "window_length" in flat
        assert len(flat["flux_flat"]) == len(lc["flux"])

    def test_handles_short_data(self):
        """Should handle data shorter than window_length without crashing."""
        lc = self._make_trended_lc(n=50)
        flat = flatten_lightcurve(lc, window_length=1001)
        assert "flux_flat" in flat


# -- Network tests (require MAST access) ------------------------------------

@pytest.mark.network
class TestSearchLightcurve:

    def test_search_kepler_target(self):
        """Should find Kepler light curves for KIC 6922244."""
        from src.lightcurve import search_lightcurve

        result = search_lightcurve("KIC 6922244", mission="Kepler")
        assert result is not None
        assert result["n_available"] > 0
        assert "Kepler" in result["mission"]

    def test_search_nonexistent_target(self):
        """Should return None for a target with no data."""
        from src.lightcurve import search_lightcurve

        result = search_lightcurve("NONEXISTENT_STAR_XYZ_12345")
        assert result is None


@pytest.mark.network
class TestRetrieveLightcurve:

    def test_full_retrieve_kepler(self):
        """Full retrieval pipeline for a known Kepler target (limited to 2 quarters)."""
        from src.lightcurve import retrieve_lightcurve

        lc = retrieve_lightcurve("KIC 6922244", mission="Kepler", max_sectors=2)
        assert lc is not None
        assert len(lc["time"]) > 100
        assert len(lc["flux_flat"]) == len(lc["time"])
        assert lc["n_sectors"] <= 2
        assert lc["time_baseline_days"] > 0
