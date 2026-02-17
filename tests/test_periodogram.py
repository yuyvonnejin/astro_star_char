"""Unit tests for Module 4b: period detection and variability classification.

These tests use synthetic signals and do NOT require network access.
"""

import numpy as np
import pytest
from src.periodogram import detect_period, classify_variability, analyze_lightcurve


# -- Synthetic signal helpers ------------------------------------------------

def make_sinusoidal(period=3.5, amplitude=0.01, n_points=5000, baseline=30.0,
                    noise_level=0.001):
    """Create a synthetic sinusoidal light curve."""
    time = np.linspace(0, baseline, n_points)
    flux = 1.0 + amplitude * np.sin(2 * np.pi * time / period)
    flux += np.random.default_rng(42).normal(0, noise_level, n_points)
    return time, flux


def make_noise(n_points=5000, baseline=30.0, noise_level=0.001):
    """Create a pure noise light curve (no periodic signal)."""
    time = np.linspace(0, baseline, n_points)
    flux = 1.0 + np.random.default_rng(42).normal(0, noise_level, n_points)
    return time, flux


# -- Period detection tests --------------------------------------------------

class TestDetectPeriod:

    def test_recovers_known_period(self):
        """Lomb-Scargle should recover a strong sinusoidal signal."""
        true_period = 3.5
        time, flux = make_sinusoidal(period=true_period, amplitude=0.01, noise_level=0.001)
        result = detect_period(time, flux, min_period=0.5, max_period=15.0)

        assert result["period_days"] is not None
        assert abs(result["period_days"] - true_period) < 0.1  # within 0.1 day
        assert result["period_fap"] < 0.001  # confident detection
        assert result["detection_method"] == "lomb_scargle"

    def test_recovers_short_period(self):
        """Should detect a short-period signal (~0.5 day)."""
        true_period = 0.5
        time, flux = make_sinusoidal(period=true_period, amplitude=0.005,
                                      n_points=10000, baseline=10.0, noise_level=0.0005)
        result = detect_period(time, flux, min_period=0.1, max_period=5.0)

        assert result["period_days"] is not None
        assert abs(result["period_days"] - true_period) < 0.05

    def test_noise_has_high_fap(self):
        """Pure noise should not produce a confident period detection."""
        time, flux = make_noise(n_points=5000, baseline=30.0, noise_level=0.001)
        result = detect_period(time, flux, min_period=0.5, max_period=15.0)

        # FAP should be high (not a confident detection)
        assert result["period_fap"] > 0.01

    def test_amplitude_estimation(self):
        """Amplitude should be roughly correct for a clean sinusoid."""
        amplitude = 0.005  # 5 ppt peak
        time, flux = make_sinusoidal(period=3.0, amplitude=amplitude, noise_level=0.0001)
        result = detect_period(time, flux, min_period=0.5, max_period=15.0)

        # Peak-to-peak is ~2*amplitude = 10 ppt
        expected_ptp = 2 * amplitude * 1000  # 10 ppt
        assert result["amplitude_ppt"] is not None
        assert abs(result["amplitude_ppt"] - expected_ptp) < 3.0  # within 3 ppt

    def test_empty_input(self):
        """Should handle degenerate input gracefully."""
        result = detect_period(np.array([0.0]), np.array([1.0]))
        assert result["period_days"] is None

    def test_max_period_capped_at_half_baseline(self):
        """max_period should be capped at half the time baseline."""
        time, flux = make_sinusoidal(period=3.0, baseline=10.0)
        result = detect_period(time, flux, min_period=0.5, max_period=100.0)
        # Should still work -- max_period capped to 5.0 internally
        assert result["period_days"] is not None

    def test_returns_frequency_and_power_arrays(self):
        """Result should include frequency and power arrays for plotting."""
        time, flux = make_sinusoidal(period=3.0)
        result = detect_period(time, flux, min_period=0.5, max_period=15.0)

        assert len(result["frequency"]) > 0
        assert len(result["power"]) == len(result["frequency"])


# -- Classification tests ---------------------------------------------------

class TestClassifyVariability:

    def test_periodic_classification(self):
        """FAP < 0.001 should classify as periodic."""
        result = classify_variability(
            {"period_days": 3.5, "period_fap": 1e-10, "period_power": 0.9, "amplitude_ppt": 5.0}
        )
        assert result["variability_class"] == "periodic"
        assert result["variability_flag"] == "ok"

    def test_possible_periodic_classification(self):
        """0.001 <= FAP < 0.01 should classify as possible_periodic."""
        result = classify_variability(
            {"period_days": 3.5, "period_fap": 0.005, "period_power": 0.3, "amplitude_ppt": 1.0}
        )
        assert result["variability_class"] == "possible_periodic"

    def test_non_variable_classification(self):
        """FAP >= 0.01 should classify as non_variable."""
        result = classify_variability(
            {"period_days": 3.5, "period_fap": 0.5, "period_power": 0.05, "amplitude_ppt": 0.1}
        )
        assert result["variability_class"] == "non_variable"

    def test_few_points_flag(self):
        """Should flag when n_points is low."""
        result = classify_variability(
            {"period_days": 3.5, "period_fap": 1e-5, "period_power": 0.7, "amplitude_ppt": 5.0},
            n_points=50,
        )
        assert result["variability_flag"] == "few_points"

    def test_short_baseline_flag(self):
        """Should flag when baseline is very short."""
        result = classify_variability(
            {"period_days": 3.5, "period_fap": 1e-5, "period_power": 0.7, "amplitude_ppt": 5.0},
            time_baseline_days=0.5,
        )
        assert result["variability_flag"] == "short_baseline"

    def test_no_data_result(self):
        """Should handle empty period result."""
        result = classify_variability(
            {"period_days": None, "period_fap": None, "period_power": None, "amplitude_ppt": None}
        )
        assert result["variability_class"] == "undetermined"


# -- Integrated analysis test ------------------------------------------------

class TestAnalyzeLightcurve:

    def test_full_analysis_pipeline(self):
        """analyze_lightcurve should chain detection + classification."""
        true_period = 5.0
        time, flux = make_sinusoidal(period=true_period, amplitude=0.01,
                                      n_points=8000, baseline=50.0, noise_level=0.001)
        lc_data = {
            "time": time,
            "flux": flux,
            "flux_flat": flux,
            "flux_err": np.full_like(flux, 0.001),
            "mission": "Synthetic",
            "author": "Test",
            "n_sectors": 1,
            "n_points_clean": len(time),
            "cadence_s": 120.0,
            "time_baseline_days": 50.0,
        }

        result = analyze_lightcurve(lc_data)

        assert result["lightcurve_available"] is True
        assert result["variability_class"] == "periodic"
        assert abs(result["period_days"] - true_period) < 0.2
        assert result["period_fap"] < 0.001
        assert result["lc_mission"] == "Synthetic"
        # Should NOT contain large arrays
        assert "frequency" not in result
        assert "power" not in result
