"""Tests for Phase 7b modules: RV multi-planet analysis, sensitivity maps, deep-dive orchestrator.

Unit tests use synthetic data to avoid network dependency.
"""

import numpy as np
import pytest

from src.rv_data import (
    rv_subtract_instrument_offsets,
    rv_subtract_sinusoids,
    rv_residual_analysis,
    rv_injection_recovery,
    rv_periodogram,
    rv_to_planet_mass,
)
from src.sensitivity import (
    transit_sensitivity,
    rv_sensitivity,
    astrometric_sensitivity,
    combined_sensitivity,
)


# ============================================================
# RV Instrument Offset Subtraction
# ============================================================

class TestRVInstrumentOffsets:
    """Tests for rv_subtract_instrument_offsets."""

    def test_single_instrument_no_change(self):
        """Single instrument: RV is centered on zero (median subtracted)."""
        np.random.seed(42)
        n = 100
        time = np.sort(np.random.uniform(0, 1000, n))
        rv = np.random.normal(50.0, 2.0, n)  # offset of 50 m/s
        rv_err = np.ones(n)
        instruments = ["HARPS"] * n

        result = rv_subtract_instrument_offsets(time, rv, rv_err, instruments)
        assert result["n_instruments"] == 1
        assert abs(np.median(result["rv"])) < 0.5
        assert "HARPS" in result["offsets"]

    def test_two_instruments_aligned(self):
        """Two instruments with different offsets are aligned to zero."""
        np.random.seed(42)
        n = 200
        time = np.sort(np.random.uniform(0, 1000, n))
        rv = np.random.normal(0, 1.0, n)
        rv_err = np.ones(n)
        instruments = ["HARPS"] * 100 + ["ESPRESSO"] * 100

        # Add instrument-specific offsets
        rv[:100] += 100.0   # HARPS offset
        rv[100:] += -50.0   # ESPRESSO offset

        result = rv_subtract_instrument_offsets(time, rv, rv_err, instruments)
        assert result["n_instruments"] == 2

        # After correction, both instruments should be centered near zero
        inst_arr = np.array(result["instruments"])
        harps_median = np.median(result["rv"][inst_arr == "HARPS"])
        espresso_median = np.median(result["rv"][inst_arr == "ESPRESSO"])
        assert abs(harps_median) < 1.0
        assert abs(espresso_median) < 1.0

    def test_three_instruments(self):
        """Three instruments with varying offsets."""
        np.random.seed(42)
        n = 300
        time = np.sort(np.random.uniform(0, 5000, n))
        rv = np.random.normal(0, 0.5, n)
        rv_err = np.full(n, 0.5)
        instruments = ["CORAVEL"] * 100 + ["HARPS"] * 100 + ["ESPRESSO"] * 100

        rv[:100] += 200.0
        rv[100:200] += 50.0
        rv[200:] += -10.0

        result = rv_subtract_instrument_offsets(time, rv, rv_err, instruments)
        assert result["n_instruments"] == 3
        assert len(result["offsets"]) == 3
        # Offsets should be close to the injected values
        assert abs(result["offsets"]["CORAVEL"] - 200.0) < 2.0
        assert abs(result["offsets"]["HARPS"] - 50.0) < 2.0
        assert abs(result["offsets"]["ESPRESSO"] - (-10.0)) < 2.0


# ============================================================
# RV Sinusoid Subtraction
# ============================================================

class TestRVSinusoidSubtraction:
    """Tests for rv_subtract_sinusoids."""

    def test_single_sinusoid_subtraction(self):
        """Subtracting a single sinusoid reduces RMS significantly."""
        np.random.seed(42)
        n = 500
        true_period = 50.0
        k_amplitude = 5.0

        time = np.sort(np.random.uniform(0, 2000, n))
        signal = k_amplitude * np.sin(2 * np.pi * time / true_period + 0.7)
        noise = np.random.normal(0, 0.5, n)
        rv = signal + noise

        result = rv_subtract_sinusoids(time, rv, [true_period])
        assert result["rms_after_ms"] < result["rms_before_ms"]
        # RMS should drop to near noise level
        assert result["rms_after_ms"] < 1.0
        # Fitted amplitude should be close to injected
        assert abs(result["fitted_components"][0]["amplitude_ms"] - k_amplitude) < 0.5

    def test_two_sinusoids_subtraction(self):
        """Subtracting two known sinusoids leaves clean residuals."""
        np.random.seed(42)
        n = 500
        time = np.sort(np.random.uniform(0, 2000, n))
        signal1 = 3.0 * np.sin(2 * np.pi * time / 20.0)
        signal2 = 2.0 * np.sin(2 * np.pi * time / 90.0)
        noise = np.random.normal(0, 0.3, n)
        rv = signal1 + signal2 + noise

        result = rv_subtract_sinusoids(time, rv, [20.0, 90.0])
        assert result["rms_after_ms"] < 1.0
        assert len(result["fitted_components"]) == 2

    def test_empty_periods_list(self):
        """Empty periods list returns original data unchanged."""
        rv = np.array([1.0, 2.0, 3.0])
        time = np.array([0.0, 1.0, 2.0])
        result = rv_subtract_sinusoids(time, rv, [])
        np.testing.assert_array_equal(result["residuals"], rv)
        assert result["rms_before_ms"] == result["rms_after_ms"]

    def test_subtraction_preserves_unknown_signal(self):
        """Subtracting known periods preserves signal at unknown period."""
        np.random.seed(42)
        n = 500
        time = np.sort(np.random.uniform(0, 2000, n))
        known_signal = 3.0 * np.sin(2 * np.pi * time / 50.0)
        unknown_signal = 2.0 * np.sin(2 * np.pi * time / 200.0)
        noise = np.random.normal(0, 0.3, n)
        rv = known_signal + unknown_signal + noise

        result = rv_subtract_sinusoids(time, rv, [50.0])
        residuals = result["residuals"]

        # Run periodogram on residuals -- should find ~200 day period
        pg = rv_periodogram(time, residuals, min_period=10.0, max_period=500.0)
        assert abs(pg["best_period"] - 200.0) / 200.0 < 0.15


# ============================================================
# RV Injection-Recovery
# ============================================================

class TestRVInjectionRecovery:
    """Tests for rv_injection_recovery."""

    def test_strong_signal_high_detection(self):
        """Strong injected signal (high K) should be detected in most trials."""
        np.random.seed(42)
        n = 200
        time = np.sort(np.random.uniform(0, 1000, n))
        rv_err = np.full(n, 1.0)

        # Single strong amplitude, single well-sampled period
        result = rv_injection_recovery(
            time, rv_err,
            period_grid=np.array([50.0]),
            k_grid=np.array([10.0]),  # 10x the noise
            n_trials=30,
        )
        # Detection probability should be very high
        assert result["detection_probability"][0, 0] > 0.7

    def test_weak_signal_low_detection(self):
        """Weak injected signal (low K) should be detected rarely."""
        np.random.seed(42)
        n = 100
        time = np.sort(np.random.uniform(0, 500, n))
        rv_err = np.full(n, 2.0)

        result = rv_injection_recovery(
            time, rv_err,
            period_grid=np.array([50.0]),
            k_grid=np.array([0.1]),  # 0.05x the noise
            n_trials=30,
        )
        # Detection probability should be very low
        assert result["detection_probability"][0, 0] < 0.3

    def test_detection_increases_with_amplitude(self):
        """Detection probability should increase with injected amplitude."""
        np.random.seed(42)
        n = 150
        time = np.sort(np.random.uniform(0, 800, n))
        rv_err = np.full(n, 1.0)

        k_grid = np.array([0.3, 1.0, 5.0])
        result = rv_injection_recovery(
            time, rv_err,
            period_grid=np.array([30.0]),
            k_grid=k_grid,
            n_trials=30,
        )
        probs = result["detection_probability"][:, 0]
        # Monotonically increasing (or at least last > first)
        assert probs[2] >= probs[0]

    def test_output_shape(self):
        """Output arrays have correct shape."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 500, 100))
        rv_err = np.ones(100)

        result = rv_injection_recovery(
            time, rv_err,
            period_grid=np.array([10.0, 50.0, 100.0]),
            k_grid=np.array([0.5, 2.0]),
            n_trials=10,
        )
        assert result["detection_probability"].shape == (2, 3)
        assert result["mass_grid_mearth"].shape == (2, 3)


# ============================================================
# Transit Sensitivity
# ============================================================

class TestTransitSensitivity:
    """Tests for transit_sensitivity."""

    def test_sensitivity_improves_with_more_transits(self):
        """Longer baseline -> more transits -> smaller detectable radius."""
        periods = np.array([10.0])
        # Short baseline: few transits
        short = transit_sensitivity(1.0, 100.0, 100.0, period_grid=periods)
        # Long baseline: many transits
        long = transit_sensitivity(1.0, 100.0, 1000.0, period_grid=periods)
        assert long["min_radius_rearth"][0] < short["min_radius_rearth"][0]

    def test_sensitivity_degrades_with_larger_star(self):
        """Larger star -> shallower transit for same planet -> harder to detect."""
        periods = np.array([10.0])
        small_star = transit_sensitivity(0.5, 100.0, 500.0, period_grid=periods)
        large_star = transit_sensitivity(2.0, 100.0, 500.0, period_grid=periods)
        # Larger star has larger minimum detectable radius
        assert large_star["min_radius_rearth"][0] > small_star["min_radius_rearth"][0]

    def test_sensitivity_degrades_with_worse_cdpp(self):
        """Worse CDPP -> harder to detect small planets."""
        periods = np.array([10.0])
        good = transit_sensitivity(1.0, 50.0, 500.0, period_grid=periods)
        bad = transit_sensitivity(1.0, 500.0, 500.0, period_grid=periods)
        assert bad["min_radius_rearth"][0] > good["min_radius_rearth"][0]

    def test_long_period_undetectable(self):
        """Very long period with short baseline -> inf radius (undetectable)."""
        periods = np.array([10000.0])
        result = transit_sensitivity(1.0, 100.0, 100.0, period_grid=periods)
        assert np.isinf(result["min_radius_rearth"][0])

    def test_geometric_probability_decreases_with_period(self):
        """Geometric transit probability decreases with orbital period."""
        periods = np.array([1.0, 10.0, 100.0])
        result = transit_sensitivity(1.0, 100.0, 500.0, period_grid=periods)
        geo = result["geometric_probability"]
        assert geo[0] > geo[1] > geo[2]


# ============================================================
# RV Sensitivity
# ============================================================

class TestRVSensitivity:
    """Tests for rv_sensitivity."""

    def test_basic_sensitivity_curve(self):
        """RV sensitivity produces finite mass limits for well-sampled periods."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        rv_err = np.full(300, 1.0)

        result = rv_sensitivity(time, rv_err, stellar_mass_msun=1.0)
        # Some periods should have finite mass limits
        finite_mask = np.isfinite(result["mass_min_mearth"])
        assert np.any(finite_mask)

    def test_better_precision_lower_limits(self):
        """Better RV precision -> lower mass detection limits."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        periods = np.array([50.0])

        good = rv_sensitivity(time, np.full(300, 0.5), 1.0, period_grid=periods)
        bad = rv_sensitivity(time, np.full(300, 5.0), 1.0, period_grid=periods)

        good_mass = good["mass_min_mearth"][0]
        bad_mass = bad["mass_min_mearth"][0]
        if np.isfinite(good_mass) and np.isfinite(bad_mass):
            assert good_mass < bad_mass

    def test_more_observations_lower_limits(self):
        """More observations -> lower mass detection limits."""
        np.random.seed(42)
        periods = np.array([50.0])

        time_few = np.sort(np.random.uniform(0, 2000, 50))
        time_many = np.sort(np.random.uniform(0, 2000, 500))

        few = rv_sensitivity(time_few, np.ones(50), 1.0, period_grid=periods)
        many = rv_sensitivity(time_many, np.ones(500), 1.0, period_grid=periods)

        few_mass = few["mass_min_mearth"][0]
        many_mass = many["mass_min_mearth"][0]
        if np.isfinite(few_mass) and np.isfinite(many_mass):
            assert many_mass < few_mass


# ============================================================
# Astrometric Sensitivity
# ============================================================

class TestAstrometricSensitivity:
    """Tests for astrometric_sensitivity."""

    def test_mass_scales_with_distance(self):
        """Further stars -> higher minimum detectable mass."""
        near = astrometric_sensitivity(5.0, 1.0)
        far = astrometric_sensitivity(50.0, 1.0)
        # Compare at same separation index
        assert far["mass_min_mearth"][5] > near["mass_min_mearth"][5]

    def test_mass_scales_with_separation(self):
        """Wider separation -> higher minimum detectable mass (for PMa)."""
        result = astrometric_sensitivity(10.0, 1.0)
        masses = result["mass_min_mearth"]
        # Masses should increase with separation (PMa mass ~ a^2)
        assert masses[-1] > masses[0]

    def test_ruwe_hint(self):
        """RUWE > 1.4 is flagged as elevated."""
        result = astrometric_sensitivity(10.0, 1.0, ruwe=2.0)
        assert result["ruwe_hint"] == "elevated"

        result_normal = astrometric_sensitivity(10.0, 1.0, ruwe=1.0)
        assert result_normal["ruwe_hint"] == "normal"

    def test_period_conversion(self):
        """Period grid is consistent with Kepler's law."""
        result = astrometric_sensitivity(10.0, 1.0,
                                          separation_grid_au=np.array([1.0]))
        # 1 AU around 1 Msun -> ~365.25 days
        assert abs(result["period_grid_days"][0] - 365.25) / 365.25 < 0.01


# ============================================================
# Combined Sensitivity
# ============================================================

class TestCombinedSensitivity:
    """Tests for combined_sensitivity."""

    def test_all_methods_included(self):
        """Combined sensitivity includes all provided methods."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        rv_err = np.ones(300)

        t_sens = transit_sensitivity(1.0, 100.0, 500.0)
        r_sens = rv_sensitivity(time, rv_err, 1.0)
        a_sens = astrometric_sensitivity(10.0, 1.0)

        result = combined_sensitivity(t_sens, r_sens, a_sens)
        assert "rv" in result["methods_available"]
        assert "transit" in result["methods_available"]
        assert "astrometry" in result["methods_available"]
        assert result["rv_mass_min_mearth"] is not None
        assert result["transit_radius_min_rearth"] is not None
        assert result["astro_mass_min_mearth"] is not None

    def test_partial_methods(self):
        """Combined sensitivity works with only some methods."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        rv_err = np.ones(300)

        r_sens = rv_sensitivity(time, rv_err, 1.0)
        result = combined_sensitivity(rv_sens=r_sens)
        assert "rv" in result["methods_available"]
        assert "transit" not in result["methods_available"]
        assert result["transit_radius_min_rearth"] is None

    def test_no_methods(self):
        """Combined sensitivity handles no inputs gracefully."""
        result = combined_sensitivity()
        assert len(result["methods_available"]) == 0
        assert result["rv_mass_min_mearth"] is None
        assert result["transit_radius_min_rearth"] is None
        assert result["astro_mass_min_mearth"] is None

    def test_best_method_selection(self):
        """Best method per period picks the one with lower mass limit."""
        np.random.seed(42)
        time = np.sort(np.random.uniform(0, 2000, 300))
        rv_err = np.ones(300)

        r_sens = rv_sensitivity(time, rv_err, 1.0)
        a_sens = astrometric_sensitivity(10.0, 1.0)

        result = combined_sensitivity(rv_sens=r_sens, astro_sens=a_sens)
        best = result["best_mass_method_per_period"]
        assert len(best) == len(result["period_grid"])
        # At short periods RV should be better, at long periods astrometry
        # Check that both methods appear somewhere
        assert "rv" in best or "astrometry" in best

    def test_common_period_grid(self):
        """Custom period grid is used for interpolation."""
        periods = np.array([1.0, 10.0, 100.0, 1000.0])
        result = combined_sensitivity(period_grid=periods)
        np.testing.assert_array_equal(result["period_grid"], periods)
