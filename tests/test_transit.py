"""Unit tests for Module 5: Transit detection and planet characterization.

Tests marked with @pytest.mark.network require MAST network access.
Run with: pytest -m network tests/test_transit.py
Skip with: pytest -m "not network"
"""

import math

import numpy as np
import pytest

from src.transit import (
    classify_planet_size,
    classify_transit_shape,
    compute_equilibrium_temp,
    compute_habitable_zone,
    compute_hz_period_range,
    compute_planet_properties,
    detect_transit,
    validate_even_odd,
    R_SUN_REARTH,
    R_SUN_AU,
    R_EARTH_RJUP,
)


# -- Offline tests (no network) -----------------------------------------------

class TestClassifyPlanetSize:
    """Test planet size classification (Fulton et al. 2017 categories)."""

    def test_sub_earth(self):
        assert classify_planet_size(0.5) == "sub_earth"

    def test_earth_sized(self):
        assert classify_planet_size(1.0) == "earth_sized"

    def test_super_earth(self):
        assert classify_planet_size(1.5) == "super_earth"

    def test_sub_neptune(self):
        assert classify_planet_size(3.0) == "sub_neptune"

    def test_neptune_sized(self):
        assert classify_planet_size(5.0) == "neptune_sized"

    def test_sub_jupiter(self):
        assert classify_planet_size(8.0) == "sub_jupiter"

    def test_jupiter_sized(self):
        # Jupiter ~ 11.2 R_earth
        assert classify_planet_size(11.2) == "jupiter_sized"

    def test_super_jupiter(self):
        assert classify_planet_size(20.0) == "super_jupiter"

    def test_boundary_earth_lower(self):
        assert classify_planet_size(0.8) == "earth_sized"

    def test_boundary_super_earth_lower(self):
        assert classify_planet_size(1.25) == "super_earth"


class TestComputeHabitableZone:
    """Test habitable zone boundaries (Kopparapu et al. 2013)."""

    def test_sun_hz_conservative(self):
        """Sun-like star: conservative HZ should be roughly 0.95-1.67 AU."""
        hz = compute_habitable_zone(1.0)
        assert 0.90 < hz["hz_conservative_inner_AU"] < 1.00
        assert 1.60 < hz["hz_conservative_outer_AU"] < 1.75

    def test_sun_hz_optimistic(self):
        """Sun-like star: optimistic HZ should be roughly 0.75-1.77 AU."""
        hz = compute_habitable_zone(1.0)
        assert 0.70 < hz["hz_optimistic_inner_AU"] < 0.80
        assert 1.70 < hz["hz_optimistic_outer_AU"] < 1.80

    def test_luminous_star_wider_hz(self):
        """More luminous star should have wider, more distant HZ."""
        hz_sun = compute_habitable_zone(1.0)
        hz_bright = compute_habitable_zone(4.0)
        # 4x luminosity -> boundaries scale as sqrt(L)
        assert hz_bright["hz_conservative_inner_AU"] > hz_sun["hz_conservative_inner_AU"]
        assert hz_bright["hz_conservative_outer_AU"] > hz_sun["hz_conservative_outer_AU"]

    def test_dim_star_narrower_hz(self):
        """Dim star should have closer-in HZ."""
        hz = compute_habitable_zone(0.01)  # 1% solar luminosity
        assert hz["hz_conservative_inner_AU"] < 0.15
        assert hz["hz_conservative_outer_AU"] < 0.20


class TestComputeEquilibriumTemp:
    """Test equilibrium temperature calculation."""

    def test_earth_like(self):
        """Earth at 1 AU from Sun should be ~255 K with albedo 0.3."""
        # T_star=5772 K, R_star=1 Rsun, a=1 AU
        T_eq = compute_equilibrium_temp(5772, 1.0, 1.0, albedo=0.3)
        assert 245 < T_eq < 265

    def test_hot_jupiter(self):
        """Hot Jupiter at 0.05 AU should be very hot."""
        T_eq = compute_equilibrium_temp(5772, 1.0, 0.05, albedo=0.1)
        assert T_eq > 1000

    def test_zero_albedo_hotter(self):
        """Zero albedo should give higher temperature."""
        T_eq_0 = compute_equilibrium_temp(5772, 1.0, 1.0, albedo=0.0)
        T_eq_3 = compute_equilibrium_temp(5772, 1.0, 1.0, albedo=0.3)
        assert T_eq_0 > T_eq_3

    def test_farther_is_cooler(self):
        """Planet at 2 AU should be cooler than at 1 AU."""
        T_eq_1 = compute_equilibrium_temp(5772, 1.0, 1.0)
        T_eq_2 = compute_equilibrium_temp(5772, 1.0, 2.0)
        assert T_eq_2 < T_eq_1


class TestComputePlanetProperties:
    """Test planet property derivation from transit + stellar parameters."""

    def _solar_stellar_props(self):
        return {
            "radius_Rsun": 1.0,
            "mass_Msun": 1.0,
            "teff_K": 5772,
            "luminosity_Lsun": 1.0,
        }

    def test_planet_radius_from_depth(self):
        """Known depth + solar-type star -> correct planet radius."""
        # Transit depth = (Rp/R*)^2
        # If Rp = 1 R_earth, R* = 1 R_sun: depth = (1/109.076)^2 = 8.41e-5
        rp_rstar = 1.0 / R_SUN_REARTH
        depth = rp_rstar ** 2

        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 365.25,  # 1 year
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert abs(result["planet_radius_Rearth"] - 1.0) < 0.1

    def test_keplers_third_law_earth(self):
        """Earth-like orbit: P=365.25 d, M*=1 Msun -> a ~ 1 AU."""
        depth = (1.0 / R_SUN_REARTH) ** 2
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 365.25,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert abs(result["orbital_semi_major_axis_AU"] - 1.0) < 0.01

    def test_keplers_third_law_hot_jupiter(self):
        """Hot Jupiter: P=3 d, M*=1 Msun -> a ~ 0.04 AU."""
        depth = (11.2 / R_SUN_REARTH) ** 2  # Jupiter-sized
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 3.0,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert 0.03 < result["orbital_semi_major_axis_AU"] < 0.05

    def test_habitable_zone_check(self):
        """Earth at 1 AU from Sun should be in habitable zone."""
        depth = (1.0 / R_SUN_REARTH) ** 2
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 365.25,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert result["in_habitable_zone"] is True

    def test_hot_planet_not_in_hz(self):
        """Planet at 0.05 AU should not be in habitable zone."""
        depth = (11.2 / R_SUN_REARTH) ** 2
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 3.0,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert result["in_habitable_zone"] is False

    def test_insolation_earth(self):
        """Earth insolation should be ~1 S_earth."""
        depth = (1.0 / R_SUN_REARTH) ** 2
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 365.25,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert 0.9 < result["insolation_Searth"] < 1.1

    def test_no_transit_detected(self):
        """When transit_detected is False, should return planet_flag='no_transit'."""
        transit = {"transit_detected": False}
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert result["planet_flag"] == "no_transit"

    def test_missing_stellar_params(self):
        """Missing stellar params -> limited output with flag."""
        depth = 0.001
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 5.0,
        }
        stellar = {
            "radius_Rsun": None,
            "mass_Msun": None,
            "teff_K": 5000,
            "luminosity_Lsun": 0.5,
        }
        result = compute_planet_properties(transit, stellar)
        assert result["planet_flag"] == "stellar_params_missing"
        assert result["planet_radius_ratio"] is not None  # ratio always computable
        assert result["planet_radius_Rearth"] is None

    def test_planet_flag_ok(self):
        """Full stellar params -> flag should be 'ok'."""
        depth = 0.001
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 5.0,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert result["planet_flag"] == "ok"

    def test_size_classification_in_result(self):
        """Result should include planet_size_class."""
        depth = (3.0 / R_SUN_REARTH) ** 2  # 3 R_earth sub-Neptune
        transit = {
            "transit_detected": True,
            "transit_depth": depth,
            "transit_period_days": 10.0,
        }
        result = compute_planet_properties(transit, self._solar_stellar_props())
        assert result["planet_size_class"] == "sub_neptune"


class TestDetectTransitSynthetic:
    """Test BLS detection on synthetic transit signals (offline)."""

    def _make_synthetic_transit(self, period=5.0, depth=0.01, duration_days=0.15,
                                 n_points=5000, baseline=90.0, noise=0.001):
        """Generate a synthetic light curve with box-shaped transit dips."""
        rng = np.random.default_rng(42)
        time = np.sort(rng.uniform(0, baseline, n_points))
        flux = np.ones(n_points) + rng.normal(0, noise, n_points)

        # Insert box transits
        phase = (time % period) / period
        in_transit = phase < (duration_days / period)
        flux[in_transit] -= depth

        flux_err = np.full(n_points, noise)
        return time, flux, flux_err

    def test_recovers_period(self):
        """BLS should recover the correct transit period from synthetic data."""
        period_true = 5.0
        time, flux, flux_err = self._make_synthetic_transit(
            period=period_true, depth=0.01, n_points=10000, baseline=100.0
        )
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        assert result["transit_detected"] is True
        # Allow 2% tolerance
        assert abs(result["transit_period_days"] - period_true) / period_true < 0.02

    def test_recovers_depth(self):
        """BLS should recover approximate transit depth."""
        depth_true = 0.01
        time, flux, flux_err = self._make_synthetic_transit(
            period=5.0, depth=depth_true, n_points=10000, baseline=100.0
        )
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        assert result["transit_detected"] is True
        # Depth within factor of 2
        assert 0.005 < result["transit_depth"] < 0.02

    def test_no_transit_in_noise(self):
        """Pure noise should produce low SDE (below typical detection thresholds).

        Note: BLS can occasionally find marginal 'signals' in pure noise at low
        SDE. We use a high threshold to ensure robust non-detection.
        """
        rng = np.random.default_rng(123)
        time = np.sort(rng.uniform(0, 90, 5000))
        flux = 1.0 + rng.normal(0, 0.001, 5000)
        flux_err = np.full(5000, 0.001)

        # Use default SDE threshold (6.0) and check that any 'detection' has low SDE
        result = detect_transit(time, flux, flux_err, sde_threshold=6.0)
        if result["transit_detected"]:
            # If BLS did trigger, SDE should be marginal (not a strong signal)
            assert result["transit_sde"] < 15, "Noise should not produce high SDE"
        # With a stricter threshold, should not detect
        result_strict = detect_transit(time, flux, flux_err, sde_threshold=15.0)
        assert result_strict["transit_detected"] is False

    def test_short_data_returns_gracefully(self):
        """Very short data should return no detection, not crash."""
        time = np.linspace(0, 1, 30)
        flux = np.ones(30)
        flux_err = np.full(30, 0.001)

        result = detect_transit(time, flux, flux_err)
        assert result["transit_detected"] is False

    def test_deep_transit_high_sde(self):
        """Deep transit should produce high SDE."""
        time, flux, flux_err = self._make_synthetic_transit(
            period=3.0, depth=0.05, n_points=10000, baseline=60.0, noise=0.001
        )
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        assert result["transit_detected"] is True
        assert result["transit_sde"] > 10

    def test_result_has_expected_keys(self):
        """Detection result should have all expected keys."""
        time, flux, flux_err = self._make_synthetic_transit()
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        assert "transit_detected" in result
        assert "detection_method" in result
        assert result["detection_method"] == "bls"
        if result["transit_detected"]:
            assert "transit_period_days" in result
            assert "transit_depth" in result
            assert "transit_depth_ppm" in result
            assert "transit_duration_hours" in result
            assert "transit_sde" in result
            assert "n_transits_observed" in result

    def test_candidates_list_present(self):
        """Detection result should include transit_candidates list."""
        time, flux, flux_err = self._make_synthetic_transit(
            period=5.0, depth=0.01, n_points=10000, baseline=100.0
        )
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        assert result["transit_detected"] is True
        assert "transit_candidates" in result
        candidates = result["transit_candidates"]
        assert isinstance(candidates, list)
        assert len(candidates) >= 1
        # Best candidate should match top-level result
        assert candidates[0]["rank"] == 1
        assert candidates[0]["transit_period_days"] == result["transit_period_days"]

    def test_candidates_span_period_range(self):
        """Stratified extraction should produce candidates across the period range."""
        time, flux, flux_err = self._make_synthetic_transit(
            period=5.0, depth=0.01, n_points=10000, baseline=100.0
        )
        result = detect_transit(time, flux, flux_err, sde_threshold=5.0)
        candidates = result.get("transit_candidates", [])
        periods = [c["transit_period_days"] for c in candidates]
        if len(periods) >= 2:
            # Candidates should span more than a single period decade
            p_range = max(periods) / min(periods)
            assert p_range > 2.0, f"Candidates too clustered: range ratio {p_range:.1f}"


class TestComputeHZPeriodRange:
    """Test HZ period range computation."""

    def test_sun_like_range(self):
        """Sun-like star: HZ period range should bracket ~365 days."""
        result = compute_hz_period_range(1.0, 1.0, broadening_factor=2.0)
        assert result is not None
        min_p, max_p = result
        # Earth's period (365.25 d) should be within the range
        assert min_p < 365.25
        assert max_p > 365.25

    def test_m_dwarf_shorter(self):
        """M-dwarf: HZ periods should be much shorter than Sun-like."""
        result_sun = compute_hz_period_range(1.0, 1.0)
        result_mdwarf = compute_hz_period_range(0.12, 0.0017)
        assert result_sun is not None
        assert result_mdwarf is not None
        assert result_mdwarf[1] < result_sun[0]

    def test_brighter_wider(self):
        """Brighter star should have wider (longer period) HZ."""
        result_sun = compute_hz_period_range(1.0, 1.0, broadening_factor=1.0)
        result_bright = compute_hz_period_range(2.0, 10.0, broadening_factor=1.0)
        assert result_sun is not None
        assert result_bright is not None
        assert result_bright[1] > result_sun[1]

    def test_missing_mass_returns_none(self):
        """Missing mass should return None."""
        assert compute_hz_period_range(None, 1.0) is None

    def test_missing_luminosity_returns_none(self):
        """Missing luminosity should return None."""
        assert compute_hz_period_range(1.0, None) is None

    def test_min_floor(self):
        """Minimum period should be floored at 0.5 days."""
        # Very bright, massive star -- inner HZ very close -> short period
        result = compute_hz_period_range(1.0, 1.0, broadening_factor=1.0)
        assert result is not None
        assert result[0] >= 0.5

    def test_no_broadening(self):
        """With broadening_factor=1.0, range should be narrower."""
        result_broad = compute_hz_period_range(1.0, 1.0, broadening_factor=2.0)
        result_narrow = compute_hz_period_range(1.0, 1.0, broadening_factor=1.0)
        assert result_broad is not None
        assert result_narrow is not None
        assert result_narrow[0] >= result_broad[0]
        assert result_narrow[1] <= result_broad[1]


class TestValidateEvenOdd:
    """Test even/odd transit validation."""

    def _make_transit_lc(self, period=5.0, depth=0.01, duration_days=0.15,
                          n_points=10000, baseline=100.0, noise=0.0005,
                          depth_ratio=1.0):
        """Make synthetic LC with optional even/odd depth difference."""
        rng = np.random.default_rng(42)
        time = np.sort(rng.uniform(0, baseline, n_points))
        flux = np.ones(n_points) + rng.normal(0, noise, n_points)

        epoch_num = np.round(time / period).astype(int)
        phase = time - epoch_num * period
        in_transit = np.abs(phase) < duration_days / 2.0

        is_even = (epoch_num % 2) == 0
        flux[in_transit & is_even] -= depth
        flux[in_transit & ~is_even] -= depth * depth_ratio

        return time, flux

    def test_consistent_depths_pass(self):
        """Equal even/odd depths should pass validation."""
        time, flux = self._make_transit_lc(depth_ratio=1.0)
        result = validate_even_odd(time, flux, 5.0, 0.0, 0.15)
        assert result["even_odd_validation_pass"] is True
        assert result["even_odd_flag"] == "ok"

    def test_inconsistent_depths_fail(self):
        """Very different even/odd depths should fail (possible EB)."""
        time, flux = self._make_transit_lc(depth=0.01, depth_ratio=5.0)
        result = validate_even_odd(time, flux, 5.0, 0.0, 0.15)
        assert result["even_odd_validation_pass"] is False
        assert result["even_odd_flag"] == "even_odd_depth_mismatch"

    def test_too_few_points(self):
        """Very short data should report too few points."""
        time = np.linspace(0, 2, 30)
        flux = np.ones(30)
        result = validate_even_odd(time, flux, 5.0, 0.0, 0.15)
        assert result["even_odd_validation_pass"] is None
        assert "too_few" in result["even_odd_flag"] or "insufficient" in result["even_odd_flag"]

    def test_invalid_period(self):
        """Invalid period should return gracefully."""
        time = np.linspace(0, 100, 1000)
        flux = np.ones(1000)
        result = validate_even_odd(time, flux, 0, 0.0, 0.15)
        assert result["even_odd_flag"] == "invalid_parameters"

    def test_result_keys(self):
        """Result should have all expected keys."""
        time, flux = self._make_transit_lc()
        result = validate_even_odd(time, flux, 5.0, 0.0, 0.15)
        expected_keys = {"depth_even_ppm", "depth_odd_ppm", "depth_ratio_even_odd",
                         "even_odd_validation_pass", "n_even", "n_odd", "even_odd_flag"}
        assert expected_keys.issubset(result.keys())


class TestClassifyTransitShape:
    """Test transit shape classification (V vs U)."""

    def _make_box_transit(self, period=5.0, depth=0.01, duration_days=0.2,
                           n_points=20000, baseline=200.0, noise=0.0003):
        """Make synthetic box transit (U-shaped)."""
        rng = np.random.default_rng(42)
        time = np.sort(rng.uniform(0, baseline, n_points))
        flux = np.ones(n_points) + rng.normal(0, noise, n_points)

        phase = (time % period) / period
        in_transit = phase < (duration_days / period)
        flux[in_transit] -= depth

        return time, flux

    def _make_v_transit(self, period=5.0, depth=0.01, duration_days=0.2,
                         n_points=20000, baseline=200.0, noise=0.0003):
        """Make synthetic V-shaped (triangular) transit."""
        rng = np.random.default_rng(42)
        time = np.sort(rng.uniform(0, baseline, n_points))
        flux = np.ones(n_points) + rng.normal(0, noise, n_points)

        phase = (time % period) / period
        dur_phase = duration_days / period
        half_dur = dur_phase / 2.0

        # Triangular shape: linear ingress and egress, no flat bottom
        in_first_half = phase < half_dur
        in_second_half = (phase >= half_dur) & (phase < dur_phase)
        flux[in_first_half] -= depth * (phase[in_first_half] / half_dur)
        # Flip: going from max depth back to 0
        flux[in_second_half] -= depth * (1.0 - (phase[in_second_half] - half_dur) / half_dur)

        return time, flux

    def test_u_shape_detected(self):
        """Box transit should be classified as U-shaped."""
        time, flux = self._make_box_transit()
        result = classify_transit_shape(time, flux, 5.0, 0.0, 0.2)
        assert result["shape_class"] == "U_shape"
        assert result["shape_flag"] == "ok"

    def test_v_shape_detected(self):
        """Triangular transit should be classified as V-shaped."""
        time, flux = self._make_v_transit()
        result = classify_transit_shape(time, flux, 5.0, 0.0, 0.2)
        assert result["shape_class"] == "V_shape"
        assert result["shape_flag"] == "ok"

    def test_too_few_points(self):
        """Very few points should report too_few_points."""
        time = np.linspace(0, 2, 15)
        flux = np.ones(15)
        result = classify_transit_shape(time, flux, 5.0, 0.0, 0.2)
        assert result["shape_class"] is None
        assert "too_few" in result["shape_flag"]

    def test_invalid_params(self):
        """Invalid period/duration should return gracefully."""
        time = np.linspace(0, 100, 1000)
        flux = np.ones(1000)
        result = classify_transit_shape(time, flux, 0, 0.0, 0.2)
        assert result["shape_flag"] == "invalid_parameters"

    def test_result_keys(self):
        """Result should have all expected keys."""
        time, flux = self._make_box_transit()
        result = classify_transit_shape(time, flux, 5.0, 0.0, 0.2)
        expected_keys = {"shape_class", "flat_bottom_fraction", "n_in_transit", "shape_flag"}
        assert expected_keys.issubset(result.keys())


# -- Network tests (require MAST access) --------------------------------------

@pytest.mark.network
class TestTransitNetworkKIC6922244:
    """Test BLS transit detection on KIC 6922244 (Kepler-410A b)."""

    def test_bls_finds_transit(self):
        """BLS should detect transit signal at P ~ 3.52 d for KIC 6922244."""
        from src.lightcurve import retrieve_lightcurve

        lc = retrieve_lightcurve("KIC 6922244", mission="Kepler", max_sectors=4)
        assert lc is not None

        result = detect_transit(
            lc["time"], lc["flux_flat"], lc.get("flux_err"),
            min_period=0.5, max_period=50.0, sde_threshold=5.0,
        )
        # Should detect the transit
        assert result["transit_detected"] is True
        # Period should be near 3.5225 d (not the 0.88 harmonic that LS found)
        detected_p = result["transit_period_days"]
        assert 3.0 < detected_p < 4.0, f"Expected ~3.52 d, got {detected_p:.4f} d"
