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
    compute_equilibrium_temp,
    compute_habitable_zone,
    compute_planet_properties,
    detect_transit,
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
