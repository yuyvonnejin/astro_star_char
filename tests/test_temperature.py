"""Tests for Module 2: temperature and luminosity computation."""

import pytest
import numpy as np
from src.temperature import compute_temperature_luminosity, _compute_bc_g


class TestEffectiveTemperature:
    """Test Teff from color-temperature polynomial."""

    def test_solar_type_dwarf(self):
        """A solar-type star with (BP-RP)_0 ~ 0.82 should yield Teff near 5772 K."""
        # Solve for (BP-RP)_0 that gives Teff ~ 5772 K:
        # theta = 5040/5772 ~ 0.8732
        # For solar feh=0: theta = 0.4929 + 0.5092*C - 0.0353*C^2
        # Numerically, C ~ 0.82 gives theta ~ 0.873
        star = {
            "bp_rp": 0.82,
            "phot_g_mean_mag": 5.0,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.0,
            "feh": 0.0,
            "logg": 4.44,
        }
        result = compute_temperature_luminosity(star, distance_pc=10.0)

        assert result["teff_flag"] == "ok"
        assert result["teff_K"] == pytest.approx(5772, abs=200)

    def test_cool_dwarf(self):
        """A cool dwarf with (BP-RP)_0 = 1.5 should have Teff < 4500 K."""
        star = {
            "bp_rp": 1.5,
            "phot_g_mean_mag": 10.0,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.0,
            "feh": 0.0,
            "logg": 4.5,
        }
        result = compute_temperature_luminosity(star, distance_pc=50.0)

        assert result["teff_K"] < 4500
        assert result["teff_flag"] == "ok"

    def test_outside_color_range_flagged(self):
        """Very red star (BP-RP)_0 > 1.5 should be flagged for dwarfs."""
        star = {
            "bp_rp": 3.5,
            "phot_g_mean_mag": 11.13,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.005,
            "feh": 0.0,
            "logg": 4.5,
        }
        result = compute_temperature_luminosity(star, distance_pc=1.3)

        assert result["teff_flag"] == "outside_valid_range"

    def test_giant_coefficients_used(self):
        """A star with logg < 3.0 should use giant coefficients."""
        star = {
            "bp_rp": 1.0,
            "phot_g_mean_mag": 3.0,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.0,
            "feh": 0.0,
            "logg": 2.0,
        }
        result_giant = compute_temperature_luminosity(star, distance_pc=100.0)

        star["logg"] = 4.0
        result_dwarf = compute_temperature_luminosity(star, distance_pc=100.0)

        # Giant and dwarf should give somewhat different Teff for same color
        assert result_giant["teff_K"] != result_dwarf["teff_K"]


class TestBolometricCorrection:
    """Test BC_G computation."""

    def test_solar_teff(self):
        """BC_G for solar Teff (5772 K) should be approximately 0.06 (the a0 coefficient)."""
        bc = _compute_bc_g(5772.0)
        assert bc == pytest.approx(0.06, abs=0.001)

    def test_cool_star_bc(self):
        """BC_G for Teff=3500 K should use the cool-star coefficients."""
        bc = _compute_bc_g(3500.0)
        assert bc is not None
        # Cool star polynomial evaluated at 3500K (dt = -2272)
        assert bc != _compute_bc_g(5772.0)  # different from solar value

    def test_outside_range_returns_none(self):
        """BC_G for Teff > 8000 K or < 3300 K should return None."""
        assert _compute_bc_g(9000.0) is None
        assert _compute_bc_g(3000.0) is None


class TestLuminosity:
    """Test luminosity computation end-to-end."""

    def test_sun_at_10pc(self):
        """Sun-like star at 10 pc: G~4.83, distance=10 -> M_G~4.83, L~1 Lsun."""
        star = {
            "bp_rp": 0.82,
            "phot_g_mean_mag": 4.83,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.0,
            "feh": 0.0,
            "logg": 4.44,
            "lum_gspphot": None,
        }
        result = compute_temperature_luminosity(star, distance_pc=10.0)

        # M_G should be ~4.83 for a star at 10 pc with apparent mag 4.83
        assert result["M_G"] == pytest.approx(4.83, abs=0.01)
        # Luminosity should be roughly solar
        assert result["luminosity_Lsun"] == pytest.approx(1.0, rel=0.3)

    def test_validation_ratio_computed(self):
        """When lum_gspphot is provided, the validation ratio should be computed."""
        star = {
            "bp_rp": 0.82,
            "phot_g_mean_mag": 4.83,
            "ag_gspphot": 0.0,
            "ebpminrp_gspphot": 0.0,
            "feh": 0.0,
            "logg": 4.44,
            "lum_gspphot": 1.0,
        }
        result = compute_temperature_luminosity(star, distance_pc=10.0)

        assert result["luminosity_validation_ratio"] is not None
        assert 0.5 < result["luminosity_validation_ratio"] < 2.0
