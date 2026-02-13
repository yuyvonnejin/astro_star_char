"""Tests for Module 3: mass estimation."""

import pytest
from src.mass import compute_mass


class TestMassLuminosity:
    """Test piecewise mass-luminosity relation."""

    def test_solar_mass(self):
        """L = 1.0 Lsun should give M ~ 1.0 Msun (alpha=4.0 segment)."""
        star = {"logg": 4.44}
        result = compute_mass(star, teff_K=5772, luminosity_Lsun=1.0)

        assert result["is_main_sequence"] is True
        assert result["mass_flag"] == "ok"
        assert result["mass_Msun"] == pytest.approx(1.0, abs=0.01)

    def test_low_mass_segment(self):
        """L = 0.01 Lsun -> alpha=2.3, M = 0.01^(1/2.3) ~ 0.12 Msun."""
        star = {"logg": 5.0}
        result = compute_mass(star, teff_K=3500, luminosity_Lsun=0.01)

        assert result["mass_flag"] == "ok"
        expected = 0.01 ** (1.0 / 2.3)
        assert result["mass_Msun"] == pytest.approx(expected, rel=0.01)

    def test_intermediate_mass_segment(self):
        """L = 10 Lsun -> alpha=4.0, M = 10^(1/4.0) ~ 1.78 Msun."""
        star = {"logg": 4.0}
        result = compute_mass(star, teff_K=7000, luminosity_Lsun=10.0)

        assert result["mass_flag"] == "ok"
        expected = 10.0 ** (1.0 / 4.0)
        assert result["mass_Msun"] == pytest.approx(expected, rel=0.01)

    def test_high_mass_segment(self):
        """L = 1000 Lsun -> alpha=3.5, M = 1000^(1/3.5)."""
        star = {"logg": 4.0}
        result = compute_mass(star, teff_K=7000, luminosity_Lsun=1000.0)

        assert result["mass_flag"] == "ok"
        expected = 1000.0 ** (1.0 / 3.5)
        assert result["mass_Msun"] == pytest.approx(expected, rel=0.01)


class TestMainSequenceCheck:
    """Test main-sequence classification."""

    def test_not_main_sequence_low_logg(self):
        """logg < 3.5 -> not main-sequence."""
        star = {"logg": 2.0}
        result = compute_mass(star, teff_K=5000, luminosity_Lsun=100.0)

        assert result["is_main_sequence"] is False
        assert result["mass_flag"] == "not_main_sequence"
        assert result["mass_Msun"] is None

    def test_not_main_sequence_hot(self):
        """Teff > 8000 K -> not main-sequence for this relation."""
        star = {"logg": 4.0}
        result = compute_mass(star, teff_K=9000, luminosity_Lsun=25.0)

        assert result["is_main_sequence"] is False
        assert result["mass_Msun"] is None

    def test_not_main_sequence_cool(self):
        """Teff < 3300 K -> not main-sequence for this relation."""
        star = {"logg": 5.0}
        result = compute_mass(star, teff_K=3000, luminosity_Lsun=0.001)

        assert result["is_main_sequence"] is False

    def test_no_luminosity(self):
        """Main-sequence star but luminosity is None."""
        star = {"logg": 4.5}
        result = compute_mass(star, teff_K=5000, luminosity_Lsun=None)

        assert result["is_main_sequence"] is True
        assert result["mass_flag"] == "no_luminosity"
        assert result["mass_Msun"] is None
