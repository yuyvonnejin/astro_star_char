"""Tests for Phase 8: survey infrastructure.

8.1: MCMC parameter uncertainties.
Later steps (reference data files, survey driver, catalog shells)
add their tests here.
"""

import numpy as np
import pytest


# ============================================================
# MCMC uncertainties (8.1)
# ============================================================

class TestChainErr:
    """Tests for the percentile-based error helper."""

    def test_gaussian_chain(self):
        from src.rv_keplerian import _chain_err

        rng = np.random.default_rng(2)
        chains = {"k1": rng.normal(5.0, 0.3, 20000)}
        err = _chain_err(chains, "k1")
        assert err == pytest.approx(0.3, rel=0.05)

    def test_missing_key_returns_none(self):
        from src.rv_keplerian import _chain_err

        assert _chain_err({"k1": [1, 2, 3]}, "per1") is None
        assert _chain_err(None, "k1") is None


class TestMCMCUncertainties:
    """MCMC on a synthetic planet returns sensible error bars."""

    @pytest.fixture(scope="class")
    def mcmc_result(self):
        radvel = pytest.importorskip("radvel")
        from src.rv_keplerian import fit_keplerian

        rng = np.random.default_rng(7)
        n = 120
        t = np.sort(rng.uniform(50000, 51000, n))
        params = radvel.Parameters(1, basis='per tc e w k')
        params['per1'] = radvel.Parameter(value=20.0)
        params['tc1'] = radvel.Parameter(value=50500.0)
        params['e1'] = radvel.Parameter(value=0.0)
        params['w1'] = radvel.Parameter(value=0.0)
        params['k1'] = radvel.Parameter(value=5.0)
        params['dvdt'] = radvel.Parameter(value=0.0)
        params['curv'] = radvel.Parameter(value=0.0)
        mod = radvel.RVModel(params, time_base=50500.0)
        rv = mod(t) + rng.normal(0, 0.5, n)
        err = np.full(n, 0.5)

        return fit_keplerian(
            t, rv, err, ["HARPS"] * n,
            planet_params=[{"period": 20.0, "tc": 50500.0,
                            "e": 0.0, "w": 0.0, "k": 4.0}],
            run_mcmc=True, mcmc_nwalkers=30, mcmc_nsteps=400,
        )

    def test_errors_populated(self, mcmc_result):
        assert mcmc_result["status"] == "ok"
        p = mcmc_result["planets"][0]
        for key in ("k_err", "period_err", "tc_err"):
            assert p[key] is not None and p[key] > 0
        inst = mcmc_result["instruments"]["HARPS"]
        assert inst["gamma_err"] is not None

    def test_k_error_magnitude(self, mcmc_result):
        """K error near the analytic photon-limit estimate.

        sigma_K ~ noise * sqrt(2/N) = 0.5 * sqrt(2/120) = 0.065 m/s.
        Allow a factor 3 band; short chains are not fully converged.
        """
        k_err = mcmc_result["planets"][0]["k_err"]
        assert 0.02 < k_err < 0.20

    def test_k_consistent_with_truth(self, mcmc_result):
        p = mcmc_result["planets"][0]
        assert abs(abs(p["k"]) - 5.0) < 4 * p["k_err"]

    def test_no_mcmc_no_errors(self):
        """Without run_mcmc, error fields are None (not missing)."""
        pytest.importorskip("radvel")
        from src.rv_keplerian import fit_keplerian

        rng = np.random.default_rng(3)
        n = 80
        t = np.sort(rng.uniform(50000, 50500, n))
        rv = 5.0 * np.sin(2 * np.pi * t / 20.0) + rng.normal(0, 0.5, n)
        result = fit_keplerian(
            t, rv, np.full(n, 0.5), ["HARPS"] * n,
            planet_params=[{"period": 20.0, "tc": 50250.0,
                            "e": 0.0, "w": 0.0, "k": 4.0}],
        )
        assert result["status"] == "ok"
        p = result["planets"][0]
        assert p["k_err"] is None
        assert p["period_err"] is None


# ============================================================
# Reference data loader (8.2)
# ============================================================

class TestReferenceData:
    """data/known_planets.json loader."""

    def test_hd20794_four_planets_one_disputed(self):
        from src.reference_data import load_known_planets

        planets = load_known_planets("HD 20794")
        assert len(planets) == 4
        disputed = [p for p in planets if p.get("confirmed") is False]
        assert len(disputed) == 1
        assert disputed[0]["period_days"] == pytest.approx(147.025)

    def test_hd115617_three_confirmed(self):
        from src.reference_data import load_known_planets

        planets = load_known_planets("HD 115617")
        assert len(planets) == 3
        assert all(p.get("confirmed") for p in planets)
        assert all(p.get("k_ms") for p in planets)

    def test_unknown_target_none(self):
        from src.reference_data import load_known_planets

        assert load_known_planets("HD 99999") is None

    def test_config_lookup(self):
        from src.reference_data import load_target_config

        cfg = load_target_config("HD 115617")
        assert cfg.get("activity_cycle_days") == 5911
        assert load_target_config("HD 99999") == {}

    def test_fallback_wrapper_uses_loader(self):
        from src.deep_dive import _known_planets_fallback

        planets = _known_planets_fallback({"hd": "HD 20794"})
        assert planets is not None and len(planets) == 4
        assert _known_planets_fallback({"hd": None}) is None

    def test_all_planets_have_period(self):
        from src.reference_data import _load

        for hd, entry in _load().items():
            for p in entry.get("planets", []):
                assert p.get("period_days"), f"{hd} planet missing period"
