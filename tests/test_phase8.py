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

    def test_hd10700_zero_confirmed(self):
        """Tau Ceti curation: 4 candidates, none confirmed (Figueira+2025)."""
        from src.reference_data import (load_known_planets,
                                        load_target_config)

        planets = load_known_planets("HD 10700")
        assert len(planets) == 4
        assert all(p.get("confirmed") is False for p in planets)
        periods = sorted(p["period_days"] for p in planets)
        assert periods == pytest.approx([20.0, 49.41, 162.87, 636.13])
        cfg = load_target_config("HD 10700")
        assert cfg.get("rotation_days") == 46
        assert cfg.get("activity_cycle_days") == 4000


# ============================================================
# Survey driver (8.3)
# ============================================================

class TestVetCandidatePeriod:
    """Candidate vetoes: alias families, rotation, cycle beats."""

    def test_one_day_alias_family(self):
        from src.survey import vet_candidate_period

        for p in (0.997, 1.0, 0.5, 2.0):
            ok, reason = vet_candidate_period(p)
            assert not ok, f"{p} d should be vetoed"
            assert "alias" in reason

    def test_annual_alias(self):
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(370.0)
        assert not ok and "annual" in reason
        ok, reason = vet_candidate_period(180.0)
        assert not ok and "annual" in reason

    def test_hd20794_414d_beat_vetoed(self):
        """The 414 d cycle-annual beat (Phase 7b lesson) is auto-vetoed."""
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(
            414.0, activity_cycle_days=3000.0)
        assert not ok
        assert "beat" in reason

    def test_rotation_harmonics_vetoed(self):
        from src.survey import vet_candidate_period

        for p in (39.0, 19.5, 78.0):
            ok, reason = vet_candidate_period(p, rotation_days=39.0)
            assert not ok and "rotation" in reason

    def test_known_planet_not_candidate(self):
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(
            89.5, known_periods=[18.315, 89.68, 647.6])
        assert not ok and "known planet" in reason

    def test_daily_alias_of_known_planet(self):
        """61 Vir b (4.215 d) daily-aliases to 1.31 d; auto-vetoed."""
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(
            1.31, known_periods=[4.21498, 38.079, 123.2])
        assert not ok and "daily alias" in reason
        # The other alias branch: 1/(1 + 1/4.215) = 0.808 d
        ok, reason = vet_candidate_period(
            0.81, known_periods=[4.21498])
        assert not ok and "daily alias" in reason

    def test_insufficient_coverage(self):
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(3000.0, baseline_days=4000.0)
        assert not ok and "coverage" in reason

    def test_genuine_candidate_passes(self):
        from src.survey import vet_candidate_period

        ok, reason = vet_candidate_period(
            55.0, known_periods=[18.315, 89.68],
            rotation_days=39.0, activity_cycle_days=3000.0,
            baseline_days=7000.0)
        assert ok and reason == "pass"


class TestBuildScorecard:
    """Scorecard construction from canned chain results (no network)."""

    @staticmethod
    def _target():
        return {"name": "Test Star", "hd": "HD 99999",
                "distance_pc": 5.0, "spectral_type": "G8V"}

    @staticmethod
    def _planets_ref():
        return [
            {"pl_name": "Test b", "period_days": 20.0, "k_ms": 1.0,
             "confirmed": True, "reference": "Test+2026"},
            {"pl_name": "Test c (candidate)", "period_days": 55.0,
             "k_ms": 0.5, "confirmed": False, "reference": "Test+2020"},
        ]

    @staticmethod
    def _rv_result(with_fit=True):
        n = 100
        rng = np.random.default_rng(11)
        time = np.sort(rng.uniform(55000, 58000, n))
        rv_raw = {
            "time": time,
            "rv": rng.normal(0, 2.0, n),
            "rv_err": np.full(n, 0.5),
            "instrument": np.array(["HARPS"] * n),
        }
        rv_data = {
            "status": "ok",
            "n_measurements": n, "n_measurements_raw": 3 * n,
            "n_measurements_filtered": 2 * n,
            "time_baseline_days": 3000.0,
            "instruments": ["HARPS"],
            "excluded_instruments": [],
            "periodogram": {"best_period": 55.2, "fap": 1e-6,
                            "peaks": [{"period_days": 55.2, "power": 0.4},
                                      {"period_days": 1.0, "power": 0.3}]},
        }
        residual = None
        if with_fit:
            residual = {
                "method": "keplerian",
                "keplerian_planets": [
                    {"period": 20.05, "k": 0.9, "k_err": 0.1,
                     "period_err": 0.01, "e": 0.05}],
                "keplerian_rms_before_ms": 5.0,
                "keplerian_rms_after_ms": 1.5,
                "residual_peaks": [
                    {"period_days": 414.0, "power": 0.2},
                    {"period_days": 120.0, "power": 0.15},
                ],
                "activity_decorrelation": {"indicators_used": ["fwhm"]},
                "gp_activity_model": {"hyperparameters": {}},
            }
        return {"rv_data": rv_data, "residual_analysis": residual,
                "rv_raw": rv_raw}

    def test_recovered_planet_matched(self):
        from src.survey import build_scorecard

        card = build_scorecard(
            self._target(), self._planets_ref(),
            {"activity_cycle_days": 3000, "rotation_days": 39},
            self._rv_result())
        assert card["status"] == "ok"
        assert card["n_planets_confirmed_ref"] == 1
        assert card["n_planets_recovered"] == 1
        rec = card["planets"][0]
        assert rec["recovered"] is True
        assert rec["k_ms_fit"] == pytest.approx(0.9)
        assert rec["k_n_sigma_vs_ref"] == pytest.approx(1.0)

    def test_beat_vetoed_genuine_survives(self):
        from src.survey import build_scorecard

        card = build_scorecard(
            self._target(), self._planets_ref(),
            {"activity_cycle_days": 3000, "rotation_days": 39},
            self._rv_result())
        vetoed = card["candidates"]["vetoed"]
        surviving = card["candidates"]["surviving"]
        assert any(v["period_days"] == 414.0 for v in vetoed)
        assert any(c["period_days"] == 120.0 for c in surviving)

    def test_no_fit_target_uses_raw_peaks(self):
        """Zero-confirmed-planet target: candidates from binned RV."""
        from src.survey import build_scorecard

        refs = [dict(self._planets_ref()[1])]  # only the unconfirmed one
        card = build_scorecard(
            self._target(), refs, {}, self._rv_result(with_fit=False))
        assert card["status"] == "ok"
        assert card["n_planets_confirmed_ref"] == 0
        assert card["candidates"]["peak_source"] == "binned_rv_no_fit"
        surviving = card["candidates"]["surviving"]
        # 55.2 d peak survives and cross-matches the literature candidate
        match = [c for c in surviving
                 if c["period_days"] == pytest.approx(55.2)]
        assert match
        assert match[0].get("matches_literature_candidate") \
            == "Test c (candidate)"
        # 1.0 d alias vetoed
        assert any(v["period_days"] == 1.0
                   for v in card["candidates"]["vetoed"])
        assert card["detection_limit"] is not None

    def test_detection_limit_summary(self):
        from src.survey import build_scorecard

        card = build_scorecard(
            self._target(), self._planets_ref(), {}, self._rv_result())
        limit = card["detection_limit"]
        assert limit["sigma_eff_ms"] == pytest.approx(1.5)
        k100 = limit["k_min_ms_at_period"]["100_d"]
        assert k100 is not None and 0 < k100 < 5

    def test_detection_limit_unconstrained_beyond_baseline(self):
        from src.survey import _detection_limit_summary

        rng = np.random.default_rng(4)
        time = np.sort(rng.uniform(55000, 56000, 80))  # 1000 d baseline
        limit = _detection_limit_summary(time, 2.0, 1000.0)
        assert limit["k_min_ms_at_period"]["1000_d"] is None
        assert limit["k_min_ms_at_period"]["100_d"] is not None


class TestRunSurveyRobustness:
    """One target failure never kills the survey."""

    def test_failure_injection(self, monkeypatch, tmp_path):
        import src.survey as survey

        def fake_survey_target(name, **kwargs):
            if name == "HD 20794":
                raise RuntimeError("injected failure")
            return {"target": {"name": name, "hd": name,
                               "distance_pc": 5.0},
                    "status": "ok", "data": {"status": "ok"},
                    "planets": [], "n_planets_confirmed_ref": 0,
                    "n_planets_recovered": 0,
                    "candidates": {"surviving": [], "vetoed": []},
                    "detection_limit": None, "residual_rms_ms": None}

        monkeypatch.setattr(survey, "survey_target", fake_survey_target)
        result = survey.run_survey(
            ["HD 20794", "HD 115617", "HD 10700"],
            out_dir=str(tmp_path), save=True)

        assert result["n_targets"] == 3
        assert result["n_ok"] == 2
        assert result["n_failed"] == 1
        failed = [c for c in result["targets"]
                  if c.get("status") == "error"]
        assert failed and "injected failure" in failed[0]["error"]
        # Outputs written despite the failure
        assert (tmp_path / "survey_results.json").exists()
        assert (tmp_path / "survey_summary.md").exists()

    def test_distance_ordering(self, monkeypatch, tmp_path):
        """Targets run nearest-first (Tau Ceti 3.65 pc before 61 Vir)."""
        import src.survey as survey

        order = []

        def fake_survey_target(name, **kwargs):
            order.append(name)
            return {"target": {"name": name}, "status": "ok",
                    "data": {"status": "ok"}, "planets": [],
                    "n_planets_confirmed_ref": 0,
                    "n_planets_recovered": 0,
                    "candidates": {"surviving": [], "vetoed": []},
                    "detection_limit": None, "residual_rms_ms": None}

        monkeypatch.setattr(survey, "survey_target", fake_survey_target)
        survey.run_survey(["HD 115617", "HD 10700", "HD 20794"],
                          save=False)
        assert order == ["HD 10700", "HD 20794", "HD 115617"]
