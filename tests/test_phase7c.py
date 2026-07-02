"""Tests for Phase 7c: RV convergence improvements.

7c.1: Nightly binning (rv_bin_nightly).
Later steps (activity decorrelation, GP) add their tests here.
"""

import numpy as np
import pytest

from src.rv_data import rv_bin_nightly


# ============================================================
# Nightly Binning
# ============================================================

class TestRVBinNightly:
    """Tests for rv_bin_nightly."""

    def _make_rv_data(self, time, rv, rv_err, instruments):
        return {
            "time": np.asarray(time, dtype=float),
            "rv": np.asarray(rv, dtype=float),
            "rv_err": np.asarray(rv_err, dtype=float),
            "instrument": list(instruments),
            "instruments": sorted(set(instruments)),
            "n_measurements": len(time),
            "time_baseline_days": float(time[-1] - time[0]) if len(time) > 1 else 0.0,
        }

    def test_bins_same_night_same_instrument(self):
        """Three points in one night collapse to one."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2, 50000.3],
            rv=[1.0, 2.0, 3.0],
            rv_err=[1.0, 1.0, 1.0],
            instruments=["HARPS"] * 3,
        )
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 1
        assert binned["n_measurements_unbinned"] == 3
        assert binned["binned_nightly"] is True
        # Equal errors: weighted mean = plain mean
        assert binned["rv"][0] == pytest.approx(2.0)

    def test_different_nights_not_merged(self):
        """Points on different nights stay separate."""
        data = self._make_rv_data(
            time=[50000.1, 50001.1, 50002.1],
            rv=[1.0, 2.0, 3.0],
            rv_err=[1.0, 1.0, 1.0],
            instruments=["HARPS"] * 3,
        )
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 3

    def test_different_instruments_not_merged(self):
        """Same night, different instruments -> separate bins."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2],
            rv=[1.0, 5.0],
            rv_err=[1.0, 1.0],
            instruments=["HARPS", "ESPRESSO"],
        )
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 2
        assert sorted(binned["instruments"]) == ["ESPRESSO", "HARPS"]

    def test_weighted_mean(self):
        """Lower-error point dominates the bin."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2],
            rv=[0.0, 10.0],
            rv_err=[0.1, 10.0],
            instruments=["HARPS"] * 2,
        )
        binned = rv_bin_nightly(data)
        # Weights 100^2 : 0.01^2 -> mean very close to 0.0
        assert abs(binned["rv"][0]) < 0.01

    def test_error_shrinks_with_n(self):
        """Binned error is below single-point error for consistent points."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2, 50000.3, 50000.4],
            rv=[1.0, 1.0, 1.0, 1.0],
            rv_err=[1.0, 1.0, 1.0, 1.0],
            instruments=["HARPS"] * 4,
        )
        binned = rv_bin_nightly(data)
        # Zero scatter -> propagated error = 1/sqrt(4) = 0.5
        assert binned["rv_err"][0] == pytest.approx(0.5)

    def test_scatter_inflates_error(self):
        """Intra-night scatter beyond formal errors inflates the bin error."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2],
            rv=[-5.0, 5.0],
            rv_err=[0.1, 0.1],
            instruments=["HARPS"] * 2,
        )
        binned = rv_bin_nightly(data)
        # Propagated error would be 0.07; scatter error = std/sqrt(2) ~ 5.0
        assert binned["rv_err"][0] > 1.0

    def test_output_sorted_by_time(self):
        """Output is time-sorted even with interleaved instruments."""
        data = self._make_rv_data(
            time=[50002.1, 50000.1, 50001.1],
            rv=[3.0, 1.0, 2.0],
            rv_err=[1.0, 1.0, 1.0],
            instruments=["A", "B", "A"],
        )
        binned = rv_bin_nightly(data)
        assert np.all(np.diff(binned["time"]) >= 0)

    def test_preserves_extra_keys(self):
        """Extra keys like excluded_instruments pass through."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2],
            rv=[1.0, 2.0],
            rv_err=[1.0, 1.0],
            instruments=["HARPS"] * 2,
        )
        data["excluded_instruments"] = ["CORAVEL-S"]
        binned = rv_bin_nightly(data)
        assert binned["excluded_instruments"] == ["CORAVEL-S"]

    def test_empty_input(self):
        """Empty data returns empty data, no crash."""
        data = self._make_rv_data(
            time=np.array([]), rv=np.array([]), rv_err=np.array([]),
            instruments=[],
        )
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 0
        assert binned["binned_nightly"] is True

    def test_instrument_summary_rebuilt(self):
        """Instrument summary reflects binned counts."""
        data = self._make_rv_data(
            time=[50000.1, 50000.2, 50001.1],
            rv=[1.0, 2.0, 3.0],
            rv_err=[1.0, 1.0, 1.0],
            instruments=["HARPS"] * 3,
        )
        binned = rv_bin_nightly(data)
        assert binned["instrument_summary"]["HARPS"]["n_measurements"] == 2

    def test_realistic_compression(self):
        """Dense HARPS-like sampling compresses roughly to nights."""
        rng = np.random.default_rng(42)
        # 200 nights, 1-8 obs per night at fractional times 0.55-0.95
        times = []
        for night in range(50000, 50200):
            n_obs = rng.integers(1, 9)
            times.extend(night + 0.55 + 0.4 * rng.random(n_obs))
        times = np.sort(np.array(times))
        n = len(times)
        data = self._make_rv_data(
            time=times,
            rv=rng.normal(0, 1.0, n),
            rv_err=np.full(n, 0.5),
            instruments=["HARPS"] * n,
        )
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 200
        assert binned["n_measurements_unbinned"] == n


# ============================================================
# Keplerian fit on binned data (integration, synthetic)
# ============================================================

class TestBinnedKeplerianRecovery:
    """Binning should not break planet recovery on synthetic data."""

    def test_recovery_after_binning(self):
        """Inject one planet, bin, fit, recover K within tolerance."""
        radvel = pytest.importorskip("radvel")
        from src.rv_keplerian import fit_keplerian

        rng = np.random.default_rng(1)
        # 300 nights, 2 obs per night
        base_nights = np.arange(51000, 51300)
        time = np.sort(np.concatenate([
            base_nights + 0.6 + 0.1 * rng.random(300),
            base_nights + 0.8 + 0.1 * rng.random(300),
        ]))
        period, k_true, tc_true = 25.0, 5.0, 51150.0

        params = radvel.Parameters(1, basis='per tc e w k')
        params['per1'] = radvel.Parameter(value=period)
        params['tc1'] = radvel.Parameter(value=tc_true)
        params['e1'] = radvel.Parameter(value=0.0)
        params['w1'] = radvel.Parameter(value=0.0)
        params['k1'] = radvel.Parameter(value=k_true)
        params['dvdt'] = radvel.Parameter(value=0.0)
        params['curv'] = radvel.Parameter(value=0.0)
        mod = radvel.RVModel(params, time_base=float(np.median(time)))

        rv = mod(time) + rng.normal(0, 0.5, len(time))
        rv_err = np.full(len(time), 0.5)

        data = {
            "time": time, "rv": rv, "rv_err": rv_err,
            "instrument": ["HARPS"] * len(time),
            "instruments": ["HARPS"],
            "n_measurements": len(time),
            "time_baseline_days": float(time[-1] - time[0]),
        }
        binned = rv_bin_nightly(data)
        assert binned["n_measurements"] == 300

        result = fit_keplerian(
            binned["time"], binned["rv"], binned["rv_err"],
            instruments=binned["instrument"],
            planet_params=[{
                "period": period, "tc": tc_true,
                "e": 0.0, "w": 0.0, "k": 3.0,
            }],
        )
        assert result["status"] == "ok"
        k_fit = abs(result["planets"][0]["k"])
        assert k_fit == pytest.approx(k_true, rel=0.15)


# ============================================================
# Vary-flag regression (RadVel vector cache)
# ============================================================

class TestVaryFlagSync:
    """RadVel >= 1.4 caches vary flags in a vector at Posterior
    construction; without an explicit sync, .vary changes are ignored
    and every parameter stays free. Regression tests for the fix.
    """

    def _fit(self, fix_e, e_inject, e_init):
        radvel = pytest.importorskip("radvel")
        from src.rv_keplerian import fit_keplerian

        rng = np.random.default_rng(7)
        time = np.sort(rng.uniform(50000, 51000, 200))
        period, k_true, tc_true = 20.0, 8.0, 50500.0

        params = radvel.Parameters(1, basis='per tc e w k')
        params['per1'] = radvel.Parameter(value=period)
        params['tc1'] = radvel.Parameter(value=tc_true)
        params['e1'] = radvel.Parameter(value=e_inject)
        params['w1'] = radvel.Parameter(value=0.5)
        params['k1'] = radvel.Parameter(value=k_true)
        params['dvdt'] = radvel.Parameter(value=0.0)
        params['curv'] = radvel.Parameter(value=0.0)
        mod = radvel.RVModel(params, time_base=50500.0)

        rv = mod(time) + rng.normal(0, 0.5, 200)
        rv_err = np.full(200, 0.5)

        return fit_keplerian(
            time, rv, rv_err, instruments=["HARPS"] * 200,
            planet_params=[{
                "period": period, "tc": tc_true,
                "e": e_init, "w": 0.5, "k": 5.0,
            }],
            fix_eccentricities=fix_e,
        )

    def test_fixed_e_stays_at_initial_value(self):
        """With fix_eccentricities=True, output e == input e exactly."""
        result = self._fit(fix_e=True, e_inject=0.3, e_init=0.3)
        assert result["status"] == "ok"
        assert result["planets"][0]["e"] == pytest.approx(0.3, abs=1e-10)

    def test_free_e_moves_toward_truth(self):
        """With fix_eccentricities=False, e is fitted (moves from init)."""
        result = self._fit(fix_e=False, e_inject=0.4, e_init=0.05)
        assert result["status"] == "ok"
        e_fit = result["planets"][0]["e"]
        assert abs(e_fit - 0.05) > 1e-6
        assert e_fit == pytest.approx(0.4, abs=0.15)
