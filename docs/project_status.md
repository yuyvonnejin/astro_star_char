# Project Status -- Read This First When Resuming

Last updated: 2026-07-03 (end of Phase 8.4)

## What this project is

Goal: meaningfully participate in, and if possible contribute to,
exoplanet search. Detection bias means only the most obvious planets
are findable, so the strategy is impact through proximity: survey the
nearest stars first, expanding outward in distance shells. The
deliverable is infrastructure -- a uniform RV analysis chain with
honest detection limits -- not single-target precision records.

Realistic contribution ladder:
1. Data quality findings (have: DACE corruption catalog)
2. Uniform survey with per-target detection limits (value even with
   zero discoveries)
3. Candidate flagging for community follow-up (not discovery claims)
4. Positioned for new data (Gaia DR4 late 2026, DACE/ESO releases)
   -- highest-probability path to a real find

## Where we are (Phase 8.4 complete)

The validated single-target chain
(DACE fetch -> filter -> nightly bin -> Keplerian fit -> activity
decorrelation -> quasi-periodic GP -> optional MCMC error bars)
is wrapped in survey infrastructure:

- `src/survey.py` -- survey driver. Scorecard per target: data
  quality, recovered planets vs literature (with error bars),
  vetoed candidate peaks, analytic detection limits. One target
  failing never kills the run. Outputs to `output/survey/`.
- `data/known_planets.json` + `src/reference_data.py` -- literature
  reference tables drive fitting; adding a curated target = editing
  data, not code.
- `data/shell_catalog.json` + `src/targets.py` -- generated catalog
  of 176 FGK dwarfs within 15 pc (SIMBAD TAP query), distance shells
  0-6 / 6-10 / 10-15 pc. Curated 10-star catalog merges on top.
- `tests/test_phase8.py` (43 tests) + Phase 7 suites = 136 offline
  regression tests, all green as of 5ebd852.

## Validation results (why we trust the chain)

| Target | Regime | Result |
|---|---|---|
| HD 20794 | K ~ 0.5 m/s BELOW activity | 3/3 planets, K within 0.6 sigma of Nari+2025 (MCMC) |
| 61 Vir | K 1.5-3.6 m/s ABOVE activity | 3/3 planets, K within 3% of Laliotis+2023 |
| Tau Ceti | zero confirmed planets | honest null: no forced fit, no false candidates, K_lim reported |

First 3-target survey run: 3/3 ok, ~6 min without MCMC.

## How to run things

```
# All offline tests (~10 s)
./venv/Scripts/python.exe -m pytest tests/test_phase7a.py tests/test_phase7b.py tests/test_phase7b2.py tests/test_phase7c.py tests/test_phase8.py -q

# Survey over named targets (GP on, no MCMC: ~2 min/target)
./venv/Scripts/python.exe -m src.survey --targets "HD 10700" "HD 20794" "HD 115617"

# Survey over a distance shell (needs data/shell_catalog.json)
./venv/Scripts/python.exe -m src.survey --shell 0 6

# MCMC error bars (~35 min/target; cap BLAS threads first)
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  ./venv/Scripts/python.exe -m src.survey --targets "HD 20794" --mcmc

# Regenerate the shell catalog (SIMBAD TAP, ~30 s)
./venv/Scripts/python.exe -c "from src.targets import build_shell_catalog; build_shell_catalog()"
```

## Next steps (in order)

1. **8.5 vetting depth** -- what stands between "peak list" and
   "credible candidate flags":
   - Per-peak false-alarm probability (current FAP is best-peak only)
   - GP-only null model (0 planets) to lower the noise floor on
     zero-planet targets like Tau Ceti (currently limited by
     unmodeled activity: RMS 3.96 m/s -> K_lim 0.76 m/s at 100 d,
     above the Feng+2017 candidate amplitudes)
   - Alias-of-candidate machinery (61 Vir 4.7/5.4/6.5 d cluster is
     plausibly inter-alias)
   - Multi-instrument confirmation requirement
2. **First full 0-6 pc shell run** (18 targets) -- after 8.5 so the
   candidate lists are credible. Many shell targets are uncurated;
   scorecards will mostly report data availability + limits.
3. **Optional**: juliet nested sampling for model comparison (8.1b);
   per-target process parallelism; raise mcmc_nsteps for
   publication-grade convergence.
4. **Gaia DR4 (late 2026)** -- epoch astrometry; revisit
   proper-motion-anomaly module when released.

## Critical gotchas (cost real debugging time)

- RadVel >= 1.4 caches vary flags at Posterior construction; setting
  `params[...].vary` later is silently ignored. Always go through
  `_sync_vary_flags()` in src/rv_keplerian.py.
- radvel.mcmc: pass `ensembles=3` explicitly (default 8 run
  SEQUENTIALLY with serial=True) and restore chain medians after
  (mcmc leaves the posterior at an arbitrary walker position).
  Both handled inside fit_keplerian.
- Cap BLAS threads (OMP/MKL/OPENBLAS_NUM_THREADS=4) for GP/MCMC:
  ~9x less CPU at the same speed.
- DACE needs dace-query >= 3.0 (2.x endpoints are dead); absolute
  RVs must be pre-centered; catastrophic outliers are normal
  (handled by rv_filter_instruments).
- NASA archive carries stale/disputed planet rows; curated
  data/known_planets.json is the authority for fitted planets.
- tc must be in the data's time system (RJD for DACE).

## Document map

- This file -- session entry point
- docs/architecture_and_data_sources.md -- master reference: every
  src module, what data is online vs static vs generated vs curated,
  and what each test file exercises (offline/mocked/network)
- docs/phase8_survey_plan.md -- Phase 8 design + per-step results
- docs/phase7c_plan.md -- chain upgrades (binning, decorrelation, GP)
- docs/thesis_final_summary.md -- project thesis at the March 2026
  pause (pre-survey; still correct on methodology)
- MEMORY.md (assistant memory, outside repo) -- condensed findings
  per phase with commit hashes
