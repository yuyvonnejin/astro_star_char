# Phase 8: From Validated Chain to Proximity-Ordered Survey

## Status: Draft (2026-07-02)

## Vision

Project goal (owner, 2026-07-02): meaningfully participate in, and if
possible contribute to, exoplanet search. Strategy: since detection
bias means only the most obvious planets are findable in survey data,
prioritize impact -- search nearest stars first, expanding gradually
outward in distance shells.

Phase 7c validated the single-target chain on two regimes:
- HD 20794 (K ~ 0.5 m/s below activity): K to 8-12% of literature
- 61 Vir (K ~ 1.5-3.6 m/s above activity): K to 3% of literature

Phase 8 turns that chain into survey infrastructure. Realistic
contribution ladder this enables:
1. Data quality findings (already have: DACE corruption catalog)
2. Uniform survey with per-target detection limits (value even with
   zero discoveries)
3. Candidate flagging for community follow-up (not discovery claims)
4. Positioned for new data (Gaia DR4 late 2026, ongoing DACE/ESO
   releases) -- the highest-probability path to a real find

## Step Ordering

Uncertainties come first: every survey product (significant residual,
detection limit) is a statistical claim; MAP point estimates cannot
support them.

### 8.1 Parameter uncertainties (MCMC)

- Exercise the existing `run_mcmc=True` path in fit_keplerian()
  (written in 7b.2, never validated; vary-flag fix now makes it
  meaningful). radvel.mcmc supports GPLikelihood.
- Wire MCMC into the deep-dive final fit stage (flag, default off for
  speed; on for survey scorecards).
- Report k_err, period_err per planet from chain percentiles;
  convergence via radvel's Gelman-Rubin statistics.
- Check jitter/K priors behave in chains (jitter can wander negative
  at MAP boundary; add priors if chains misbehave).
- Verify: 61 Vir K_err magnitudes comparable to Laliotis+2023
  published errors (+/-0.11 to 0.17 m/s); HD 20794 K values
  consistent with Nari within joint error bars.
- Optional 8.1b: juliet nested sampling for Bayesian evidence
  (model comparison: N vs N+1 planets). Defer unless needed.

### 8.2 Reference data as data, not code

- Move `_known_planets_fallback()` tables from Python dicts to
  `data/known_planets.json`, keyed by HD identifier; keep schema
  (period, k_ms, e, confirmed, reference, notes).
- Add per-target config in the same file: GP hyperparameter
  overrides (activity cycle period), instrument exclusions,
  binning options.
- Loader with schema validation; deep_dive.py reads via loader.
- Adding a survey target then requires editing data, not code.

### 8.3 Survey driver

- New `src/survey.py`: `run_survey(targets, rv_only=True, ...)`
  - Iterate catalog in distance order
  - Per target: RV chain (filter -> bin -> Keplerian -> decorrelate
    -> GP -> optional MCMC); optionally TESS/PMa later
  - Per-target try/except: one failure never kills the survey
- Scorecard per target:
  - Data: n_measurements, instruments, baseline, quality flags
  - Planets: recovered P/K/e vs reference (when known), with errors
  - Candidates: residual peaks above FAP threshold, minus known
    aliases (1 d family, annual) and activity periods
  - Detection limit: K sensitivity vs period from residuals
- Aggregate outputs: `output/survey/survey_results.json` +
  markdown summary table across targets.
- First run: HD 20794, 61 Vir (regression anchors) + Tau Ceti
  (HD 10700, new hard case: 4 disputed sub-m/s planets).

### 8.4 Distance-shell catalog expansion

- `targets.py`: generate catalog by query instead of hardcoded list.
  - Gaia/SIMBAD TAP: main-sequence FGK (Teff 4000-6300 K, include
    K dwarfs -- better RV targets than G), distance shells
    (< 6 pc, 6-10 pc, 10-15 pc, ...), ordered by distance
  - Cross-match HD/HIP/TIC identifiers for archive queries
- Survey runs shell by shell; a target without usable RV data gets
  a scorecard saying so (data availability is itself survey output).

### 8.5 Candidate vetting rules (as they emerge)

- Codify the vetoes learned so far: 1-day alias family, annual
  aliases and beats (414 d lesson), stellar rotation period,
  magnetic cycle and harmonics, instrument-change epochs.
- A "candidate" must survive all vetoes AND have FAP below threshold
  AND appear in more than one instrument when coverage allows.

## Verification Criteria

| Step | Check |
|------|-------|
| 8.1 | 61 Vir K_err within ~2x of published errors; chains converge (G-R < 1.03) |
| 8.1 | Recovered K consistent with literature within joint 2-sigma |
| 8.2 | All Phase 7 tests pass with JSON-backed reference data |
| 8.3 | Survey run over 3 targets produces 3 scorecards, no crash on failure injection |
| 8.4 | Shell query reproduces the 10 hardcoded targets (sanity) plus K dwarfs |
| 8.5 | HD 20794 414 d beat and 61 Vir 1 d alias both auto-vetoed |

## Files

| File | Action |
|------|--------|
| src/rv_keplerian.py | Validate/repair MCMC path; chain-based errors |
| src/deep_dive.py | MCMC flag; read reference data via loader |
| data/known_planets.json | New: reference tables + per-target config |
| src/survey.py | New: survey driver + scorecard + aggregation |
| src/targets.py | Shell-based catalog generation |
| tests/test_phase8.py | New test suite per step |
| docs/phase8_survey_plan.md | This file; results appended per step |

## Test Commands

- `./venv/Scripts/python.exe -m pytest tests/test_phase8.py -v`
- Regression: all Phase 7 suites must stay green
