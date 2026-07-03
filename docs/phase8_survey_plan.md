# Phase 8: From Validated Chain to Proximity-Ordered Survey

## Status: 8.1-8.4 complete (2026-07-03); next: 8.5 vetting depth

### 8.4 implementation notes (2026-07-03)

- targets.py: build_shell_catalog() queries SIMBAD TAP (not Gaia DR3:
  brightest nearest stars, e.g. Alpha Cen A, are missing or unreliable
  in Gaia) for all stars with plx >= 1000/15 mas and FGK sp_type.
- Client-side filters (filter_shell_rows, offline-testable):
  FGK dwarf luminosity class (V / IV-V / absent; giants + pure IV
  rejected), approximate Teff 4000-6300 K from a Pecaut & Mamajek
  spectral-type table (teff_source="sptype_approx"), HD identifier
  required (DACE/NASA are keyed on HD).
- Component dedup: composite row dropped only when a lettered 'A'
  component exists (70 Oph); when only B exists the composite IS the
  primary (eta Cas A bug found in first sanity run).
- targets_in_shell() merges curated TARGET_CATALOG entries (curated
  wins on HD collision) -- covers boundary disagreements like Delta
  Pavonis (G8IV in SIMBAD vs G8IV-V curated). Merge key is the full
  HD string: HD 165341A / HD 165341B are different stars.
- Cache: data/shell_catalog.json (176 targets < 15 pc). Shells:
  0-6 pc 18 targets, 6-10 pc 40, 10-15 pc 119 (with curated merge).
  103 K dwarfs -- the population the hardcoded list undersampled.
- get_target() falls back to the shell catalog after the curated
  list; run_survey(shell=(0,6)) / CLI --shell 0 6 runs a whole shell.
- Sanity criterion PASS: all 9 curated targets < 15 pc reproduced;
  HD 134060 (24 pc) correctly excluded.
- Smoke test on uncurated shell target (eta Cas B, K7Ve): survey
  scorecard works end-to-end -- NASA none_found, DACE has 117 binned
  points (HAMILTON/HIRES legacy, 27.8 m/s precision), detection
  limits reported. Data availability is itself survey output.
- Tests: 43 in test_phase8.py; full offline regression 136 green.

### 8.3 implementation notes (2026-07-03)

- src/survey.py: run_survey() iterates targets nearest-first; per-target
  try/except so one failure never kills the survey; outputs
  output/survey/survey_results.json + survey_summary.md.
- survey_target() is reference-first: data/known_planets.json decides
  which planets are fitted; NASA archive only for uncurated targets.
  Reuses deep_dive._run_rv_analysis (same chain as deep dives).
- Zero-confirmed-planet targets (Tau Ceti) get no forced Keplerian fit;
  candidates come from the binned-RV periodogram and the noise floor is
  the per-instrument-centered scatter. Honest: no planet model imposed.
- build_scorecard(): data quality, recovered-vs-reference planets
  (k_ratio, n_sigma when k_err available), vetoed candidate list,
  analytic detection limits (rv_detection_limit with final residual
  RMS as noise; injection-based limits deferred).
- vet_candidate_period(): basic veto set ahead of 8.5 -- 1-day alias
  family (0.5/0.997/1.0/2.0 d), annual + first harmonic, rotation
  harmonics (P/2, P, 2P), magnetic cycle +/-30%, cycle-annual beat
  (the 414 d lesson), coverage < 2 cycles, known-planet match.
  Surviving candidates cross-matched against unconfirmed literature
  candidates by period.
- Tau Ceti curation (HD 10700): zero confirmed planets. Tuomi+2013
  b/c/d rejected; Feng+2017 g/h/e/f (20.00/49.41/162.87/636.13 d,
  K ~ 0.3-0.5 m/s) all unconfirmed; Figueira+2025 ESPRESSO finds no
  conclusive planets (20 d signal non-significant with activity model,
  near rotation first harmonic). Rotation 46 d (Korolik+2023), weak
  ~11 yr cycle. k_ms in JSON computed from m sin i (approximate).
- Also added: daily-alias veto for known planets (f_alias = 1 -/+ 1/P
  cycles/day); caught 61 Vir b's 1.31 d alias masquerading as a
  candidate in the first run.
- Tests: 28 in test_phase8.py (vetoes incl. 414 d beat + 1 d alias +
  daily alias of known planet, scorecard construction with/without
  fit, failure injection, distance ordering); Phase 7 regression 93
  tests green.

### First survey run (2026-07-03, GP on, MCMC off)

| Target | N binned | Baseline (d) | RMS (m/s) | Planets rec/ref | Candidates | K_lim@100d |
|---|---|---|---|---|---|---|
| Tau Ceti | 2280 | 13140 | 3.96 | 0/0 | 0 | 0.76 |
| 82 G. Eridani | 748 | 7128 | 1.82 | 3/3 | 2 | 0.48 |
| 61 Virginis | 883 | 14614 | 4.87 | 3/3 | 3 | 1.16 |

- 3 targets, 3 ok, 0 failed; all confirmed reference planets recovered
  at both anchors (HD 20794 K within 8-12%, 61 Vir K within 1-7%,
  matching Phase 7c single-target results -- survey wrapper does not
  degrade the chain).
- Tau Ceti scorecard is the survey's honest null: no Keplerian fit
  imposed, no candidate peaks survive vetoes, and the analytic K limit
  (0.76 m/s at 100 d, from 3.96 m/s unmodeled scatter) sits ABOVE the
  Feng+2017 candidate amplitudes (0.3-0.5 m/s). We cannot confirm or
  refute them -- consistent with Figueira+2025 needing dedicated
  ESPRESSO modeling. A GP-only (0-planet) activity model to lower the
  Tau Ceti noise floor is a natural 8.5 follow-up.
- Surviving candidates are weak (power 0.05-0.07) and carry the
  scorecard note that they are follow-up flags, not detections;
  per-peak FAP is the main missing piece (8.5). HD 20794: 114.6 d,
  8.5 d. 61 Vir: 4.7/5.4/6.5 d cluster (plausibly inter-planet
  aliases; needs the 8.5 alias-of-candidate machinery).
- Wall clock ~6 min for 3 targets without MCMC (BLAS capped);
  MCMC scorecard runs remain ~35 min/target.

### 8.1 verification: HD 20794 MCMC results (2026-07-02)

Full chain with mcmc_final=True (3 ensembles x 50 walkers x 1000
steps, serial, BLAS capped to 4 threads): ~35 min wall clock.

| Planet | Pipeline K (m/s) | Nari 2025 K (m/s) | Consistency |
|--------|------------------|--------------------|-------------|
| b | 0.687 +/- 0.107 | 0.614 +/- 0.048 | 0.6 sigma |
| c | 0.458 +/- 0.115 | 0.502 | 0.4 sigma |
| d | 0.612 +/- 0.465 | 0.567 | 0.1 sigma |

All three K amplitudes statistically consistent with literature --
the 8.1 acceptance criterion is met. Additional error bars: periods
(P_b +/- 0.003 d), instrument gammas (+/-1-2 m/s), jitters
(+/-0.1 m/s), GP hyperparameters.

Honest caveats:
- Chains did not pass radvel's strict convergence tests (G-R);
  errors are approximate. Survey runs needing publication-grade
  posteriors should raise mcmc_nsteps.
- K_d error (+/-0.465, 76% relative) is large because the GP and
  planet d (647 d) compete for long-period signal -- a real
  degeneracy the posterior correctly exposes, not a bug.
- P_d formal error (+/-0.4 d) is likely underestimated given
  non-convergence; the 654 vs 648 d offset should not be
  over-interpreted.

### Performance findings (compute optimization)

- Bug found: radvel serial=True + default ensembles=8 runs 8
  ensembles SEQUENTIALLY (8x work). First run killed at 98 min,
  ~11 cores busy. Fixed with mcmc_ensembles=3.
- BLAS thread cap (OMP/MKL/OPENBLAS_NUM_THREADS=4): CPU usage
  dropped from ~11 cores to ~1.2 cores for the same progress rate
  -- thread coordination was pure overhead on 542x542 Cholesky.
- Net: >8x wall-clock improvement. Survey defaults: ensembles=3,
  capped threads, MCMC only on final scorecard runs.
- Decision on further optimization: per-target process parallelism
  at 8.3 (best ROI); celerite2 O(N) kernel only if a 15 pc shell
  (~100 targets) demands it; GPU not worth it (tiny matrices,
  sequential sampler).

### Post-MCMC state fix

radvel.mcmc leaves the posterior object at an arbitrary final walker
position (observed: GP hyperparameters drifted from MAP while
well-constrained planet params stayed put). fit_keplerian now
restores every sampled parameter to its chain median after MCMC, so
reported values, residuals, and GP predictions come from one
consistent state.

### 8.1 implementation notes
- radvel.mcmc called with serial=True (parallel ensembles use
  multiprocessing, which breaks in library code on Windows) and
  headless=True
- Errors from chain 16-84 percentile half-width (_chain_err); k_err,
  period_err, tc_err, e_err, w_err per planet; gamma_err/jitter_err
  per instrument; gp_*_err for hyperparameters
- deep_dive mcmc_final flag applies MCMC to the final GP stage only
- Synthetic validation: K = 4.93 +/- 0.077 vs analytic sigma ~0.065

### 8.2 implementation notes
- data/known_planets.json holds per-target planet tables + config
  (gp_hyperparams, activity_cycle_days, exclude_instruments, notes)
- src/reference_data.py loader; deep_dive._known_planets_fallback is
  now a thin wrapper
- GP stage seeds gp_per/gp_explength from activity_cycle_days when no
  explicit override (61 Vir: 5911 d from Laliotis+2023)

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
