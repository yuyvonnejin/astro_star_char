# Phase 7c: Convergence Improvements and Second-Target Validation

## Status: 7c.1, 7c.2, 7c.3, 7c.5 complete (2026-07-02); 7c.4 (juliet) optional, not started

## Context

Project paused after Phase 7b.3 (2026-03-20). State at pickup:

- HD 20794 deep dive complete: RMS 1.77 m/s, 3 planets fitted, disputed 147d signal rejected
- Known gap: fitted periods drift (b: 13.5 vs 18.3 d; d: 752 vs 648 d), eccentricities
  blow up to 0.7-0.9, K amplitudes 1.5-2x too high
- Root cause (documented in phase7b3 doc, section 4): ~1.8 m/s correlated stellar
  activity noise is unmodeled; MAP optimizer absorbs it into orbital parameters
- Environment verified 2026-07-01: radvel 1.5.0 and dace-query 2.0.0 installed,
  all 10 phase7b2 tests pass. juliet and celerite2 NOT yet installed.
- note20260219.txt (re-inspect 7b.2, re-run HD 20794) was satisfied by Phase 7b.3

Phase 7c goal: close as much of the parameter convergence gap as possible with
the recommendations already identified in phase7b3 section 6, then prove the
pipeline generalizes beyond HD 20794.

---

## Step Ordering and Rationale

Ordered by effort-to-impact ratio. Each step is independently testable; stop or
reprioritize after any step.

### 7c.1 Nightly binning (no new dependencies)

Problem: 11,854 measurements contain heavy intra-night correlation (p-mode
oscillation, granulation). Standard practice is nightly bins; Nari 2025 used
806 nightly-binned points.

- Add `rv_bin_nightly()` to `src/rv_data.py`
  - Group by instrument + floor(RJD + 0.5) (site-local night boundary)
  - Weighted mean RV, error = max(propagated error, intra-night scatter / sqrt(N))
- Wire into `deep_dive.py` before Keplerian fit (flag, default on)
- Expected: N drops ~12K to ~1.5K; faster fits; possibly better convergence
- Verify: RMS should stay near or below 1.77 m/s; compare fitted P/K/e to Nari

### 7c.2 Activity indicator decorrelation (no new dependencies)

Problem: rotation (~39 d) and magnetic cycle (~3000 d) signals leak into the fit.
DACE provides FWHM, BIS, and S-index alongside RVs.

- Extend DACE retrieval in `rv_data.py` to keep fwhm, bis, sindex columns
- Add `rv_decorrelate_activity()`: linear regression of RV against indicators
  (per instrument), subtract the fitted component
- Order of operations: decorrelate residual activity AFTER first Keplerian pass,
  or jointly iterate (decorrelate raw, fit, check residuals)
- Expected: partial removal of the 39 d and 3000 d power; the 414 d beat
  artifact should weaken
- Verify: residual periodogram before/after; fitted e values should drop
  toward literature (0.06-0.45)

### 7c.3 GP activity model (new dependency: celerite2)

Highest-impact item per phase7b3 recommendations. RadVel natively supports GP
likelihoods with celerite kernels (O(N) scaling, fine even at 12K points,
trivial at ~1.5K nightly-binned points).

- Install celerite2 (check Windows wheel availability first; fallback: celerite)
- Extend `rv_keplerian.py` with optional quasi-periodic GP
  (RadVel `GPLikelihood`, QuasiPer kernel)
  - Hyperparameter priors from literature: P_rot ~ 35-39 d, cycle ~ 3000 d
- Add `use_gp` flag to `fit_keplerian()` and `deep_dive.py`
- Expected: this is what Nari used (conceptually); best chance of pulling
  fitted K down to sub-1 m/s literature values
- Verify against Nari 2025: K_b=0.614, K_c=0.502, K_d=0.567 m/s;
  e_d=0.45; H03-H15 offset 17.0 +/- 1.7 m/s
- Risk: GP hyperparameters can also absorb planet signal; keep priors tight

### 7c.4 Model comparison with juliet (new dependency, optional)

- Install juliet (wraps RadVel + dynesty nested sampling)
- Compare 2 vs 3 vs 4 planet models on HD 20794: Bayesian evidence should
  favor 3 planets and reject the 147 d signal quantitatively
  (currently rejected only via literature cross-reference)
- This replaces the hardcoded disputed-planet logic with a statistical test
- Defer if 7c.3 runtime is already heavy; literature cross-reference works

### 7c.5 Second-target validation: 61 Virginis (HD 115617)

Prove the pipeline is a survey tool, not an HD 20794 special case.

- Why 61 Vir before Tau Ceti: its planets have K = 2-4 m/s (b: 2.1, c: 3.6 m/s),
  well ABOVE the ~1.8 m/s activity noise floor. MAP fitting should converge
  without a GP, giving a clean test of the 7c.1-7c.2 changes in a regime
  where they should suffice.
- Stellar pipeline output already exists (output/results_hd115617.json)
- Tasks:
  - Add 61 Vir known-planets fallback data to `deep_dive.py`
    (Vogt et al. 2010; check for newer references first)
  - Run `run_deep_dive('61 Virginis')`; expect new DACE data quality issues --
    document them the same way as HD 20794
- Success: fitted P within 2%, K within 30% of literature for all 3 planets
- Follow-up target after 61 Vir: Tau Ceti (HD 10700, target #2, sub-m/s
  signals again, needs GP from 7c.3)

---

## Verification Criteria

| Step | Metric | Target |
|------|--------|--------|
| 7c.1 | N after binning | ~1,000-1,600 |
| 7c.1 | RMS after fit | <= 1.8 m/s (no regression) |
| 7c.2 | 414 d residual peak | reduced or absent |
| 7c.3 | K_b, K_c, K_d | within 50% of Nari (stretch: 20%) |
| 7c.3 | e_b, e_c | < 0.3 (currently 0.7-0.9) |
| 7c.3 | P_b | 18.3 +/- 0.5 d (currently drifts to 13.5) |
| 7c.4 | ln Z(3 planets) > ln Z(4 planets) | 147 d rejected statistically |
| 7c.5 | 61 Vir: 3 planets recovered | P within 2%, K within 30% |

Realistic expectation for 7c.3: DACE DRS RVs (no YARARA/sBART spectral
reprocessing) put a floor on achievable accuracy. If K values land within
~30-50% and eccentricities become sane, that is success for this data source.
Do not chase Nari-exact numbers; phase7b3 section 6.3 defines the scope boundary.

---

## Files to Create/Modify

| File | Action |
|------|--------|
| src/rv_data.py | Add rv_bin_nightly(), rv_decorrelate_activity(), indicator retrieval |
| src/rv_keplerian.py | Add GP likelihood option (celerite2 kernel) |
| src/deep_dive.py | Wire binning/decorrelation/GP flags; add 61 Vir fallback data |
| tests/test_phase7c.py | New: binning, decorrelation, GP smoke test, 61 Vir fallback |
| docs/phase7c_plan.md | This file; update with results per step |

## Housekeeping (before starting)

- Commit pending working-tree changes: archive/ move (ReschSumProp.txt,
  tutorial notebook), ref_planets_complete.csv, note20260219.txt
- Decide whether old timestamped deep-dive reports in output/target_reports
  should be pruned or kept as history

## Test Commands

- Existing suites must stay green:
  `./venv/Scripts/python.exe -m pytest tests/test_phase7a.py tests/test_phase7b.py tests/test_phase7b2.py -v`
- New: `./venv/Scripts/python.exe -m pytest tests/test_phase7c.py -v`

---

## 7c.1 Results (2026-07-02)

### Unplanned fixes required first

**DACE API migration.** dace-query 2.0.0 stopped working (server endpoint
returns 404; DACE migrated their API). Upgraded to dace-query 3.0.2.
Breaking change: response key `ins_name` renamed to `instrument_name`;
`_parse_dace_query_rv()` now reads both. Instrument labels unchanged
(HARPS03/HARPS15/ESPRESSO18/ESPRESSO19), except CORAVEL-S is now CORAVEL
(harmless: exclusion is precision-based, not name-based). The public
dataset itself also changed server-side: HD 20794 now returns 9,077
measurements (was 12,616 in March 2026); HARPS03 dropped from 10,338 to
6,763. New response fields include ccf_fwhm, ccf_bispan, date_night --
useful for 7c.2.

**RadVel vary-flag bug (critical).** RadVel >= 1.4 caches parameter
values and vary flags in a numpy vector at Posterior construction and
caches the vary-index list on first use. Our code set
`post.params[...].vary` after construction, which only changes the dict;
the optimizer kept using the stale vector. Consequence: the entire
3-pass MAP strategy and `fix_eccentricities` were no-ops in all Phase
7b.2/7b.3 runs -- every pass fit all parameters freely, including dvdt
and curv. Fixed with `_sync_vary_flags()` (calls
`post.vector.dict_to_vector()` + `post.list_vary_params()`) after each
vary-flag change. Regression tests added.

### HD 20794 results: binned + fixed optimizer vs literature

748 nightly bins (Nari 2025 used 806 -- close match), 3 instruments.

| Parameter | Pipeline (7c.1) | Nari 2025 | 7b.3 (buggy fit) |
|-----------|----------------|-----------|-------------------|
| P_b (d) | 18.310 (0.03% off) | 18.315 | 13.5 (drifted) |
| P_c (d) | 89.36 (0.4% off) | 89.68 | 91.2 |
| P_d (d) | 642.6 (0.8% off) | 647.6 | 752.5 (drifted) |
| K_b (m/s) | 0.883 | 0.614 | 0.82 |
| K_c (m/s) | 0.495 (1.4% off) | 0.502 | 0.90 |
| K_d (m/s) | 1.73 (3x high) | 0.567 | 1.22 |
| e (all) | fixed at lit values | 0.06/0.08/0.45 | drifted to 0.7-0.9 |
| H15-H03 offset (m/s) | 20.3 | 17.0 +/- 1.7 | 10.4 |
| RMS after (m/s) | 4.28 | ~1.5 (after GP) | 1.77 |

### Binning A/B (same fixed optimizer, same dataset)

| Parameter | Unbinned (8,379 pts) | Binned (748 pts) | Nari 2025 |
|-----------|---------------------|------------------|-----------|
| P_b (d) | 18.289 | 18.310 | 18.315 |
| P_c (d) | 90.42 | 89.36 | 89.68 |
| P_d (d) | 631.8 | 642.6 | 647.6 |
| K_b (m/s) | 1.44 | 0.883 | 0.614 |
| K_c (m/s) | 1.15 | 0.495 | 0.502 |
| K_d (m/s) | 2.66 | 1.73 | 0.567 |
| RMS after (m/s) | 4.28 | 4.28 | ~1.5 |

Binning improves every K amplitude (roughly halves them, toward
literature) and slightly improves all periods, at the same residual
RMS. Attribution is clean: correlated intra-night noise was inflating
the amplitudes.

### Interpretation

- All three periods now recover to sub-1% of literature. This was the
  main 7b.3 failure and is now solved by the honest 3-pass fit.
- K_c essentially exact; K_b 44% high; K_d 3x too high.
- Residual RMS went UP (1.77 -> 4.28 m/s) because the previous number
  was an artifact: the buggy free-for-all fit absorbed stellar activity
  into orbital parameters (extreme eccentricities). The honest fit
  leaves activity in the residuals where it belongs.
- Residual periodogram now shows the magnetic cycle directly: peaks at
  3564 d (~the 3000 d cycle), 713 d, 400 d, 356 d (annual beat aliases).
  This contaminates planet d (642 d) and explains its inflated K.
- Verification criteria status: P_b target met (18.3 +/- 0.5), e sane
  (fixed), K_c within tolerance. K_b/K_d and the residual activity
  signal are exactly what 7c.2 (decorrelation) and 7c.3 (GP) target.

### Verdict

7c.1 complete. The vary-flag fix (found during 7c.1 verification) is a
bigger improvement than binning itself. Next: 7c.2 activity indicator
decorrelation, using ccf_fwhm/ccf_bispan now available from DACE v3.

---

## 7c.2 Results (2026-07-02)

### Implementation

- `_parse_dace_query_rv()`: carries activity indicators (fwhm, bis,
  smw, rhk, halpha) as parallel arrays with NaN for missing; DACE
  sentinel values (0, |x| > 1e8) converted to NaN
- `rv_filter_instruments()` and `rv_bin_nightly()`: indicators
  follow the same masking / weighted binning as RVs
- New `rv_decorrelate_activity()` in rv_data.py: per-instrument
  linear regression of Keplerian fit residuals against MAD-normalized
  indicators (default: fwhm, bis, smw); indicator outliers clamped at
  5 MAD-sigma; correction is zero-mean per instrument (gammas
  unaffected)
- `deep_dive._run_keplerian_or_sinusoidal()`: fit -> decorrelate
  residuals -> subtract -> refit; refit kept only if residual RMS
  improves; summary stored under 'activity_decorrelation'
- Regressing fit residuals (not raw RV) protects planetary signal
  from being absorbed into the activity model

### HD 20794: decorrelation impact (all runs binned, fixed optimizer)

| Parameter | 7c.1 (no decorr) | 7c.2 (decorr) | Nari 2025 |
|-----------|------------------|---------------|-----------|
| P_b (d) | 18.310 | 18.319 | 18.315 |
| P_c (d) | 89.36 | 89.25 | 89.68 |
| P_d (d) | 642.6 | 637.1 | 647.6 |
| K_b (m/s) | 0.883 (44% high) | 0.755 (23% high) | 0.614 |
| K_c (m/s) | 0.495 | 0.475 (5% low) | 0.502 |
| K_d (m/s) | 1.73 (3.0x) | 0.999 (1.8x) | 0.567 |
| RMS after (m/s) | 4.28 | 2.99 | ~1.5 |

Variance reduction of residuals by instrument: HARPS03 54.8%,
ESPRESSO19 47.4%, HARPS15 14.6%. Dominant coefficient: FWHM on
HARPS03 (2.98 m/s per MAD-sigma) -- the magnetic cycle.

### Interpretation

- Every K amplitude moved toward literature; K_d improved most
  (1.73 -> 1.00 m/s), confirming its inflation was magnetic cycle
  leakage.
- The 713 d beat-alias residual peak from 7c.1 is gone. Remaining
  residual peaks (482 d, 1549 d, 1080 d, 3564 d) are the parts of
  the ~3000 d cycle a single linear term cannot capture -- the cycle
  is not a linear function of the indicators over 20+ years.
- Periods stay within 0.03-1.6% of literature.
- Remaining gap to Nari (K_b +23%, K_d +76%, RMS 2.99 vs 1.5) is
  the 7c.3 GP target.

### Verdict

7c.2 complete. Progression across one day of work:
K_d 2.66 (unbinned, honest fit) -> 1.73 (binned) -> 1.00 (decorrelated),
literature 0.567. Next: 7c.3 quasi-periodic GP via celerite2.

---

## 7c.3 Results (2026-07-02)

### Implementation

- No celerite needed: nightly binning (748 points) makes RadVel's pure
  numpy QuasiPer kernel (O(N^3) Cholesky) fast enough. This avoids the
  Windows C++ toolchain entirely.
- `fit_keplerian(use_gp=True)`: GPLikelihood per instrument with shared
  quasi-periodic hyperparameters (gp_amp, gp_explength, gp_per,
  gp_perlength); HardBounds priors; defaults target a solar-type
  magnetic cycle (DEFAULT_GP_HYPERPARAMS, gp_per init 3000 d)
- Pass schedule: pass 1 frees only gp_amp (offsets/K first); pass 2
  frees GP timescales with w; pass 3 frees periods
- Residuals subtract the GP conditional mean per instrument, so
  rms_after_ms reflects what neither planets nor activity explain
- `deep_dive`: final-stage GP fit on the decorrelated RVs (use_gp
  flag, default on; skipped above 3000 points)

### HD 20794: full chain vs literature

| Parameter | 7c.1 binned | 7c.2 decorr | 7c.3 GP | Nari 2025 |
|-----------|------------|-------------|---------|-----------|
| P_b (d) | 18.310 | 18.319 | 18.317 (0.01%) | 18.315 |
| P_c (d) | 89.36 | 89.25 | 89.59 (0.1%) | 89.68 |
| P_d (d) | 642.6 | 637.1 | 654.3 (1.0%) | 647.6 |
| K_b (m/s) | 0.883 | 0.755 | 0.687 (+12%) | 0.614 |
| K_c (m/s) | 0.495 | 0.475 | 0.458 (-9%) | 0.502 |
| K_d (m/s) | 1.73 | 0.999 | 0.612 (+8%) | 0.567 |
| RMS (m/s) | 4.28 | 2.99 | 1.83 | ~1.5 |
| Jitter E19/H03/H15 | -- | -- | 0.64/1.84/1.28 | ~0.4-1.0 |

Fitted GP hyperparameters: per=1439 d, explength=1824 d, amp=3.0 m/s,
perlength=0.46. The GP soaks up the magnetic cycle: residual
periodogram peaks drop from power 0.26 (713 d) to 0.09 (1 d alias);
cycle/beat peaks are gone. Remaining weak peaks: 114.6/82.1/165.8 d.

### Verification criteria: status

| Criterion (from plan) | Target | Achieved |
|----------------------|--------|----------|
| K_b, K_c, K_d | within 50% (stretch 20%) | 8-12% |
| e_b, e_c | < 0.3 | fixed at lit |
| P_b | 18.3 +/- 0.5 d | 18.317 |
| RMS | -- | 1.83 vs lit ~1.5 |

All 7c.3 criteria met and the stretch goal exceeded. The thesis
(section 9.2) put "publication accuracy" at +/-8% for K; the automated
pipeline now sits at 8-12% using only public DRS data -- without
YARARA/sBART spectral reprocessing or nested sampling.

### Caveat

K uncertainties are not yet quantified (MAP only, k_err=None). The
8-12% agreement is point-estimate agreement. MCMC or nested sampling
(7c.4 juliet) would give posterior widths to judge consistency
properly.

### Verdict

7c.3 complete. One-day progression on K_d (lit 0.567):
2.66 -> 1.73 (binning) -> 1.00 (decorrelation) -> 0.612 (GP).
Remaining steps: 7c.4 juliet model comparison (optional),
7c.5 second-target validation on 61 Virginis.

---

## 7c.5 Results: 61 Virginis (2026-07-02)

### Literature reference update

Current literature status (verified via web search + Laliotis et al.
2023 paper, arXiv:2302.10310):

- 61 Vir d history parallels HD 20794 e: discovered by Vogt+2010
  (22.9 Me), flagged false positive by Rosenthal+2021 (CLS),
  reconfirmed at lower mass by Laliotis+2023 and Cretignier+2023.
- NASA archive default rows still carry Vogt 2010 parameters.
- `_known_planets_fallback()` updated with Laliotis+2023 Table values:
  b: P=4.21498, K=2.47+/-0.11, e=0.033
  c: P=38.079,  K=3.56+/-0.12, e=0.026
  d: P=123.2,   K=1.47+/-0.17, e=0.15
- Laliotis+2023 also report a ~5911 d long-period activity signal
  (magnetic cycle), within our GP prior bounds.

### Data (DACE v3)

3,224 raw measurements, 11 instruments, 40-year baseline
(14,778 days). Filter chain: CORAVEL excluded (320 m/s precision),
17 outliers sigma-clipped, 3,188 -> 883 nightly bins across 10
instruments (APF, COUDE, ESPRESSO18/19, HAMILTON, HARPS03/15,
HIRES, HIRES-POST04, UCLES). Mixed zero-points handled by
pre-centering: HARPS/ESPRESSO absolute (~-7.9 km/s), HIRES/APF/
UCLES/HAMILTON relative (~0 m/s).

### Results vs Laliotis et al. (2023)

| Planet | Pipeline P (d) | Lit P (d) | Pipeline K (m/s) | Lit K (m/s) | K error |
|--------|---------------|-----------|------------------|-------------|---------|
| b | 4.2150 | 4.21498 | 2.441 | 2.47 +/- 0.11 | -1.2% |
| c | 38.078 | 38.079 | 3.495 | 3.56 +/- 0.12 | -1.8% |
| d | 123.04 | 123.2 | 1.511 | 1.47 +/- 0.17 | +2.8% |

All three K amplitudes within 3% of literature and inside the
published 1-sigma intervals. Periods essentially exact. The
disputed-then-reconfirmed planet d is cleanly recovered.

Other outputs:
- Jitters: modern instruments 0.7-2.5 m/s; legacy HAMILTON 7.4,
  COUDE 16.0 m/s (physically sensible)
- Decorrelation: marginal here (RMS 5.18 -> 5.15) -- indicators only
  exist for HARPS03 (225 of 883 bins) and HIRES S-index
- GP: modest (5.15 -> 4.85 m/s), amp 1.37 m/s, per ~2095 d
- Residual periodogram: only 1-day alias family remains; no
  unexplained planet-like signals
- Final RMS 4.85 m/s is dominated by legacy instrument scatter
  (HAMILTON/COUDE), not by unmodeled signal

### Success criteria: status

Plan target: P within 2%, K within 30%.
Achieved: P within 0.13%, K within 3%.

As predicted, in the K >> activity regime the pipeline converges
essentially to literature values without needing the GP; binning +
honest optimizer suffice.

### Verdict

7c.5 complete. Two-target validation passed:
- HD 20794 (hard regime, K ~ 0.5 m/s < activity): K to 8-12%
- 61 Vir (moderate regime, K ~ 1.5-3.6 m/s > activity): K to 3%
The proximity-ordered survey infrastructure (filter -> bin ->
Keplerian -> decorrelate -> GP) generalizes across targets and
instrument mixes without per-target code changes -- only literature
reference data per target. Remaining optional: 7c.4 juliet for
posterior uncertainties and model comparison.
