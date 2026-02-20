# Phase 7b.2: Literature Review, Methodology Gaps, and Refinement Plan

## Status: Implemented (Core)

Implementation completed: RadVel Keplerian fitting (rv_keplerian.py), CORAVEL exclusion
(rv_filter_instruments in rv_data.py), deep_dive.py integration, 10 unit tests passing.
Remaining: activity GP modeling (7b.3), model comparison via juliet (7b.3), alias analysis.

---

## 1. Literature Review: HD 20794 Planetary System

### 1.1 Discovery and Confirmation Timeline

| Paper | Year | Instruments | N_obs | Planets claimed | Key methodology |
|-------|------|-------------|-------|-----------------|-----------------|
| Pepe et al. (A&A 534, A58) | 2011 | HARPS | 173 | b(18.3d), c(40.1d), d(89.7d) | GLS periodogram, iterative Keplerian, e=0 fixed |
| Feng, Tuomi & Jones (A&A 605, A103) | 2017 | HARPS (TERRA) | 713 nightly | b(18d), d(89d), e(147d), f(330d tentative) | MA(4) noise model, AM-MCMC, BF>150 |
| Cretignier et al. (A&A 678, A2) | 2023 | HARPS (YARARA V2) | ~700 nightly | HZ candidate ~600d | Line-by-line PCA decorrelation |
| Nari et al. (A&A 693, A297) | 2025 | HARPS(YARARA)+ESPRESSO(sBART) | 806 nightly | b(18.3d), c(89.7d), d(647.6d) | GP activity model, FIP, nested sampling |
| Lovis et al. (A&A, 2025) | 2025 | ESPRESSO | 63 nights | -- | 10 cm/s short-term precision demonstrated |

### 1.2 Planet Parameters: Literature vs Our Pipeline

| Planet | Nari 2025 K (m/s) | Nari 2025 M (Me) | Our pipeline K (fitted) | Our pipeline M (Me) | Problem |
|--------|-------------------|------------------|------------------------|---------------------|---------|
| b (18.3 d) | 0.614 +/- 0.048 | 2.15 +/- 0.17 | 402.5 m/s | 2.15 (from NASA) | Fitted K is 656x too high |
| c/d (89.7 d) | 0.502 | 2.98 +/- 0.29 | 388.4 m/s | 2.98 (from NASA) | Fitted K is 774x too high |
| e (147.0 d) | NOT CONFIRMED | -- | 748.7 m/s | 4.77 (from NASA) | Planet likely not real |
| d/f (647.6 d) | 0.567, e=0.45 | 5.82 +/- 0.57 | 628.6 m/s | 5.82 (from NASA) | Fitted K is 1109x too high |

**Critical finding**: The sinusoidal fitter in `rv_subtract_sinusoids()` returns amplitudes
hundreds of m/s because it is fitting the residual instrument offset structure, not planet
signals. The true K amplitudes are all sub-1 m/s, buried under ~8500 m/s residual RMS.

### 1.3 Planet Status: What We Got Wrong

Our NASA Exoplanet Archive query returned 4 planets (b, d, e, f). Per Nari et al. (2025):

| Signal | NASA Archive | Nari 2025 Status | Notes |
|--------|-------------|------------------|-------|
| 18.3 d (b) | Confirmed | **Confirmed** | K=0.614 m/s, circular |
| 40.1 d (c) | Listed | **Rejected** | Stellar rotation (~39 d), false positive |
| 89.7 d (d) | Confirmed | **Confirmed** (renamed "c") | K=0.502 m/s, circular |
| 147.0 d (e) | Listed | **Not confirmed** | "worsened the evidence" when added to model |
| 330 d (f) | Not listed | **Unconfirmed** | Near annual alias, activity indicators present |
| 647.6 d | Listed as "f" | **Confirmed** (renamed "d") | K=0.567 m/s, e=0.45, HZ planet |

**Our pipeline used 4 periods for subtraction (18.3, 89.7, 147.0, 647.6 d), but the
147.0-day signal is likely not a real planet.** The correct set is 3 planets: 18.3, 89.7, 647.6 d.

### 1.4 The 414-Day Period: Solved

Our dominant periodogram peak at 413.98 days is almost certainly a **beat frequency between
the annual observing cadence and the stellar magnetic activity cycle**:

- Stellar magnetic cycle: ~3000 days (Nari 2025, detected in FWHM, S-index, BIS, CCF contrast)
- Annual observing cadence: 365.25 days
- Beat period: 1/(1/365.25 - 1/3000) = **416 days**

This matches our observed 413.98 days within uncertainty. It is NOT a planet, NOT a simple
alias of any known planet, and NOT a new signal.

Supporting evidence:
- Nari et al. (2025) found no additional Keplerian signals with P > 2000 days
- No linear RV trend detected (rules out massive outer companion causing the 414d signal)
- The 414-day peak disappeared in Nari's analysis which used GP modeling for the activity cycle

### 1.5 Stellar Activity Parameters (from Nari 2025)

| Parameter | Value | Source |
|-----------|-------|--------|
| log R'HK | -4.98 +/- 0.02 | Inactive regime |
| Rotation period | 35-39 days | FWHM: 35.0d, BIS: 38.8d |
| Magnetic cycle | ~3000 days (~8.2 yr) | FWHM, S-index, BIS, CCF contrast |
| RV jitter (ESPRESSO) | ~0.40 m/s | Lovis et al. 2025 |
| RV jitter (HARPS) | ~0.7-1.0 m/s | Pepe et al. 2011 |

### 1.6 Instrument Offset Handling in Nari 2025

Nari et al. split the data into **3 instrument epochs** (not our 5):

| Instrument | Our label | Nari label | N (Nari) | Notes |
|------------|-----------|-----------|----------|-------|
| CORAVEL-S | CORAVEL-S | Excluded | 0 | Too low precision (280 m/s median err) |
| HARPS pre-upgrade | HARPS03 | H03 | 512 nights | Pre-2015 fiber change |
| HARPS post-upgrade | HARPS15 | H15 | 231 nights | Post-2015 fiber change |
| ESPRESSO pre-intervention | ESPRESSO18 | E18 | ~34 nights | Pre-June 2019 |
| ESPRESSO post-intervention | ESPRESSO19 | E19 | ~661 obs | Post-June 2019 |

Key: Nari fitted H03-to-H15 offset as a free parameter (result: 17.0 +/- 1.7 m/s).
ESPRESSO epochs were also treated with separate offsets.

**Our pipeline included CORAVEL (16 measurements, 280 m/s median error) which is
~2000x worse precision than ESPRESSO. This dominates the RMS and corrupts the fit.**

---

## 2. Methodology Gap Analysis

### Gap 1: Instrument Offset Model (CRITICAL)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Offset model | Median RV per instrument | Joint Keplerian + per-instrument gamma + jitter | Fundamental |
| CORAVEL included | Yes (16 pts, 280 m/s err) | No (excluded by Nari 2025) | Must exclude |
| HARPS split | No (single "HARPS03"+"HARPS15") | Yes (2015 fiber upgrade creates 17 m/s offset) | Must split |
| ESPRESSO split | Yes (ESPRESSO18/19) | Yes | OK |
| Jitter terms | None | Per-instrument jitter added in quadrature | Missing |
| Offset fitting | Sequential (offset then planets) | Simultaneous (offset + planets in one model) | Fundamental |

**Impact**: This is the root cause of the 8500 m/s residual RMS and the 656x-overestimated
K amplitudes. Fixing this is the single highest priority.

### Gap 2: Orbital Model (CRITICAL)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Planet model | Sinusoidal (A*sin + B*cos) | Keplerian (P, K, e, omega, Tp) | Fundamental |
| Eccentricity | Fixed e=0 (implicit in sinusoid) | Free (Nari: e=0.45 for planet d) | Significant for d |
| N planets | 4 (including disputed 147d) | 3 (18.3, 89.7, 647.6 d) | Wrong planet count |
| Fitting | Sequential (subtract one at a time) | Simultaneous (all planets + offsets jointly) | Fundamental |

**Impact**: Sinusoidal model cannot represent planet d's eccentric orbit (e=0.45).
Even for circular orbits, sequential subtraction introduces cross-talk between signals.

### Gap 3: Activity Modeling (MODERATE)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Magnetic cycle | Not modeled | GP with SHO/MEP kernel, or sinusoid at ~3000d | Missing |
| Rotation | Not modeled | GP quasi-periodic kernel at ~37d | Missing |
| Activity proxies | Not used | S-index, BIS, FWHM as GP inputs or linear decorrelation | Missing |
| Impact | 414d beat artifact appears | Absorbed by GP or decorrelation | Moderate |

**Impact**: The 414-day spurious peak is a direct consequence. With proper activity modeling,
it would be absorbed. However, for a first-pass validation of known planets, the Keplerian
model alone (Gap 2) matters more.

### Gap 4: Detection Statistics (MODERATE)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Criterion | FAP from Lomb-Scargle | False Inclusion Probability (FIP) via nested sampling | Significant |
| Model comparison | None | Bayesian evidence (ln Z) comparing N vs N+1 planet models | Missing |
| Noise model | White noise | Correlated noise (MA, GP, or red noise) | Moderate |

**Impact**: FAP alone cannot distinguish real planets from activity-induced signals.
The 147-day signal (Feng 2017) would likely fail FIP analysis, matching Nari's finding.

### Gap 5: Injection-Recovery Realism (MINOR)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Noise source | Gaussian draws from rv_err | Residuals after planet subtraction | Optimistic |
| Signal model | Sinusoidal | Keplerian (with eccentricity) | Minor |
| Offset structure | Not included | Should include realistic instrument switching | Missing |

**Impact**: Our 90% detection at K=0.05 m/s is unrealistically optimistic.
Nari 2025 reports sensitivity to ~0.30 m/s across orbital periods -- 6x less sensitive
than our estimate, which is more realistic for multi-instrument data.

### Gap 6: Window Function / Alias Analysis (MINOR)

| Aspect | Our pipeline | Literature standard | Gap |
|--------|-------------|--------------------|----|
| Window function | Not computed | Spectral window + Dawson & Fabrycky (2010) alias test | Missing |
| Alias checking | Manual inspection | Systematic: f_alias = f_true +/- n*f_obs | Missing |

**Impact**: We cannot distinguish true periods from aliases. The 652.54-day residual
could be verified/rejected with alias analysis.

---

## 3. Available Tools for Fixing the RV Analysis

### 3.1 RadVel (Recommended Primary Tool)

- **Repo**: github.com/California-Planet-Search/radvel
- **Install**: `pip install radvel`
- **Python**: 3.8-3.12, actively maintained (v1.5.0, Sep 2025)
- **Features**: Joint Keplerian fitting, per-instrument gamma + jitter, GP for activity,
  MCMC via emcee, publication-quality plots
- **Scalability**: Keplerian evaluation is O(N) -- 12K points fine. GP is O(N^3) --
  problematic at 12K, but celerite kernels reduce to O(N).
- **Example**: HD 164922 config -- 2 planets, 3 instruments, exactly our use case

### 3.2 juliet (For Model Comparison)

- **Repo**: github.com/nespinoza/juliet
- **Install**: `pip install juliet`
- **Python**: >=2.7, actively maintained (v2.2.10, Feb 2026)
- **Features**: Wraps RadVel for RV model, adds nested sampling (dynesty) for Bayesian
  evidence. RV-only mode supported.
- **Use case**: Compare 2-planet vs 3-planet vs 4-planet models

### 3.3 kima (Alternative)

- **Repo**: github.com/kima-org/kima
- **Install**: `pip install kima`
- **Features**: Trans-dimensional Bayesian analysis -- automatically determines number
  of Keplerian signals. C++ core for speed. DNest4 sampler.
- **Use case**: Agnostic planet search without pre-specifying N_planets

### 3.4 AliasFinder (Diagnostic)

- **Repo**: github.com/JonasKemmer/AliasFinder
- **Install**: `pip install git+https://github.com/JonasKemmer/AliasFinder.git`
- **Features**: Dawson & Fabrycky (2010) alias analysis, Monte Carlo noise evaluation
- **Use case**: Verify/reject the 652.54-day residual and 414-day peak

---

## 4. Phase 7b.2 Refinement Plan

### Objective

Upgrade the RV analysis to reproduce the Nari et al. (2025) 3-planet solution as
validation, then assess residual signals and true detection sensitivity.

### Step 1: Data Preparation (modify `rv_data.py`)

1. **Exclude CORAVEL**: Filter out CORAVEL-S measurements (16 pts, 280 m/s err).
   These are 2000x worse precision than ESPRESSO and dominate the RMS.
2. **Verify HARPS epoch split**: Confirm HARPS03/HARPS15 correspond to pre/post-2015
   fiber upgrade. DACE labels already separate these.
3. **Verify ESPRESSO epoch split**: Confirm ESPRESSO18/ESPRESSO19 correspond to
   pre/post-June 2019 intervention.
4. **Compute per-instrument summary statistics**: N, median RV, median error, time span,
   to verify instrument labels match expectations.

Expected result: ~12,287 measurements across 4 instruments (H03, H15, E18, E19).

### Step 2: RadVel Joint Keplerian Fit (new `src/rv_keplerian.py`)

1. **Install RadVel**: `pip install radvel`
2. **Implement `fit_keplerian_model()`**:
   - Input: DACE RV data (time, rv, rv_err, instruments)
   - Define 3-planet model: planets at 18.3d, 89.7d, 647.6d
   - Per-instrument free parameters: gamma (offset), jitter
   - Per-planet free parameters: P, Tp, e, omega, K (5 per planet = 15)
   - Total: 15 planet params + 4*2 instrument params = 23 free parameters
   - Use RadVel's built-in MAP optimization followed by MCMC
3. **Validation targets**: Reproduce Nari 2025 values:
   - K_b = 0.614 +/- 0.048 m/s
   - K_c = 0.502 m/s
   - K_d = 0.567 m/s, e = 0.45
   - H03-H15 offset ~ 17 m/s
4. **Compute residual RMS**: Should drop from ~8500 m/s to ~1-2 m/s after joint fit

### Step 3: Residual Analysis (modify `rv_data.py`)

1. **Subtract best-fit 3-planet Keplerian model** (from RadVel)
2. **Run Lomb-Scargle on residuals**: Should NOT show 414-day peak
3. **Check for 147-day signal**: Should be absent or weak (confirming Nari 2025)
4. **Window function analysis**: Compute spectral window to identify alias frequencies
5. **Compare residual periodogram to Nari 2025 Fig. 7**: Their residuals showed
   peaks at 85.6d and 111.7d (mutual 1-year aliases)

### Step 4: Alias Analysis (new function or standalone)

1. Use Dawson & Fabrycky (2010) method or AliasFinder
2. For each residual peak: compute expected alias frequencies from window function
3. Verify 652.54-day signal is an incomplete subtraction of planet d (647.6d, e=0.45)
4. Check if 414-day peak vanishes after proper Keplerian subtraction

### Step 5: Realistic Injection-Recovery (modify `rv_data.py`)

1. **Inject into post-Keplerian residuals** (not raw noise)
2. **Include instrument offset structure**: Synthetic signal samples should reflect
   actual time sampling and instrument switching
3. **Use Keplerian signal model** (not sinusoidal) for injected signals
4. **Compare to Nari 2025 sensitivity**: They report ~0.30 m/s amplitude threshold

### Step 6: Update Known Planets Reference Data

1. **Update `_known_planets_fallback()`** in deep_dive.py:
   - Use Nari 2025 values (3 planets, not 4)
   - Add K amplitudes and eccentricities
   - Flag 147-day as "disputed / not confirmed"
2. **Update NASA query handling**: Cross-reference with Nari 2025 to filter
   disputed signals

### Step 7: Integration Test and Documentation

1. Re-run full deep-dive with upgraded RV analysis
2. Compare results to Nari 2025 Table 3
3. Document methodology gaps that remain vs state-of-the-art (YARARA, GP, FIP)
4. Update `docs/phase7b_deep_dive.md` with 7b.2 results

---

## 5. Verification Criteria for Phase 7b.2

| Criterion | Target | Source |
|-----------|--------|--------|
| Planet b K amplitude | 0.614 +/- 0.10 m/s | Nari 2025 |
| Planet c K amplitude | 0.502 +/- 0.10 m/s | Nari 2025 |
| Planet d K amplitude | 0.567 +/- 0.10 m/s | Nari 2025 |
| Planet d eccentricity | 0.45 +/- 0.15 | Nari 2025 |
| H03-H15 offset | ~17 +/- 5 m/s | Nari 2025 |
| Residual RMS (all) | < 3 m/s | Nari: 0.93 m/s (YARARA), expect ~2 m/s from DRS |
| 414-day peak | Absent from residual periodogram | Beat artifact removed by Keplerian fit |
| 147-day signal | Absent or weak | Not confirmed by Nari 2025 |
| Injection-recovery threshold | K ~ 0.2-0.5 m/s | Nari: ~0.3 m/s |

Note: We use DACE-provided DRS-reduced RVs (not YARARA-reprocessed), so our residual
RMS will be worse than Nari's 0.93 m/s. A target of 1.5-3.0 m/s is realistic.

---

## 6. Methodology Gaps That Will Remain After 7b.2

Even after implementing RadVel-based fitting, our pipeline will still lack:

1. **YARARA/sBART spectral reprocessing**: We use DACE DRS RVs, not YARARA-corrected.
   This is a ~20% precision penalty but acceptable for validation.
2. **GP activity modeling**: Would require GP on 12K points (O(N^3) unless celerite).
   RadVel supports celerite GPs but adds complexity. Could implement as 7b.3.
3. **False Inclusion Probability (FIP)**: State-of-the-art detection statistic from
   nested sampling. juliet provides this but adds runtime. Could implement as 7b.3.
4. **Line-by-line RV analysis**: Requires raw spectra, not available from DACE.
   Out of scope for this project.

---

## 7. Files to Create/Modify

| File | Action | Description |
|------|--------|-------------|
| `src/rv_keplerian.py` | Create | RadVel-based joint Keplerian fitting |
| `src/rv_data.py` | Modify | Add CORAVEL exclusion, improve instrument handling |
| `src/deep_dive.py` | Modify | Use Keplerian fit instead of sinusoidal subtraction |
| `tests/test_phase7b2.py` | Create | Tests for Keplerian fitting |
| `docs/phase7b2_literature_gaps.md` | This file | Literature review and plan |

## 8. Key References

- Pepe et al. 2011, A&A 534, A58 -- discovery (arXiv:1108.3447)
- Feng, Tuomi & Jones 2017, A&A 605, A103 -- TERRA reanalysis (arXiv:1705.05124)
- Cretignier et al. 2023, A&A 678, A2 -- YARARA V2, HZ candidate (arXiv:2308.11812)
- Nari et al. 2025, A&A 693, A297 -- definitive confirmation (arXiv:2501.17092)
- Lovis et al. 2025, A&A -- ESPRESSO 10 cm/s precision (arXiv:2507.07514)
- Kervella et al. 2022, A&A 657, A7 -- PMa catalog (arXiv:2109.10912)
- Fulton et al. 2018, PASP 130, 044504 -- RadVel package
- Dawson & Fabrycky 2010, ApJ 722, 937 -- alias analysis method
- Kennedy et al. 2015, MNRAS 449, 3121 -- debris disk
- Dumusque et al. 2017, A&A 598, A133 -- RV fitting challenge (K/N >= 7.5 threshold)
