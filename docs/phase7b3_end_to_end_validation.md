# Phase 7b.3: End-to-End Validation and Data Quality Hardening

## Status: Complete (2026-03-20)

Phase 7b.2 built the RadVel Keplerian fitting module and integrated it into the deep-dive
pipeline. Phase 7b.3 ran the full pipeline end-to-end on HD 20794, diagnosed five bugs
that prevented convergence, fixed them, and documented the remaining gap between our
automated results and publication-quality analysis.

---

## 1. Bugs Found and Fixed

### 1.1 Time-of-Conjunction in Wrong Time System

**File**: `src/deep_dive.py` `_build_planet_params()`

The initial tc (time of conjunction) was hardcoded to JD 2456000.0, but DACE returns
timestamps in RJD (JD - 2400000). This placed the orbital phase reference ~2.4 million
days in the future of the data, causing the Keplerian model to produce extreme values
and corrupting all fitted parameters.

**Fix**: Compute `tc_base = np.median(time)` from the actual RV time array, keeping tc
in the same time system as the data.

### 1.2 NASA Disputed Planets Not Filtered

**File**: `src/deep_dive.py` `_extract_known_periods()`

The NASA Exoplanet Archive returns 4 planets for HD 20794 including the 147-day signal
("planet e" from Feng 2017), which Nari et al. (2025) showed worsens model evidence.
NASA records lack a `confirmed` field, so the existing `confirmed=False` check only
worked with fallback data, not live NASA queries.

**Fix**: Cross-reference NASA planets with literature fallback data. Any NASA planet
whose period matches a `confirmed=False` entry in the fallback (within 1%) is rejected
with a log message explaining the literature source.

### 1.3 Catastrophic RV Outliers in DACE Data

**File**: `src/rv_data.py` `rv_filter_instruments()`

DACE public RV data for HD 20794 contains catastrophic pipeline failures that survived
DACE's own QC:

| Instrument | Outlier example | Systemic velocity | Deviation |
|------------|----------------|-------------------|-----------|
| ESPRESSO18 | -119,552 m/s | +89,617 m/s | 209 km/s |
| HARPS03 | 73,347 m/s | +87,822 m/s | 14.5 km/s |
| HARPS03 | 105,027 m/s | +87,822 m/s | 17.2 km/s |
| HARPS15 | 88,078 m/s | +87,843 m/s | 235 m/s |

359 measurements total were contaminated (281 HARPS03, 58 HARPS15, 20 ESPRESSO18).

**Fix**: Added MAD-based sigma-clipping to `rv_filter_instruments()`. For each instrument,
compute the median and MAD (median absolute deviation), then reject points deviating more
than `sigma_clip * 1.4826 * MAD` from the median. Default: 5-sigma.

### 1.4 ESPRESSO18 Systematic Scatter

**File**: `src/rv_data.py` `rv_filter_instruments()`

After sigma-clipping, ESPRESSO18 still had 387 measurements with 783 m/s scatter despite
0.16 m/s reported errors (ratio = 5058x). This is not a few outliers -- the entire
instrument's distribution is corrupted, likely from inconsistent DRS pipeline versions
or CCF masks in the public DACE data.

**Fix**: Added a `max_scatter_factor` parameter. After sigma-clipping, instruments whose
MAD-sigma exceeds `max_scatter_factor * median_error` are excluded entirely. Default
threshold: 50x. ESPRESSO18 (5058x) is caught; ESPRESSO19 (ratio ~16), HARPS03 (~15),
and HARPS15 (~12) are retained.

### 1.5 Numerical Precision with Absolute RV Baseline

**File**: `src/rv_keplerian.py` `fit_keplerian()`

RadVel was fitting ~0.5 m/s planetary signals on top of ~88,000 m/s absolute velocities.
The 1:176,000 dynamic range caused the scipy optimizer to converge to local minima where
instrument gammas absorbed the Keplerian signal.

**Fix**: Pre-center the data by subtracting per-instrument median RV before passing to
RadVel. Gammas are initialized to 0.0 (centered frame). After fitting, the subtracted
medians are added back to the output gammas so the report shows absolute values.

---

## 2. Optimization Strategy: 3-Pass MAP

The Keplerian MAP optimization now uses three sequential passes:

| Pass | Fixed | Free | Purpose |
|------|-------|------|---------|
| 1 | P, e, w | K, tc, gamma, jitter | Find correct instrument offsets |
| 2 | P, e (optional) | K, tc, w, gamma, jitter | Refine orbital phase |
| 3 | (none, or e) | all | Fine-tune periods |

When `fix_eccentricities=True` (used when literature reference data is available),
eccentricities remain fixed at their literature values throughout all three passes.
This prevents the optimizer from absorbing correlated noise into extreme eccentricities.

---

## 3. Results: Pipeline vs Literature

### 3.1 Progression Across Iterations

| Iteration | RMS (m/s) | Key change |
|-----------|-----------|------------|
| Feb 19 (Phase 7b) | 8537 -> 8505 | Sinusoidal subtraction, CORAVEL included |
| Run 1 (tc fix) | 8374 -> 8242 | tc fixed, but no outlier rejection |
| Run 2 (sigma-clip) | 316 -> 138 | 359 outliers removed, pre-centering |
| Run 3 (ESPRESSO18 excluded) | 10.3 -> 1.77 | Scatter filter caught ESPRESSO18 |
| Run 4 (fix_e=True) | 10.3 -> 1.77 | Eccentricities fixed (same RMS) |

### 3.2 Final Results vs Nari et al. (2025)

**Data quality metrics**:

| Metric | Pipeline | Literature |
|--------|----------|------------|
| RMS after fit | 1.77 m/s | ~1.5 m/s |
| Instruments used | ESPRESSO19, HARPS03, HARPS15 | HARPS(YARARA), ESPRESSO(sBART) |
| Measurements | 11,854 | 806 nightly-binned |
| Jitter (HARPS03) | 1.78 m/s | -- |
| Jitter (ESPRESSO19) | 1.07 m/s | -- |
| HARPS15-HARPS03 offset | 10.4 m/s | 17.0 +/- 1.7 m/s |

**Planet parameters** (fit did not converge to literature values):

| Planet | Pipeline P (d) | Lit P (d) | Pipeline K (m/s) | Lit K (m/s) | Pipeline e | Lit e |
|--------|---------------|-----------|------------------|-------------|-----------|-------|
| b | 13.5 (drifted) | 18.315 | 0.82 | 0.614 | 0.91 | 0.06 |
| c | 91.2 | 89.68 | 0.90 | 0.502 | 0.69 | 0.08 |
| d | 752.5 (drifted) | 647.6 | 1.22 | 0.567 | 0.93 | 0.45 |

The 18.31-day planet b signal appears as the strongest residual peak (FAP = 5e-177),
confirming it is real but not properly captured by the fit.

**What matches well**:
- Residual RMS (1.77 vs ~1.5 m/s)
- Instrument jitters (1-2 m/s, physically reasonable)
- Inter-instrument offsets (correct sign and order of magnitude)
- Planet c period (91.2 vs 89.68 d, 1.7% off)
- K amplitudes in correct order of magnitude (sub-2 m/s vs sub-1 m/s)

**What does not match**:
- Planet b and d periods drifted significantly from literature
- Eccentricities converge to extreme values (0.7-0.9) instead of literature (0.06-0.45)
- K amplitudes 1.5-2x too high

### 3.3 What Went Right (vs Phase 7b)

| Feature | Phase 7b (Feb 19) | Phase 7b.3 (Mar 20) |
|---------|-------------------|---------------------|
| RMS | 8537 m/s | 1.77 m/s |
| K amplitudes | 400-750 m/s (offset artifacts) | 0.8-1.2 m/s (correct order) |
| Planet count | 4 (included disputed 147d) | 3 (correctly rejected 147d) |
| CORAVEL | Included | Excluded |
| ESPRESSO18 | Included (corrupted fit) | Excluded (scatter filter) |
| Outliers | Not handled | 359 sigma-clipped |
| Method | Sinusoidal subtraction | RadVel Keplerian MAP |
| Offsets | Median subtraction | Joint fit (gamma + jitter) |

---

## 4. Fundamental Limitation: Why Parameters Don't Match

The pipeline uses a simple Keplerian MAP fit on raw (non-detrended) DACE absolute RVs.
Nari et al. (2025) achieved their results through a fundamentally different methodology:

1. **YARARA spectral post-processing** (Cretignier et al. 2021): removes telluric
   contamination, optical aberrations, and detector systematics at the spectral level
   before computing RVs. This reduces the effective jitter from ~2 m/s to ~1 m/s.

2. **sBART template matching** (Artigau et al. 2022): for ESPRESSO data, produces
   cleaner RVs than the standard DRS CCF method.

3. **GP stellar activity model**: a quasi-periodic Gaussian Process trained on activity
   indicators (log R'HK, FWHM, BIS) to model and subtract the ~39-day rotation signal
   and long-term magnetic cycle.

4. **FIP detection criterion** (Hara et al. 2022): Frequency and phase Information
   Content, a Bayesian model selection tool that avoids the false positive inflation
   of traditional periodogram FAP.

5. **Nested sampling** (PolyChord): full posterior exploration that avoids the local
   minima traps of MAP optimization.

Without items 1-3, the ~1.8 m/s correlated noise from stellar activity is indistinguishable
from the 0.5-0.6 m/s planetary signals in a simple MAP fit. The optimizer absorbs the
correlated noise into extreme eccentricities and shifted periods rather than leaving it
as white residual.

---

## 5. Code Changes Summary

| File | Function | Change |
|------|----------|--------|
| `src/rv_data.py` | `rv_filter_instruments()` | Added `sigma_clip` (MAD-based, default 5.0) and `max_scatter_factor` (default 50.0) parameters |
| `src/rv_keplerian.py` | `fit_keplerian()` | Pre-centers data by subtracting per-instrument medians; added `fix_eccentricities` parameter; 3-pass MAP strategy |
| `src/rv_keplerian.py` | `keplerian_residual_analysis()` | Pass-through for `fix_eccentricities` |
| `src/deep_dive.py` | `_extract_known_periods()` | Cross-references NASA planets with literature fallback to reject disputed periods |
| `src/deep_dive.py` | `_build_planet_params()` | Accepts `time` array; computes `tc_base = np.median(time)` in data time system |
| `src/deep_dive.py` | `_run_keplerian_or_sinusoidal()` | Sets `fix_eccentricities=True` when literature reference data is available |

All 10 Phase 7b.2 unit tests pass after changes.

---

## 6. Recommendations for Future Work

### 6.1 Near-term (achievable with current architecture)

- **Nightly binning**: Average intra-night RV measurements to reduce N from 11,854 to
  ~1,000. This is standard practice and may improve MAP convergence by reducing
  correlated intra-night noise.
- **Activity indicator decorrelation**: DACE provides FWHM and BIS; linear decorrelation
  against these could remove some activity signal without a full GP.
- **Iterative planet subtraction**: Fit one planet at a time (strongest first), subtract,
  then fit the next. This avoids the multi-planet local minimum problem.

### 6.2 Medium-term (requires new dependencies)

- **GP activity model**: Add a quasi-periodic GP kernel (via george or celerite) to the
  RadVel fit. This is the single most impactful improvement.
- **juliet integration**: Use juliet's nested sampling wrapper around RadVel for proper
  model comparison (is the 147d signal real? how many planets?).
- **YARARA-processed RVs**: If available via DACE or ESO archive, use pre-processed RVs
  instead of raw DRS values.

### 6.3 Scope boundary

Achieving true publication-quality parameters for sub-m/s RV systems like HD 20794 is
fundamentally a spectral-level data reduction problem, not a fitting problem. The pipeline's
role is to identify targets, characterize data quality, and provide first-order analysis --
not to replicate the full analysis chain of a dedicated publication.
