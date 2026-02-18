# Transit Detection Diagnostic Guide

When transit detection works on easy targets (deep, single-planet, quiet star) but fails on harder cases, use this guide to identify where the signal is lost.

## Root Causes

Transit detection can fail at three stages:

1. **Detrending** -- the flattening step absorbs the transit dip into the baseline, or leaves stellar variability that drowns the signal.
2. **Noise floor** -- the transit depth is below the photometric precision of the data. No algorithm can recover it.
3. **BLS configuration** -- the period grid, duration grid, or SDE threshold misses a real signal that is present in the data.

The diagnostic notebook (`diagnostic_transit.ipynb`) tests all three systematically.

## Diagnostic Steps

### Step 1-2: Raw Light Curve Inspection

Download all available Kepler/TESS quarters and plot the full time series.

**What to look for:**
- Overall variability amplitude (ppm) -- higher variability means harder transit detection
- Gaps between quarters/sectors -- affects the effective baseline for long-period planets
- Obvious systematics (jumps, slopes, outliers)

### Step 3: Outlier Clipping Check

Compare the 5-sigma clipping threshold against the expected transit depths.

**Key question:** Is the clipping threshold larger than the transit depth?

- If `5 * std > transit_depth`: transits survive clipping (good)
- If `5 * std < transit_depth`: transits are being removed as outliers (problem)

For most targets this is fine because transit depths (tens to thousands of ppm) are much smaller than 5-sigma of the overall scatter. But on very quiet stars with deep transits, check this.

### Step 4: Flatten Window Sweep

Flatten with multiple window sizes and measure residual scatter.

| Window (cadences) | Window (days, 30-min Kepler) | Use case |
|---|---|---|
| 101 | ~2.1 d | Aggressive -- removes most variability but may eat transits > few hours |
| 301 | ~6.3 d | Moderate |
| 501 | ~10.4 d | Default-ish |
| 1001 | ~20.8 d | Pipeline default -- good for periods up to ~10 d |
| 2001 | ~41.7 d | Conservative -- preserves long-term signals |

**Rule of thumb:** the window should be at least 3x the transit duration. For a typical transit duration of 2-6 hours, even the smallest window (101) is safe. The risk is that shorter windows also remove stellar variability more aggressively, which reduces scatter and improves SNR.

### Step 5: Phase-Fold on Known Periods (Key Diagnostic)

Fold the flattened light curve at each known planet's period, bin into ~200 phase bins, and measure the depth + SNR of the deepest dip.

**Interpretation:**

| Fold SNR | Status | Meaning |
|---|---|---|
| > 3.0 | DETECTED | Signal is visible -- BLS should find it |
| 1.5 - 3.0 | MARGINAL | Signal is borderline -- may need better detrending |
| < 1.5 | NOT VISIBLE | Signal lost before BLS runs |

If a planet is NOT VISIBLE even when folding on the known period, BLS cannot find it. The problem is upstream (detrending or noise floor), not in BLS.

### Step 6: Window Sensitivity Matrix

Cross-reference: deepest planet from each target x all window sizes. Shows whether any window size recovers the signal.

If no window works -> the problem is fundamental (noise floor or detrending approach, not window size).

### Step 7: Sigma-Clipping Flatten (Transit-Aware Detrending)

Compare three flattening methods:

1. **Pipeline default** (`flatten(window_length=1001)`) -- no transit masking
2. **Sigma=3 mask** (`flatten(window_length=1001, sigma=3)`) -- masks 3-sigma outliers during trend fit
3. **Sigma=2 mask** (`flatten(window_length=1001, sigma=2)`) -- more aggressive masking

**How sigma-masking works:** During the Savitzky-Golay trend fit, points that dip below `sigma` standard deviations are excluded from the fit. This prevents the trend line from being pulled down by transit dips, preserving the true transit depth in the residuals.

**If sigma-mask flatten recovers a signal that the default doesn't**, the fix is to add sigma masking to the pipeline's `flatten_lightcurve()` function.

### Step 8: Noise Floor Analysis

Compute the theoretical detection limit for each planet:

```
expected_SNR = depth_ppm * sqrt(n_transits * n_points_per_transit) / scatter_ppm
```

Where:
- `scatter_ppm` = point-to-point scatter of the flattened light curve
- `n_transits` = baseline / period
- `n_points_per_transit` = transit_duration / cadence

**If expected_SNR < 3**, the planet is not detectable with this data quality regardless of algorithm. Focus efforts on planets with expected_SNR > 6.

### Step 9: BLS Power Spectrum

Run BLS and overlay known planet periods on the power spectrum.

**What to look for:**
- Is there a peak at or near the known period? (Even sub-threshold peaks indicate the signal is present but weak)
- Are there dominant peaks at unrelated periods? (Systematics or stellar variability leaking through)
- Are there peaks at harmonics (P/2, P/3) of the known period? (BLS harmonic confusion)

### Step 10: Summary Table

Automated diagnosis per planet combining all metrics.

## Decision Tree

After running the diagnostic notebook:

```
Phase-fold on known period shows transit?
  |
  +-- YES, SNR > 3 --> BLS config issue
  |     Check: period grid resolution, SDE threshold, duration range
  |
  +-- NO --> Compare sigma-mask vs default flatten
        |
        +-- Sigma-mask recovers it --> Implement sigma-mask in pipeline
        |
        +-- Still not visible --> Check noise floor
              |
              +-- expected_SNR > 6 --> Detrending approach fundamentally wrong
              |     Consider: GP detrending, iterative clipping, per-quarter flatten
              |
              +-- expected_SNR < 3 --> Below noise floor
                    Not recoverable with current data/method
```

## Common Fixes

| Diagnosis | Fix |
|---|---|
| Signal absorbed by detrending | Add `sigma` parameter to `flatten()` call |
| Window too small eats transit | Increase `window_length` |
| Window too large leaves variability | Decrease `window_length` or pre-whiten variability |
| Depth refinement corrupts signal | Disable masking when period precision is low (> 1% error) |
| Period slightly off -> bad fold | Add fine grid search around BLS peak before depth measurement |
| Below noise floor | Accept limitation; would need better data or stacking techniques |

## Reference: Known Test Systems

### KIC 6922244 (Kepler-410A) -- Easy

Single planet, 490 ppm depth, G-type star. Pipeline handles this correctly.

### KIC 11497958 (Kepler-90) -- Hard

8 planets ranging from 51 to 4600 ppm depth around a G0V star (R ~ 1.2 Rsun). The two deepest (Kepler-90g at 2370 ppm, Kepler-90h at 4600 ppm) should be detectable. The inner planets (b, c, i at 50-63 ppm) are near or below the noise floor.

### KIC 4138008 (Kepler-186) -- Hard

5 planets around an M1 dwarf (R ~ 0.47 Rsun). Small star means relatively deep transits (340-585 ppm) for Earth-sized planets. Kepler-186f (129.9 d period, 408 ppm) is the famous HZ planet.
