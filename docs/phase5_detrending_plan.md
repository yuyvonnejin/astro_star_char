# Phase 5: Detrending Improvements for Shallow Transit Recovery

**Date:** 2026-02-18
**Status:** Planning
**Depends on:** Phase 4 (validation features, complete)

---

## Motivation

Diagnostic analysis of three Kepler systems revealed that the current Savitzky-Golay flatten step is a primary bottleneck for detecting shallow transits in real-world data:

| System | Scatter (ppm) | Best Fold SNR | Deepest Transit (ppm) | Diagnosis |
|--------|--------------|---------------|----------------------|-----------|
| Kepler-410A (KIC 6922244) | 389 | 18.3 | 1400 | Clear detection, good benchmark |
| Kepler-90 (KIC 11497958) | 2824 | 0.0 | 4600 (planet h) | Signal completely missing -- possible data/detrending issue |
| Kepler-186 (KIC 4138008) | 381 | ~marginal | 340 (planet f) | Genuine challenge, signal near noise floor |

Two problems emerged:

1. **The SG filter eats transits.** A fixed window of 1001 cadences can partially or fully absorb transit dips, especially for longer-duration transits. The filter does not know which dips are astrophysical signals vs. instrumental trends.

2. **One window size does not fit all.** Kepler long-cadence (29.4 min) vs. short-cadence (1 min) vs. TESS (2 min / 20 s) produce very different effective time spans for the same `window_length` in cadences. A window that preserves 6-hour transits at 2-min cadence will destroy them at 30-min cadence.

---

## Current Implementation

**File:** `src/lightcurve.py::flatten_lightcurve()` (lines 190-240)

```python
def flatten_lightcurve(lc_data, window_length=1001):
    lc_obj = lk.LightCurve(time=..., flux=..., flux_err=...)
    if window_length % 2 == 0:
        window_length += 1
    # clamp to data length ...
    lc_flat = lc_obj.flatten(window_length=window_length)
    out["flux_flat"] = lc_flat.flux.value
```

**Call chain:**
```
retrieve_lightcurve(target, window_length=1001)
  -> flatten_lightcurve(lc_data, window_length=1001)
     -> lightkurve LightCurve.flatten(window_length=1001)   # Savitzky-Golay
```

**Consumers of flattened flux:**
- `src/periodogram.py` (Lomb-Scargle): `flux = lc_data.get("flux_flat", lc_data["flux"])`
- `src/transit.py` (BLS detection): `flux = lc_data.get("flux_flat", lc_data["flux"])`

**Default parameters:**
| Parameter | Value | Meaning |
|-----------|-------|---------|
| `window_length` | 1001 | SG filter window in cadences |
| `bin_cadence_s` | 120.0 | Pre-flatten binning cadence |

At 120 s binned cadence, 1001 points ~ 1.4 days. At Kepler long-cadence (29.4 min) unbinned, 1001 points ~ 20.5 days. The binning step normalizes high-cadence data but leaves Kepler long-cadence untouched, creating a large asymmetry.

---

## Feature 1: Sigma-Mask Flatten

### Problem

Standard SG flatten fits a polynomial through ALL data points, including transit dips. The fit "bows" toward the dips, partially filling them in. For shallow transits (< 500 ppm), this can reduce depth below the BLS detection threshold.

### Solution

**Sigma-mask (iterative sigma-clipping) flatten:**
1. Run an initial SG flatten pass.
2. Identify outlier points that deviate by more than `N` sigma below the trend (these are likely transit dips or flare artifacts).
3. Mask those points (exclude from the polynomial fit).
4. Re-run the SG flatten on the masked data.
5. Apply the resulting trend to ALL points (including the masked ones).

This is a well-established technique. lightkurve's `flatten()` already supports it via the `sigma` parameter:

```python
lc_flat = lc_obj.flatten(window_length=window_length, sigma=3.0)
```

When `sigma` is set, lightkurve internally runs iterative sigma clipping before the SG fit. Points more than `sigma` standard deviations below the trend are excluded from the fit (but the trend is still applied to divide them out).

### Implementation

**File: `src/lightcurve.py`**

Modify `flatten_lightcurve()` signature:
```python
def flatten_lightcurve(lc_data, window_length=1001, sigma_mask=None):
```

- `sigma_mask=None` (default): current behavior, no masking.
- `sigma_mask=3.0` (recommended for transit work): mask 3-sigma outliers before fitting.

Pass through from `retrieve_lightcurve()`:
```python
def retrieve_lightcurve(target, ..., sigma_mask=None):
    ...
    lc_data = flatten_lightcurve(lc_data, window_length=window_length, sigma_mask=sigma_mask)
```

lightkurve call becomes:
```python
if sigma_mask is not None:
    lc_flat = lc_obj.flatten(window_length=window_length, sigma=sigma_mask)
else:
    lc_flat = lc_obj.flatten(window_length=window_length)
```

Log the masking:
```python
logger.info("Flattened with window_length=%d, sigma_mask=%s", window_length,
            sigma_mask if sigma_mask else "off")
```

Store in metadata:
```python
out["sigma_mask"] = sigma_mask
```

**File: `src/pipeline.py`**

Add `lc_sigma_mask=None` to `process_star()` and `_run_lightcurve_analysis()`. Pass through to `retrieve_lightcurve()`.

**File: `run_stars.py`**

Add CLI arg:
```python
parser.add_argument("--sigma-mask", type=float, default=None,
                    help="Sigma threshold for masking transit dips during flatten (e.g. 3.0)")
```

Pass through `run()` to `process_star(lc_sigma_mask=...)`.

### Testing

**File: `tests/test_lightcurve.py`** -- add `TestSigmaMaskFlatten`:

1. **`test_sigma_mask_preserves_transit_dip`**: Inject a box transit into synthetic data with a slow trend. Flatten with and without sigma_mask. Verify that sigma_mask=3.0 preserves more depth than no masking.

2. **`test_sigma_mask_none_is_default_behavior`**: Verify that `sigma_mask=None` produces identical output to the original code.

3. **`test_sigma_mask_metadata`**: Verify `lc_data["sigma_mask"]` is stored correctly.

4. **`test_sigma_mask_with_short_data`**: Edge case -- fewer points than window.

### Expected Impact

- **Kepler-186f** (340 ppm transit, 381 ppm scatter): sigma-mask should recover ~50-100 ppm of depth that the current flatten absorbs, potentially pushing fold SNR above detection threshold.
- **Kepler-90** planets: if the data quality issue is resolved, sigma-mask protects the 200-4600 ppm transits during detrending.
- **No regression risk**: `sigma_mask=None` default preserves current behavior for all existing users.

### Recommended Default

For transit-focused runs (`--transit` flag), consider auto-enabling `sigma_mask=3.0` unless explicitly overridden. This can be done in `process_star()`:
```python
if include_transit and lc_sigma_mask is None:
    lc_sigma_mask = 3.0  # auto-enable for transit detection
    logger.info("Auto-enabling sigma_mask=3.0 for transit detection")
```

---

## Feature 2: Configurable Flatten Window

### Problem

The hardcoded `window_length=1001` is a compromise that works reasonably well for binned 120 s data (~1.4 day window) but is suboptimal in several scenarios:

| Scenario | Cadence | Effective Window | Problem |
|----------|---------|-----------------|---------|
| Kepler long-cadence, unbinned | 29.4 min | ~20.5 days | Too wide -- doesn't remove instrumental trends |
| Kepler short-cadence, binned to 120 s | 120 s | ~1.4 days | OK for most transits |
| TESS 2-min, binned to 120 s | 120 s | ~1.4 days | OK for most transits |
| TESS 20 s, unbinned | 20 s | ~5.6 hours | Too narrow -- eats transits > 3 hours |
| Long-period HZ planets | any | any | Transit duration 6-13 hours, needs wider window |

The fix is to either:
- **(a)** Expose `window_length` as a tunable parameter (simple), or
- **(b)** Compute `window_length` dynamically from cadence and expected transit duration (smarter but more complex).

We implement both: expose the parameter for manual control, and add an "auto" mode that adapts to the data.

### Implementation

**File: `src/lightcurve.py`**

Modify `flatten_lightcurve()` and `retrieve_lightcurve()`:

```python
def flatten_lightcurve(lc_data, window_length="auto", sigma_mask=None):
    ...
    if window_length == "auto":
        window_length = _compute_auto_window(lc_data)
    ...
```

Add helper:
```python
def _compute_auto_window(lc_data):
    """Compute SG window length from data cadence.

    Target: window covers ~2 days of data (good for planets up to ~15 day period).
    Minimum: 101 cadences. Maximum: 2001 cadences.
    """
    time = lc_data["time"]
    if len(time) < 10:
        return 101  # minimum safe window

    # Estimate median cadence in days
    dt = np.median(np.diff(time))
    cadence_s = dt * 86400.0

    # Target window duration: 2.0 days
    target_duration_s = 2.0 * 86400.0
    window = int(target_duration_s / cadence_s)

    # Ensure odd
    if window % 2 == 0:
        window += 1

    # Clamp to [101, 2001]
    window = max(101, min(2001, window))

    logger.info("Auto window_length=%d (cadence=%.1f s, target=2.0 d)", window, cadence_s)
    return window
```

**Rationale for 2-day target:**
- A 2-day SG window removes trends slower than ~1 day.
- Planet transits last 1-13 hours (well within the window).
- Even long-duration HZ transits (~13 hours for Earth analogs) are < 1 day.
- Removes Kepler quarterly systematics (days-scale ramps) effectively.

**File: `src/pipeline.py`**

Add `lc_window_length="auto"` to `process_star()` and pass through.

**File: `run_stars.py`**

Add CLI arg:
```python
parser.add_argument("--flatten-window", type=str, default="auto",
                    help="SG flatten window length in cadences, or 'auto' (default: auto)")
```

Parse as int or "auto":
```python
window = args.flatten_window
if window != "auto":
    window = int(window)
```

### Testing

**File: `tests/test_lightcurve.py`** -- add `TestAutoWindow`:

1. **`test_auto_window_kepler_longcadence`**: Synthetic data at 29.4 min cadence. Auto window should be ~98 cadences (2 days / 29.4 min ~ 98), clamped and made odd -> ~99.

2. **`test_auto_window_tess_2min`**: Synthetic data at 120 s cadence. Auto window should be ~1441 cadences.

3. **`test_auto_window_tess_20s`**: Synthetic data at 20 s cadence. Auto window should be capped at 2001.

4. **`test_auto_window_short_data`**: Fewer than 10 points -> returns 101.

5. **`test_explicit_window_overrides_auto`**: Passing `window_length=501` should use 501, not auto.

6. **`test_auto_preserves_existing_behavior`**: For 120 s cadence data, auto window ~ 1441 (close to current 1001 but wider -- verify no regression on existing test targets).

### Migration

- Change default from `1001` to `"auto"` in `flatten_lightcurve()` and `retrieve_lightcurve()`.
- Existing tests that explicitly pass `window_length=1001` continue to work.
- Tests that rely on the default will get the auto-computed value instead. Run full test suite to verify no regressions.

---

## Files Modified (Summary)

| File | Changes |
|------|---------|
| `src/lightcurve.py` | Modify `flatten_lightcurve()` and `retrieve_lightcurve()` signatures; add `_compute_auto_window()` |
| `src/pipeline.py` | Add `lc_sigma_mask`, `lc_window_length` params to `process_star()` and `_run_lightcurve_analysis()` |
| `run_stars.py` | Add `--sigma-mask`, `--flatten-window` CLI args; pass through |
| `tests/test_lightcurve.py` | Add `TestSigmaMaskFlatten` (4 tests), `TestAutoWindow` (6 tests) |

## Implementation Order

1. **Sigma-mask flatten** -- highest impact, smallest code change (leverages existing lightkurve `sigma` parameter).
2. **Configurable flatten window** -- broader improvement, requires auto-compute helper and more testing.
3. **Integration test** -- re-run diagnostic notebook on all three Kepler targets with new defaults (`sigma_mask=3.0`, `window_length="auto"`).

## Verification Plan

1. Unit tests: `.\venv\Scripts\python -m pytest tests/test_lightcurve.py -v -k "SigmaMask or AutoWindow"`
2. Full regression: `.\venv\Scripts\python -m pytest tests/ -v` -- all existing tests pass
3. Integration: re-run diagnostic transit notebook for KIC 6922244, KIC 11497958, KIC 4138008 with `--sigma-mask 3.0`
4. Compare fold SNR and transit depth recovery before/after

## Open Questions

1. **Should sigma_mask auto-enable when `--transit` is used?** Proposed yes, with opt-out via `--sigma-mask 0` (disabled). Needs user input.
2. **Should auto window target be configurable?** Currently hardcoded to 2.0 days. Could expose as `--flatten-target-days` but this may be over-engineering. Start with fixed 2.0 day target.
3. **Kepler-90 data quality**: Before attributing all failures to detrending, need to verify that the raw data is actually present and not corrupted. A separate diagnostic step should check total number of quarters, cadence, and raw scatter before flatten.
