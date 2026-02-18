# Phase 4: Transit Validation & HZ-Targeted Detection

## Problem Statement

Transit detection works on textbook cases (KIC 6922244) but fails on real-world
targets (KIC 11497958, KIC 4138008). Diagnostic analysis shows signals are lost
before BLS runs -- variable stars and aggressive detrending drown shallow transits.

Rather than fixing detrending for all stars, the smarter approach is:
1. Focus on quiet Sun-like stars where shallow transits are detectable
2. Add validation to distinguish real planets from false positives

## Feature 1: Quiet-Star Filter (Gate Module 5)

**Purpose:** Skip transit detection on highly variable stars where BLS cannot
reliably find shallow transits. Saves compute and avoids false positives.

**Implementation:**
- `_check_quiet_star(result, max_amplitude_ppt)` in `pipeline.py`
- Uses `amplitude_ppt` from Module 4 periodogram analysis
- Default threshold: 10 ppt (~Jupiter-transit depth around Sun)
- Configurable via `max_variability_ppt` parameter
- Override with `force_transit=True` for research/testing

**Behavior:**
- Quiet star (amplitude < threshold): proceed to Module 5
- Variable star (amplitude >= threshold): skip Module 5, set `planet_flag="too_variable"`
- Missing amplitude data: proceed (conservative -- do not block)

## Feature 2: HZ-Targeted BLS Mode

**Purpose:** Narrow the BLS period search to the habitable zone, improving
sensitivity for the most scientifically interesting planets.

**Implementation:**
- `compute_hz_period_range(mass_Msun, luminosity_Lsun, broadening_factor)` in `transit.py`
- Converts optimistic HZ boundaries (AU) to orbital periods via Kepler's 3rd law
- Broadens range by configurable factor (default 2.0x) to account for stellar
  property uncertainties (~20-30%)
- Floor at 0.5 days minimum period

**Behavior:**
- When `hz_targeted=True` and stellar properties available: narrow period range
- When stellar properties missing: fall back to full [0.5, 100] day range
- CLI: `--hz-only` flag, `--hz-broadening` parameter

## Feature 3: Even/Odd Transit Validation

**Purpose:** Detect eclipsing binaries by comparing transit depths in even vs
odd numbered transits. Real planets produce equal depths; EBs often show
different primary/secondary eclipse depths.

**Implementation:**
- `validate_even_odd(time, flux, period, epoch, duration_days)` in `transit.py`
- Assigns epoch numbers, splits in-transit points by even/odd
- Measures median depth independently
- Pass criterion: depth ratio within [0.5, 2.0]

**Design decisions:**
- Generous tolerance [0.5, 2.0] because shallow transits are noisy
- EBs typically show ratios of 3-10x
- Reports: `depth_even_ppm`, `depth_odd_ppm`, `depth_ratio_even_odd`, `validation_pass`

## Feature 4: Transit Shape Test (V vs U)

**Purpose:** Classify transit shape as U-shaped (planet-like flat bottom) or
V-shaped (eclipsing binary grazing geometry). U-shaped transits have a flat
ingress-egress-bottom pattern; V-shaped are triangular.

**Implementation:**
- `classify_transit_shape(time, flux, period, epoch, duration_days, n_phase_bins)` in `transit.py`
- Phase-folds and bins finely within 1.5x transit duration window
- Measures flat-bottom fraction: bins within 20% of minimum / total valid bins
- Classification: U_shape (>0.3), V_shape (<0.15), ambiguous (between)

**Design decisions:**
- Thresholds based on typical transit shapes in Kepler data
- U-shape suggests planet with well-defined ingress/egress
- V-shape suggests grazing eclipse or blended EB

## Files Modified

| File | Changes |
|------|---------|
| `src/transit.py` | +3 functions, modify `analyze_transit` signature |
| `src/pipeline.py` | +1 function, modify `process_star` and `_run_transit_analysis` |
| `run_stars.py` | +4 CLI args, update `run()` and `format_result()` |
| `tests/test_transit.py` | +3 test classes (17 new tests) |
| `tests/test_pipeline.py` | +1 test class (4 new tests) |
| `docs/detailed_pipeline_gaia_based.md` | Add Module 5c section |

## Testing Strategy

1. Unit tests for each new function with synthetic data
2. Existing 96 tests must continue passing (no regressions)
3. Integration test: `run_stars.py --name "KIC 6922244" --transit --hz-only`
