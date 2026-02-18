# Phase 4: Transit Validation & HZ-Targeted Detection

## Problem Statement

Transit detection works on textbook cases (KIC 6922244) but fails on real-world
targets (KIC 11497958, KIC 4138008). Diagnostic analysis shows signals are lost
before BLS runs -- variable stars and aggressive detrending drown shallow transits.

Rather than fixing detrending for all stars, the smarter approach is:
1. Focus on quiet Sun-like stars where shallow transits are detectable
2. Add validation to distinguish real planets from false positives

## Scientific Rationale: Starting with Sun-Like Stars

The pipeline's transit detection strategy starts with **quiet, Sun-like stars** as
the primary science case. This is a deliberate design choice, not a limitation:

- **Signal detectability:** A Sun-like star (Teff ~5000-6000 K, amplitude < 10 ppt)
  provides the cleanest photometric baseline for detecting shallow transit dips.
  Earth-sized planets produce ~84 ppm transits around Sun-like stars -- detectable
  with Kepler/TESS data, but only if the stellar noise floor is well below this.

- **Habitable zone relevance:** The habitable zone around Sun-like stars (G/K dwarfs)
  corresponds to orbital periods of ~100-600 days, well within Kepler's baseline.
  Targeting HZ periods on quiet stars is the most efficient path to finding
  potentially habitable planets.

- **False positive rejection:** Eclipsing binaries (EBs) are the dominant contaminant
  in transit surveys. EBs tend to be variable and have deep, V-shaped eclipses.
  The quiet-star filter removes most EBs before BLS runs, and the validation
  checks (even/odd depth, transit shape) catch those that remain.

- **Practical scope:** Noisy stars (active M dwarfs, pulsators, spotted giants)
  require specialized detrending techniques (Gaussian processes, pixel-level
  decorrelation) that are beyond the current pipeline. Deferring these to Phase 5+
  keeps the pipeline focused and reliable.

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

## Feature 5: Pre-Whitening (Variability Removal)

**Purpose:** Remove periodic stellar variability (starspot modulation, pulsations)
before BLS to prevent it from masking or mimicking transit signals.

**Implementation:**
- `remove_variability(time, flux, period, n_bins=50)` in `transit.py`
- Uses model-free phase-folded template subtraction
- The variability period comes from Module 4 Lomb-Scargle detection
- Only applied when `variability_class="periodic"` in Module 4 output

**Algorithm:**
1. Phase-fold light curve at detected variability period
2. Bin into 50 phase bins, compute median flux per bin
3. Interpolate template at each data point's phase (with wrap-around)
4. Subtract template and re-center on original median

**Integration:** Called in `analyze_transit()` before BLS when `variability_period`
is provided. Controlled by `pipeline.py` which passes the Module 4 detected period.

## Feature 6: Multi-Candidate BLS with Stratified Extraction

**Purpose:** Extract multiple transit candidates across different period ranges,
preventing dominant short-period systematics from suppressing real signals at
longer periods.

**Implementation:**
- `_extract_bls_candidates_stratified()` in `transit.py`
- Log-uniform period grid (50,000 points) for balanced sensitivity
- Divides period range into ~5 logarithmic bins (3 per decade, min 3)
- Per-bin (local) SDE prevents cross-contamination between period regimes
- Returns one candidate per period bin, sorted by descending local SDE

**Key innovation -- Local SDE:** Standard global SDE compares peaks to mean/std of
the full power spectrum. A dominant short-period noise peak inflates the global std,
making real signals at longer periods appear insignificant. Local SDE computes
significance within each period bin only, ensuring each period decade is evaluated
on its own noise floor.

## Feature 7: Depth Refinement for Multi-Candidate Systems

**Purpose:** In multi-planet systems, BLS depth for each candidate is contaminated
by transits from other planets. Refinement gives cleaner depth measurements.

**Implementation:**
- `refine_transit_depth(time, flux, best, others)` in `transit.py`
- Masks in-transit points from other **clean** (unflagged) candidates
- Phase-folds at best candidate's period, re-measures depth from medians
- If refined depth <= 0, demotes candidate and tries next one
- Flagged candidates are NOT masked (likely spurious; masking short-period
  false positives can corrupt >30% of data)

**Output:** `transit_depth_raw_ppm` (original BLS), `transit_depth_ppm` (refined).
Planet properties are recomputed from the refined depth.

## Feature 8: Candidate Re-Ranking with HZ Priority

**Purpose:** Promote habitable zone planets over non-HZ candidates when both pass
validation, making the "best candidate" the most scientifically interesting one.

**Implementation:**
- `_rerank_candidates()` in `transit.py`
- Three-tier ranking: (1) clean + HZ, (2) clean + non-HZ, (3) flagged
- Within tiers: sorted by descending SDE
- If re-ranking changes the #1 candidate, top-level results are updated

## Feature 9: Candidate Sanity Checks

**Purpose:** Flag physically implausible BLS detections before they contaminate
the final results.

**Transit-level flags** (`_candidate_to_result()` in `transit.py`):
- `duration_exceeds_15pct_of_period`: implausible transit geometry
- `fewer_than_3_transits`: insufficient data for reliable detection
- `very_deep_transit`: depth > 5%, likely eclipsing binary
- `negative_depth`: non-physical BLS artifact

**Planet-level flags** (`compute_planet_properties()` in `transit.py`):
- `orbit_inside_star`: orbital distance less than stellar radius
- `planet_larger_than_half_star`: likely eclipsing binary
- `extreme_temperature`: T_eq > 4000 K, planet would be destroyed

Clean candidates (`transit_flag="ok"`) are preferred over flagged ones, even at
lower SDE.

## Feature 10: Period Alias Candidate Generation

**Purpose:** BLS commonly detects at harmonic aliases (2x or 0.5x the true period),
causing ~38% planet radius error. Generate alias candidates at P/2 and 2P for each
BLS detection so that downstream ranking selects the true period.

**Implementation:**
- `_measure_phase_fold_depth(time, flux, period, epoch, duration)` in `transit.py`
- `_generate_alias_candidates(time, flux, candidates, baseline, min_period)` in `transit.py`
- Called in `detect_transit()` after stratified extraction, before candidate conversion
- Alias candidates tagged with `alias_of` field, flow through full pipeline

**Behavior:**
- P/2 alias generated only if >= min_period (default 0.5 days)
- 2P alias generated only if <= baseline
- Only aliases with positive measured depth and >= 10 in-transit points are kept

## Files Modified

| File | Changes |
|------|---------|
| `src/transit.py` | +7 functions, modify `analyze_transit` signature |
| `src/pipeline.py` | +1 function, modify `process_star` and `_run_transit_analysis` |
| `run_stars.py` | +4 CLI args, update `run()` and `format_result()` |
| `tests/test_transit.py` | +3 test classes (17 new tests) |
| `tests/test_pipeline.py` | +1 test class (4 new tests) |
| `docs/detailed_pipeline_gaia_based.md` | Add Module 5c section, expand Module 5a |

## Testing Strategy

1. Unit tests for each new function with synthetic data
2. Existing 96 tests must continue passing (no regressions)
3. Integration test: `run_stars.py --name "KIC 6922244" --transit --hz-only`
