# Phase 2: Light Curve Analysis via lightkurve/MAST -- Plan

## Goal

Add a new Module 4 (Light Curve Analysis) to the pipeline that downloads space-based time-series photometry (TESS, Kepler, K2) from MAST, detects periodic signals, and classifies stellar variability. This connects directly to the Phase 1 pipeline and enables the project's long-term vision: "specify a coordinate in the sky, resolve to object, fetch light curve data, and generate information about the object."

## What lightkurve and MAST Are

**MAST** (Mikulski Archive for Space Telescopes) is NASA's archive for Hubble, TESS, Kepler, K2, JWST, and other missions. It holds the actual time-series photometry (flux vs time) for millions of stars observed by space telescopes.

**lightkurve** is an Astropy-affiliated Python package maintained by STScI (Space Telescope Science Institute) that wraps MAST access. It handles:
- Searching for available observations by star name, coordinates, or catalog ID
- Downloading light curve (pre-extracted flux) or target pixel files (raw CCD cutouts)
- Stitching multi-sector/quarter data into continuous time series
- Detrending and flattening (removing spacecraft systematics)
- Periodogram analysis (Lomb-Scargle, Box Least Squares for transits)
- Folding, binning, and visualization

Key point: lightkurve uses MAST under the hood, so you get both in one package.

### Coverage

| Mission | Sky coverage | Cadence | Duration | Catalog ID |
|---------|-------------|---------|----------|------------|
| TESS | ~85% of sky (ongoing) | 2 min, 10 min, 30 min | 27 days/sector, many re-observed | TIC (TESS Input Catalog) |
| Kepler | ~116 sq deg (Cygnus-Lyra) | 1 min, 30 min | 4 years continuous | KIC (Kepler Input Catalog) |
| K2 | ~20 fields along ecliptic | 1 min, 30 min | ~80 days/campaign | EPIC |

Most bright stars (G < 13) have TESS data. Kepler field stars have the best precision and longest baseline.

---

## How It Connects to Phase 1

```
Phase 1 (existing):
  Gaia catalog data -> [Module 1: Distance] -> [Module 2: Teff/L/R] -> [Module 3: Mass]

Phase 2 (new):
  Star identifier
      |
      +---> [Module 4a: Light Curve Retrieval] -- search MAST, download, stitch, clean
      |
      +---> [Module 4b: Period Detection] -- Lomb-Scargle periodogram
      |
      +---> [Module 4c: Variability Classification] -- periodic/transit/stochastic/quiet
      |
      +---> Feed period back into Module 1 (Cepheid period -> distance)
```

### Integration points

1. **Identifier bridging**: lightkurve can search by common name ("Proxima Cen"), coordinates (RA/Dec), or catalog ID (TIC, KIC). We already resolve names via SIMBAD and have Gaia coordinates -- both work as lightkurve search inputs.

2. **Cepheid loop**: If Module 4 detects a period and the star is classified as a pulsator, feed `cepheid_period_days` back into Module 1 instead of requiring manual input.

3. **Activity/rotation**: Detected rotation periods inform stellar age (gyrochronology) -- a future Module 5 possibility.

---

## Module 4 Design

### 4a: Light Curve Retrieval (`src/lightcurve.py`)

**Input**: Star identifier (name string, or dict with RA/Dec coordinates)

**Steps**:
1. Search MAST for available light curves: `lk.search_lightcurve(target, author=...)`.
   - Try TESS first (broadest sky coverage), fall back to Kepler/K2.
   - Search priority: `['SPOC', 'Kepler', 'K2']` (SPOC = TESS Science Processing Operations Center).
2. Download light curve collection. For multiple sectors/quarters, download all available.
   - Configurable limit: `max_sectors` parameter (default: all available, set a sensible limit like 20 for testing).
3. Stitch into a single continuous light curve: `lc_collection.stitch()`.
4. Clean: remove NaN flux values, remove outliers (>5-sigma clips).
5. Flatten: remove long-term trends with Savitzky-Golay filter: `lc.flatten(window_length=...)`.
   - Window length depends on expected signal timescale. Default: 1001 cadences.
   - Preserve short-period signals while removing spacecraft systematics.

**Output**:
```python
{
    "mission": "TESS",           # or "Kepler", "K2"
    "author": "SPOC",
    "n_sectors": 5,              # number of sectors/quarters used
    "n_points": 18432,           # total data points after cleaning
    "time_baseline_days": 135.2, # total observation span
    "cadence_s": 120,            # exposure cadence in seconds
    "time": [...],               # array of times (BJD - offset)
    "flux": [...],               # normalized, flattened flux array
}
```

**Key lightkurve calls** (from user's working code, extended):
```python
import lightkurve as lk

# Search
search_res = lk.search_lightcurve(target, author='SPOC')  # TESS
# or: lk.search_lightcurve(target, author='Kepler')

# Download all available sectors
lc_collection = search_res.download_all()

# Stitch into single light curve
lc = lc_collection.stitch()

# Clean and flatten
lc = lc.remove_nans().remove_outliers(sigma=5)
lc_flat = lc.flatten(window_length=1001)
```

### 4b: Period Detection

**Method**: Lomb-Scargle Periodogram (lightkurve wraps astropy's implementation).

**Steps**:
1. Compute periodogram on flattened light curve.
   - Frequency range: 0.01 to 50 cycles/day (periods from 100 days down to ~30 minutes).
   - Oversample factor: 5 (standard for reliable peak detection).
2. Identify the dominant peak (highest power).
3. Compute False Alarm Probability (FAP) for the peak.
   - FAP < 0.001 (0.1%): confident detection.
   - FAP < 0.01 (1%): likely detection.
   - FAP >= 0.01: no significant period.
4. Optionally run Box Least Squares (BLS) periodogram if transit detection is desired.
   - BLS is optimized for box-shaped dips (planetary transits, eclipsing binaries).

**Output**:
```python
{
    "period_days": 3.5225,
    "period_uncertainty_days": 0.002,  # from peak width
    "period_fap": 1.2e-15,            # false alarm probability
    "period_power": 0.87,             # normalized LS power at peak
    "amplitude_ppt": 2.5,             # peak-to-peak amplitude in parts per thousand
    "n_harmonics_detected": 2,        # secondary peaks at P/2, P/3, etc.
    "detection_method": "lomb_scargle",
}
```

**Key lightkurve calls**:
```python
# Lomb-Scargle
pg = lc_flat.to_periodogram(method='lombscargle',
                             minimum_period=0.02,  # days
                             maximum_period=100)    # days

best_period = pg.period_at_max_power
best_power = pg.max_power
fap = pg.false_alarm_probability(best_power)

# Phase fold at detected period
lc_folded = lc_flat.fold(period=best_period)

# Box Least Squares (for transits)
pg_bls = lc_flat.to_periodogram(method='bls',
                                 minimum_period=0.5,
                                 maximum_period=20)
```

### 4c: Variability Classification

Simple rule-based classification (not ML -- keep it straightforward for v2.0):

| Classification | Criteria |
|---------------|----------|
| `periodic_pulsator` | Significant LS period detected (FAP < 0.001), smooth sinusoidal phase curve, amplitude > 1 ppt |
| `eclipsing_binary` | Significant LS period, phase curve shows sharp dip(s), BLS detection |
| `transit_candidate` | BLS detection with box-shaped dip, depth < 2%, LS period weak or absent |
| `rotational_variable` | Significant LS period, quasi-sinusoidal, amplitude < 10 ppt, period 1-50 days |
| `stochastic_variable` | No significant period, but RMS scatter > 3x median photometric error |
| `quiet` | No significant period, scatter consistent with photometric noise |

**Implementation**: Start with just the basics -- `periodic` (FAP < 0.001), `possible_periodic` (FAP < 0.01), `non_variable` (FAP >= 0.01). Expand to the full classification table later.

**Output**:
```python
{
    "variability_class": "periodic",
    "variability_flag": "ok",        # or "low_snr", "few_points", "short_baseline"
}
```

---

## Implementation Plan

### Step 1: Setup and Dependencies

- Add `lightkurve` to requirements.txt.
- lightkurve depends on: astropy, matplotlib, requests, beautifulsoup4, astroquery (already have this), scipy (already have this), numpy (already have this).
- Verify install: `.\venv\Scripts\pip install lightkurve`.
- Security: lightkurve is maintained by STScI (NASA), widely used in the astronomy community. No known vulnerabilities or policy concerns.

### Step 2: Implement `src/lightcurve.py` -- Light Curve Retrieval

Create `src/lightcurve.py` with:
- `search_lightcurve(target, mission_priority=None, max_sectors=20)` -- search MAST.
- `download_and_stitch(search_result, max_sectors=20)` -- download, stitch, clean.
- `flatten_lightcurve(lc, window_length=1001)` -- flatten with SG filter.
- All functions return plain Python dicts/arrays (not lightkurve objects) so the rest of the pipeline stays decoupled.
- Logging, not print statements.

### Step 3: Implement `src/periodogram.py` -- Period Detection

Create `src/periodogram.py` with:
- `detect_period(time, flux, min_period=0.02, max_period=100)` -- Lomb-Scargle.
- `compute_fap(periodogram, max_power)` -- false alarm probability.
- `classify_variability(period_result, lc_stats)` -- simple rule-based classifier.
- Returns dict matching the output schema above.

### Step 4: Integrate into Pipeline

- Update `src/pipeline.py`: add optional Module 4 call after Module 3.
  - Module 4 is optional because it requires network access to MAST and not all stars have light curves.
  - Add `include_lightcurve=False` parameter to `process_star()`.
  - When enabled and a period is detected on a Cepheid candidate, update Module 1 distance.
- Update `run_stars.py`: add `--lightcurve` flag to enable Module 4.
- Update output schema with Module 4 fields.

### Step 5: Unit Tests (`tests/test_lightcurve.py`, `tests/test_periodogram.py`)

**test_lightcurve.py**:
- Test search function returns results for known Kepler target (KIC 6922244).
- Test download and stitch produces valid time/flux arrays.
- Test NaN removal and outlier clipping.
- Test flatten preserves periodic signal.
- Note: these tests hit the network (MAST). Mark with `@pytest.mark.network` so they can be skipped in CI.

**test_periodogram.py**:
- Test period detection on synthetic sinusoidal signal (known period, no network needed).
- Test FAP computation returns low value for strong signal.
- Test FAP returns high value for pure noise.
- Test classification rules produce correct labels.
- Test with user's known example: KIC 6922244 at period 3.5225 days.

### Step 6: Validation and Demo Script

- Run on user's working example: KIC 6922244 (period 3.5225 days) and verify detection.
- Run on Delta Cephei: verify pipeline auto-detects pulsation period and improves distance estimate.
- Run on a quiet star (e.g., Sun-like with no significant variability) to test non-detection path.
- Update `run_stars.py` predefined stars with light curve results where available.

### Step 7: Documentation

- Update `docs/detailed_pipeline_gaia_based.md` with Module 4 spec.
- Update README with lightkurve usage examples.

---

## Validation Targets for Phase 2

| Star | Mission | Expected | Purpose |
|------|---------|----------|---------|
| KIC 6922244 | Kepler | P ~ 3.5225 d, eclipsing/transit | User's known working example |
| Delta Cephei | TESS | P ~ 5.37 d, pulsator | Cepheid period -> distance feedback loop |
| Proxima Centauri | TESS | Flare star, rotation ~83 d | Tests stochastic/rotation detection |
| Sun (via solar analog) | Kepler (KIC 10644253) | Quiet, rotation ~10-25 d | Tests low-amplitude / non-detection |

---

## File Changes Summary

| File | Action | Description |
|------|--------|-------------|
| `src/lightcurve.py` | NEW | Light curve search, download, stitch, clean, flatten |
| `src/periodogram.py` | NEW | Lomb-Scargle period detection, FAP, variability classification |
| `src/pipeline.py` | EDIT | Add optional Module 4 integration |
| `src/data_access.py` | EDIT (minor) | Potentially add coordinate extraction helper for lightkurve search |
| `run_stars.py` | EDIT | Add `--lightcurve` flag |
| `tests/test_lightcurve.py` | NEW | Light curve retrieval tests |
| `tests/test_periodogram.py` | NEW | Period detection and classification tests |
| `requirements.txt` | EDIT | Add `lightkurve` |
| `docs/detailed_pipeline_gaia_based.md` | EDIT | Add Module 4 specification |

---

## Risk Assessment and Scope Control

**Network dependency**: MAST queries can be slow (5-30s per download). Module 4 should be opt-in, not default.

**Data volume**: A single Kepler star (17 quarters) can be ~70,000 data points. TESS multi-sector can be similar. This is manageable in memory but downloading many stars sequentially would be slow. For Phase 2, process one star at a time -- batch optimization is a Phase 3 concern.

**Not in scope for Phase 2**:
- Target pixel file (TPF) analysis or custom aperture photometry (use pre-extracted light curves only).
- Transit fitting (depth, duration, limb darkening) -- just detection.
- Gyrochronology / age estimation from rotation periods.
- ML-based classification.
- Batch processing of many stars.

**Fallback**: If a star has no MAST data, Module 4 gracefully returns `{"lightcurve_available": false}` and the rest of the pipeline is unaffected.

---

## References

- Lightkurve Collaboration (2018). Lightkurve: Kepler and TESS time series analysis in Python. Astrophysics Source Code Library. ascl:1812.013
- Ricker, G. R. et al. (2015). Transiting Exoplanet Survey Satellite (TESS). JATIS 1, 014003.
- VanderPlas, J. T. (2018). Understanding the Lomb-Scargle Periodogram. ApJS 236, 16.
- Kovacs, G., Zucker, S. & Mazeh, T. (2002). A box-fitting algorithm in the search for periodic transits. A&A 391, 369.
