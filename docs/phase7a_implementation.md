# Phase 7a Implementation Plan

**Date: 2026-02-19**
**Prerequisite: Phase 7 PoC plan approved (docs/phase7_exo_earth_poc.md)**

---

## Scope

Phase 7a builds the foundation modules for the multi-method exo-Earth targeted search:

1. **Target catalog** (`src/targets.py`) -- hardcoded 10-target list with Gaia DR3 validation
2. **RV data retrieval** (`src/rv_data.py`) -- query NASA Exoplanet Archive + DACE for RV data
3. **Proper motion anomaly** (`src/proper_motion.py`) -- Gaia-Hipparcos PMa analysis
4. **Per-target report** (`src/target_report.py`) -- orchestrator producing structured reports

---

## 1. Target Catalog Module (`src/targets.py`)

### Design

- Hardcoded `TARGET_CATALOG` list of 10 targets from Phase 7 plan Section 3.5
- Each entry: common name, HD number, TIC ID (for TESS), Gaia DR3 source_id, HIP ID, spectral type, Teff, distance, mass, known planet count
- `get_targets()` -- returns the full catalog
- `get_target(name)` -- lookup by common name or HD number
- `validate_against_gaia(target)` -- queries Gaia DR3 via existing `data_access.query_stars_by_id()`, compares catalog values with Gaia DR3 values, returns discrepancy report
- `check_data_availability(target)` -- checks TESS light curve availability via `lightcurve.search_lightcurve()` and reports sector count

### Dependencies
- `src/data_access.py` (existing)
- `src/lightcurve.py` (existing)

### Target IDs needed
Stars require: SIMBAD name, HD number, TIC ID, Gaia DR3 source_id, HIP ID. These will be resolved during initial validation using `data_access.resolve_simbad_name()`.

---

## 2. RV Data Retrieval Module (`src/rv_data.py`)

### Design

Two data sources, both publicly accessible without authentication:

**A. NASA Exoplanet Archive (known planets)**
- TAP query via astroquery for confirmed planets around our targets
- Returns: planet name, period, mass/radius, discovery method, reference
- Endpoint: `https://exoplanetarchive.ipac.caltech.edu/TAP/sync`

**B. DACE REST API (RV time series)**
- Geneva Observatory public RV database
- REST API: `https://dace-api.unige.ch/`
- Endpoints: search by target name, retrieve RV measurements
- Returns: BJD, RV (m/s), RV_err (m/s), instrument
- No authentication for published data

### Functions
- `query_known_planets(star_name)` -- query NASA Exoplanet Archive for confirmed/candidate planets
- `query_dace_rv(star_name)` -- retrieve public RV time series from DACE
- `rv_periodogram(time, rv, rv_err)` -- Lomb-Scargle on RV data (reuses periodogram.py logic)
- `rv_detection_limit(time, rv_err, period_grid)` -- analytical K_min at each period (semi-amplitude detection floor)

### Dependencies
- `requests` (stdlib-adjacent, should be installed)
- `astroquery` (existing)
- `numpy`, `scipy` (existing)
- Reuses Lomb-Scargle from `astropy.timeseries` (already used in periodogram.py)

### No new pip packages needed

---

## 3. Proper Motion Anomaly Module (`src/proper_motion.py`)

### Design

The proper motion anomaly (PMa) is the difference between Gaia DR3 proper motion (short-term, quasi-instantaneous) and Hipparcos proper motion (measured circa 1991). A significant difference indicates an unseen companion pulling the star over the ~25-year baseline.

### Functions
- `query_gaia_hipparcos(source_id_or_hip)` -- retrieve PM from Gaia DR3 and Hipparcos-2 catalog via astroquery
- `compute_pma(pm_gaia, pm_hip, pm_gaia_err, pm_hip_err)` -- compute PMa magnitude and significance (SNR)
- `pma_companion_sensitivity(distance_pc, stellar_mass, pma_snr_threshold)` -- estimate minimum detectable companion mass as function of separation, given PMa precision

### Data Sources
- Gaia DR3: `pmra`, `pmdec`, `pmra_error`, `pmdec_error` from `gaiadr3.gaia_source`
- Hipparcos-2: via `astroquery.vizier` querying catalog `I/311/hip2` (van Leeuwen 2007)

### Dependencies
- `astroquery` (existing -- Gaia + Vizier modules)
- `numpy` (existing)

---

## 4. Per-Target Report Generator (`src/target_report.py`)

### Design

Orchestrates all available analyses for a single target star:

1. Stellar properties via existing pipeline (Modules 1-3)
2. Light curve + variability analysis (Module 4)
3. Transit search (Module 5)
4. Known planets from NASA Exoplanet Archive
5. RV data inventory (DACE query + periodogram)
6. Proper motion anomaly
7. Combined data availability summary

### Functions
- `generate_report(target_name, include_rv=True, include_pma=True, include_transit=True)` -- runs all analyses, returns structured dict
- `format_report_markdown(report_dict)` -- formats report dict as markdown string for human review
- `save_report(report_dict, output_dir)` -- saves JSON + markdown to output directory

### Output Structure
```
{
  "target": { name, hd, tic, gaia_id, ... },
  "stellar_properties": { distance, teff, radius, mass, ... },
  "lightcurve": { available, n_sectors, variability_class, ... },
  "transit": { detected, candidates, ... },
  "known_planets": [ { name, period, mass, method, ... }, ... ],
  "rv_data": { n_measurements, instruments, time_baseline, periodogram_peaks, ... },
  "proper_motion_anomaly": { pma_mas_yr, snr, significant, ... },
  "data_summary": { methods_available, sensitivity_notes }
}
```

---

## 5. Testing Plan

### Unit Tests (`tests/test_phase7a.py`)

| Test | What it verifies |
|------|-----------------|
| `test_target_catalog_completeness` | All 10 targets present with required fields |
| `test_target_lookup` | Lookup by common name and HD number works |
| `test_gaia_validation_alpha_cen` | Gaia DR3 query returns valid data for alpha Cen A |
| `test_known_planets_query` | NASA Exoplanet Archive returns known planets for 82 G. Eridani |
| `test_pma_computation` | PMa computation with synthetic data gives expected result |
| `test_rv_detection_limit` | Detection limit calculation is physically reasonable |
| `test_report_structure` | Report dict has all required top-level keys |

### Integration Test
- Run full report on **82 G. Eridani** (HD 20794) as validation
- Should find 3 confirmed planets, HARPS/ESPRESSO RV data, TESS sectors 3/4/30/31
- This is Phase 7b's first deep dive target

---

## 6. Architecture Update

New modules added to `src/`:
```
src/
  __init__.py         (existing)
  data_access.py      (existing, no changes)
  distance.py         (existing, no changes)
  temperature.py      (existing, no changes)
  mass.py             (existing, no changes)
  lightcurve.py       (existing, no changes)
  periodogram.py      (existing, no changes)
  transit.py          (existing, no changes)
  pipeline.py         (existing, no changes)
  targets.py          (NEW - Phase 7 target catalog)
  rv_data.py          (NEW - RV data retrieval + analysis)
  proper_motion.py    (NEW - Gaia-Hipparcos PMa)
  target_report.py    (NEW - per-target report orchestrator)
```

No changes to existing modules. All new code follows existing patterns:
- Logging via `logging.getLogger(__name__)`
- No print statements
- Functions return dicts
- Errors handled gracefully with None returns and log warnings
