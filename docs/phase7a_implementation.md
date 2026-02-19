# Phase 7a Implementation -- COMPLETE

**Date: 2026-02-19**
**Status: Complete** -- all modules implemented, 27 unit tests + 5 integration checks passing
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

### SIMBAD Name Resolution

Stars require: SIMBAD name, HD number, TIC ID, Gaia DR3 source_id, HIP ID. Resolution uses `data_access.resolve_simbad_name()`.

**Gotcha**: Some common names (e.g., "82 G. Eridani") fail SIMBAD lookup. The module implements HD-number fallback: `resolve_target_ids()` tries the common name first, then falls back to the HD number (e.g., "HD 20794"). This pattern is also used in `target_report.py` for stellar property and PMa queries.

---

## 2. RV Data Retrieval Module (`src/rv_data.py`)

### Design

Two data sources, both publicly accessible without authentication:

**A. NASA Exoplanet Archive (known planets + orbital solutions)**
- TAP query via astroquery for confirmed planets around our targets
- Queries the `ps` (Planetary Systems) table
- Returns: planet name, period, mass/radius, discovery method, reference
- Note: No raw RV time series table exists on the archive; only orbital solutions

**B. DACE (RV time series via `dace-query` v2.0.0)**
- Geneva Observatory public RV database
- Uses official `dace-query` Python package (not the REST API, which has DNS issues)
- `Spectroscopy.get_timeseries(target, sorted_by_instrument=False, output_format='dict')` returns flat dict with parallel arrays
- Returns: RJD (BJD-2400000), RV (km/s), RV_err (km/s), instrument, DRS QC flags
- Name resolution: generates HD-format names without spaces (e.g., "HD20794") via `_generate_dace_names()`
- Data filtering: removes NaN/infinite values and failed DRS QC; converts km/s to m/s

### Functions
- `query_known_planets(star_name)` -- query NASA Exoplanet Archive for confirmed/candidate planets
- `query_dace_rv(star_name)` -- retrieve public RV time series from DACE via dace-query
- `query_nasa_rv_data(star_name)` -- fallback: orbital solutions from `ps` table
- `rv_periodogram(time, rv, rv_err)` -- Lomb-Scargle on RV data
- `rv_detection_limit(time, rv_err, period_grid)` -- analytical K_min at each period (semi-amplitude detection floor)
- `rv_to_planet_mass(k_ms, period_days, stellar_mass_msun)` -- convert RV semi-amplitude to minimum planet mass

### Dependencies
- `dace-query` v2.0.0 (new; added to requirements.txt)
- `astroquery` (existing)
- `numpy`, `scipy` (existing)
- `astropy.timeseries.LombScargle` (existing)

---

## 3. Proper Motion Anomaly Module (`src/proper_motion.py`)

### Design

The proper motion anomaly (PMa) is the difference between Gaia DR3 proper motion (epoch ~2016, quasi-instantaneous) and Hipparcos proper motion (epoch ~1991.25). A significant difference indicates an unseen companion pulling the star over the ~24.75-year baseline.

### Functions
- `query_gaia_proper_motion(source_id)` -- retrieve PM from Gaia DR3 via TAP+
- `query_hipparcos_proper_motion(hip_id)` -- retrieve PM from Hipparcos-2 catalog (`I/311/hip2`) via VizieR
- `compute_pma(gaia_pm, hip_pm)` -- compute vector PMa with error propagation and SNR
- `pma_companion_mass(pma_mas_yr, distance_pc, stellar_mass_msun, separation_au)` -- estimate companion mass from PMa signal
- `analyze_pma(gaia_source_id, hip_id, ...)` -- full pipeline: query both catalogs, compute PMa, assess significance
- `assess_ruwe(ruwe_value)` -- interpret RUWE (>1.4 indicates possible unresolved companion)

### Companion Mass Formula

The corrected formula (initial version had M_c inversely proportional to separation, which is physically wrong):

```
M_c = (PMa_rad * d_AU * a^2) / (G_AU * Delta_t)
```

Where G_AU = 4*pi^2 in AU^3/(Msun*yr^2), Delta_t = 24.75 yr (Gaia-Hipparcos epoch difference). M_c is proportional to separation^2 because a more distant companion needs more mass to produce the same angular acceleration.

### Data Sources
- Gaia DR3: `pmra`, `pmdec`, `pmra_error`, `pmdec_error`, `ruwe` from `gaiadr3.gaia_source`
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

### Internal Runners
- `_run_stellar_properties(target)` -- Gaia DR3 via pipeline (with HD fallback for SIMBAD)
- `_run_known_planets(target)` -- NASA Exoplanet Archive TAP
- `_run_lightcurve(target)` -- TESS sector count and variability via lightkurve
- `_run_transit(target)` -- deferred (placeholder; full implementation in Phase 7b)
- `_run_rv_analysis(target)` -- DACE time series + Lomb-Scargle periodogram
- `_run_pma(target)` -- Gaia-Hipparcos PMa with RUWE

### Output Structure
```
{
  "target": { name, hd, tic, gaia_id, ... },
  "stellar_properties": { distance_pc, teff_K, radius_Rsun, mass_Msun, luminosity_Lsun, pipeline_value, catalog_value },
  "lightcurve": { available, n_sectors, variability, amplitude_ppt },
  "transit": { status: "deferred" },
  "known_planets": [ { name, period_days, mass_earth, radius_earth, method, year }, ... ],
  "rv_data": { source, n_measurements, baseline_days, instruments, best_period_days, best_period_fap },
  "proper_motion_anomaly": { pma_total_mas_yr, snr, significant, ruwe },
  "data_summary": { methods_available, n_confirmed_planets, notes }
}
```

---

## 5. Testing -- Results

### Unit Tests (`tests/test_phase7a.py`) -- 27/27 passing

| Test Class | Count | What it verifies |
|------------|-------|-----------------|
| `TestTargetCatalog` | 11 | Catalog completeness, field validation, lookup by name/HD, distance/Teff/mass ranges, tier labels, required IDs |
| `TestRVAnalysis` | 6 | Synthetic periodogram recovery, detection limits, mass conversion, period grid generation |
| `TestProperMotionAnomaly` | 7 | Zero PMa baseline, significant detection, companion mass estimation (M_c proportional to a^2), RUWE assessment, error propagation |
| `TestReportStructure` | 3 | Report dict keys, markdown formatting, JSON date serializer |

### Integration Test (`tests/test_82_eridani_integration.py`) -- 5/5 checks passing

- **Gaia ID**: Resolved to source_id 4847957293278177024 (via HD 20794 fallback)
- **Known planets**: 4 confirmed (b, d, e, f) from NASA Exoplanet Archive
- **PMa**: 3.34 mas/yr, SNR=15.5, significant=True, RUWE=1.98
- **Report**: Full markdown + JSON generated at `output/target_reports/`
- **Stellar properties**: Pipeline Teff=5381 K (catalog: 5401 K, 0.4% error)

### DACE RV Validation (82 G. Eridani / HD 20794)

- 12,303 measurements after QC filtering (from 12,356 raw)
- 5 instruments: CORAVEL-S, ESPRESSO18, ESPRESSO19, HARPS03, HARPS15
- Baseline: 42.2 years (15,400 days)
- Best Lomb-Scargle period: 414.0 days (FAP=0.0) -- does not match any known planet period; likely stellar activity or window function alias (to investigate in Phase 7b)

---

## 6. Bugs Found and Fixed

| Bug | Root Cause | Fix |
|-----|-----------|-----|
| **PMa companion mass formula wrong** | Had M_c inversely proportional to separation; physically, closer companions need *less* mass to produce same PMa | Corrected to M_c = PMa_rad * d_AU * a^2 / (G_AU * Delta_t). Caught by unit test (at 2 AU mass was 64.8 Mjup vs 25.9 at 5 AU). |
| **SIMBAD can't resolve "82 G. Eridani"** | `Simbad.query_object("82 G. Eridani")` returns NoResultsWarning; name format not in SIMBAD's primary aliases | Added HD number fallback in `resolve_target_ids()`, `_run_stellar_properties()`, and `_run_pma()` |
| **DACE REST API DNS failure** | `dace-api.unige.ch` does not resolve from this machine | Replaced REST API with `dace-query` v2.0.0 pip package (uses different backend `obs-webapp.obsuksprd2.unige.ch`) |
| **NASA `radialvelocity` TAP table doesn't exist** | The table name was assumed; NASA archive has no raw RV time series table | Changed to query orbital solutions from the `ps` (Planetary Systems) table instead |

---

## 7. Architecture Update

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

New dependencies in `requirements.txt`:
- `dace-query` (v2.0.0) -- DACE RV time series access

No changes to existing modules. All new code follows existing patterns:
- Logging via `logging.getLogger(__name__)`
- No print statements
- Functions return dicts
- Errors handled gracefully with None returns and log warnings

---

## 8. Open Questions for Phase 7b

1. **414-day RV period**: The dominant Lomb-Scargle peak at 414 days does not correspond to any of the 4 known planet periods (18.3, 89.7, 147.0, 647.6 d). Is this a stellar activity cycle, a window function alias, or an uncatalogued signal?
2. **Elevated RUWE (1.98) + significant PMa (SNR=15.5)**: Both indicators suggest an outer companion beyond the 4 confirmed planets. What mass/separation constraints can we derive?
3. **TESS transit search**: 26 sectors available but not yet analyzed. Can we constrain planet radii for the RV-detected planets?
4. **RV detection sensitivity**: With 12,303 measurements and 42-year baseline, what is the minimum detectable planet mass as a function of orbital period?
