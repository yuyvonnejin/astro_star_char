# Architecture and Data Sources -- Master Reference

Last updated: 2026-07-03 (Phase 8.4)

Purpose: one place that answers "what does each module in src/ do,
where does its input data come from (online API, static table,
generated file), and what do the tests actually exercise (real
network, mocked, or synthetic)".

Companion docs: docs/project_status.md (where we are, how to run),
docs/phase8_survey_plan.md (survey design), phase plan docs per phase.

---

## 1. Data provenance categories

Every input in this project falls into one of five categories:

| Category | What it is | Examples |
|---|---|---|
| ONLINE | Live query to an external archive at runtime | Gaia DR3 TAP, SIMBAD, MAST, DACE, NASA Exoplanet Archive, Vizier (Hipparcos-2) |
| STATIC (code) | Hardcoded tables/constants in src/, from published literature | Mucciarelli Teff coefficients, curated 10-star TARGET_CATALOG, spectral-type anchors, name aliases |
| GENERATED (data/) | JSON files produced once by a script, then read at runtime | data/shell_catalog.json (from SIMBAD TAP) |
| CURATED (data/) | JSON files hand-written from literature reading | data/known_planets.json |
| SYNTHETIC (tests) | Signals injected in test code; no external data | Sinusoid RV series, box transits, canned star dicts |

Key design rule (Phase 8): curated data beats online data. The NASA
archive carries stale/disputed planet rows, so data/known_planets.json
is the authority for which planets get fitted; NASA is only consulted
for targets with no curated entry.

---

## 2. Online API inventory

| API | Accessed via | Calling functions | Used for | Reliability notes |
|---|---|---|---|---|
| Gaia DR3 TAP | astroquery.gaia | data_access.query_stars_by_id, query_cepheids, resolve_simbad_name (cone search); proper_motion.query_gaia_proper_motion | Parallax, photometry, FLAME reference values, PM/RUWE | Works; cosmetic password warning |
| SIMBAD | astroquery.simbad | data_access.resolve_simbad_name (query_objectids); targets.resolve_target_ids (query_objectids); targets.build_shell_catalog (query_tap) | Name -> Gaia/TIC/HIP IDs; shell catalog generation | ADQL quirk: ORDER BY needs column alias |
| MAST | lightkurve | lightcurve.search_lightcurve, download_and_stitch | TESS/Kepler/K2 light curves | Filter to SPOC author; stitch() needs time sort |
| DACE | dace-query >= 3.0 | rv_data.query_dace_rv | RV time series + activity indicators (fwhm, bis, smw, rhk, halpha) | 2.x endpoints dead; datasets change server-side; catastrophic outliers normal |
| NASA Exoplanet Archive TAP | requests.get | rv_data.query_known_planets, query_nasa_rv_data | Known-planet lookup (uncurated targets only) | Intermittent timeouts; stale rows; hardcoded fallback in deep_dive |
| Vizier (Hipparcos-2, I/311/hip2) | astroquery.vizier | proper_motion.query_hipparcos_proper_motion | Historic proper motion for PMa | Stable |

Everything else in src/ is offline: pure computation on arrays passed
in, plus reads of the two data/ JSON files.

---

## 3. src/ modules

### 3.1 Pure computation (no network, no file I/O)

Deterministic given inputs; fully testable offline.

| Module | Purpose | Static data embedded (literature source) |
|---|---|---|
| distance.py | Bayesian parallax inversion; Cepheid Leavitt law | Parallax zero-point 0.017 mas, prior scale 1350 pc (Bailer-Jones 2015) |
| temperature.py | Teff from BP-RP color; BC_G; luminosity; radius | DWARF_COEFFS/GIANT_COEFFS (Mucciarelli+2021), BC_COEFFS_HOT/COOL (Andrae+2018), M_bol_sun, T_eff_sun |
| mass.py | Piecewise mass-luminosity relation | Power-law exponents per luminosity range; MS cuts (logg >= 3.5, 3300-8000 K) |
| periodogram.py | Lomb-Scargle period detection + variability class | FAP thresholds 0.001/0.01 |
| transit.py | BLS transit search + planet characterization | HZ coefficients (Kopparapu+2013), size classes (Fulton+2017), BLS grid settings |
| sensitivity.py | Analytic detection limits: transit + RV + astrometry | SNR thresholds (7.1 transit, 5.0 RV), physical constants |
| rv_keplerian.py | RadVel wrapper: Keplerian fit + GP + MCMC | DEFAULT_GP_HYPERPARAMS bounds (gp_per 500-10000 d etc.); MCMC defaults (nwalkers=50, nsteps=1000, ensembles=3). Gotcha: _sync_vary_flags() required after any vary change (RadVel >= 1.4 caches flags) |

### 3.2 Online-access modules

| Module | Purpose | Online sources | Static data embedded |
|---|---|---|---|
| data_access.py | Gaia DR3 + SIMBAD star property queries | Gaia TAP, SIMBAD | ADQL templates; parallax_over_error > 10 cut; cone radius 0.005 deg. Fallback: SIMBAD ID match -> Gaia cone search |
| lightcurve.py | MAST light curve retrieval + clean/bin/flatten | MAST via lightkurve | DEFAULT_AUTHOR_PRIORITY = SPOC, Kepler, K2 |
| rv_data.py | RV retrieval + periodogram + filtering + binning + decorrelation | DACE (dace-query), NASA TAP | 22 hardcoded name aliases (_generate_search_names); filter defaults (5-MAD clip, max_scatter_factor 50); nightly-bin convention (floor of JD). Fallback: DACE raw RVs -> NASA orbital solutions |
| proper_motion.py | Gaia-Hipparcos proper motion anomaly | Gaia TAP, Vizier I/311/hip2 | Epoch difference 24.75 yr; RUWE > 1.4 threshold; PMa SNR > 3 |
| targets.py | Target catalogs: curated + generated shell | SIMBAD (ID resolution, shell catalog TAP query) | Curated TARGET_CATALOG (10 stars, hand-picked, literature values); _SPTYPE_ANCHORS (Pecaut-Mamajek Teff/mass per spectral type); SHELLS = 0-6/6-10/10-15 pc. Reads/writes data/shell_catalog.json. Lookup order: curated catalog -> shell catalog -> (resolve_target_ids: SIMBAD live) |

### 3.3 Data-file modules

| Module | Purpose | Reads | Writes |
|---|---|---|---|
| reference_data.py | Loader for curated literature reference | data/known_planets.json (cached in-process, schema-checked) | none |

### 3.4 Orchestrators (chain other modules; write reports)

| Module | Purpose | Chains | Output files |
|---|---|---|---|
| pipeline.py | Modules 1-5 for one star (stellar props, optional LC + transit) | distance, temperature, mass, lightcurve, periodogram, transit | logs/pipeline.log; optional JSON via CLI |
| target_report.py | Phase 7a per-target report (stellar + LC + RV + PMa + known planets) | data_access, pipeline, lightcurve, rv_data, proper_motion | output/target_reports/*.json + .md |
| deep_dive.py | Phase 7b single-target deep analysis; hosts _run_rv_analysis (the validated RV chain) and _run_known_planets (NASA + hardcoded literature fallback for HD 20794) | targets, pipeline, lightcurve, transit, rv_data, rv_keplerian, proper_motion, sensitivity, reference_data | output/target_reports/*.json + .md |
| survey.py | Phase 8 multi-target survey driver: scorecard per target, candidate vetting, detection limits | targets, reference_data, deep_dive (lazy imports), rv_data | output/survey/survey_results.json + survey_summary.md |

survey.py static data: alias/veto constants (sidereal day, annual
365.25 d, veto tolerances), DEFAULT_SURVEY_TARGETS = HD 10700 /
HD 20794 / HD 115617.

### 3.5 Where the RV survey chain gets its data (end to end)

```
DACE (online, dace-query)          data/known_planets.json (curated)
        |                                   |
        v                                   v
rv_filter_instruments  ------>  which planets to fit + GP priors
        |                                   |
rv_bin_nightly                              |
        |                                   v
        +----------------->  fit_keplerian (RadVel, offline)
                                    |
                     rv_decorrelate_activity (DACE indicators)
                                    |
                             GP refit (offline)
                                    |
                     survey.build_scorecard (offline)
                                    |
              output/survey/*.json + *.md (written locally)
```

Uncurated targets (e.g. shell-catalog stars): same chain, but planet
list comes from NASA archive query (online), and if nothing is found
no fit is forced -- the scorecard reports data availability,
periodogram candidates (vetted), and detection limits only.

---

## 4. data/ files

| File | Origin | Regenerate/edit how | Consumed by |
|---|---|---|---|
| known_planets.json | CURATED by hand from literature (Nari+2025, Laliotis+2023, Feng+2017, Figueira+2025). Per-target: planets (period, K, confirmed flag, reference string) + pipeline config (rotation_days, activity_cycle_days, exclude_instruments, notes) | Edit by hand; adding a survey target = adding an entry here, no code change | reference_data.py -> deep_dive, survey |
| shell_catalog.json | GENERATED from SIMBAD TAP (2026-07-03): 176 FGK dwarfs < 15 pc, Teff 4000-6300 K. Teff/mass are APPROXIMATE (interpolated from spectral type via _SPTYPE_ANCHORS, not measured) | `python -c "from src.targets import build_shell_catalog; build_shell_catalog()"` (~30 s, network) | targets.targets_in_shell -> survey --shell |

Curated TARGET_CATALOG (in targets.py code, not data/) merges on top
of the shell catalog; curated wins on HD collision.

---

## 5. Output locations (all generated, safe to delete)

- output/survey/ -- survey_results.json + survey_summary.md
- output/target_reports/ -- per-target deep-dive JSON + markdown
- logs/pipeline.log -- runtime log

---

## 6. tests/ -- what is real, what is synthetic

Configuration: pytest.ini defines one marker, `network` (tests that
hit MAST/SIMBAD/NASA/DACE). No conftest.py. Invocation:

```
# Standard offline regression (136 tests, ~10 s)
./venv/Scripts/python.exe -m pytest tests/test_phase7a.py tests/test_phase7b.py tests/test_phase7b2.py tests/test_phase7c.py tests/test_phase8.py -q

# Everything except network
./venv/Scripts/python.exe -m pytest tests/ -m "not network" -q

# Network integration only (slow, needs internet)
./venv/Scripts/python.exe -m pytest tests/ -m network -q
```

261 test functions total across 13 files.

| File | Tests | Type | Data used |
|---|---|---|---|
| test_distance.py | 5 | OFFLINE | Canned star dicts (Proxima-like parallax, Cepheid period) |
| test_temperature.py | 9 | OFFLINE | Canned star dicts (solar BP-RP, giants, out-of-range colors) |
| test_mass.py | 8 | OFFLINE | Canned dicts per luminosity segment |
| test_pipeline.py | 16 | OFFLINE | 4 validation star dicts hardcoded from spec (PROXIMA_CEN, SUN_LIKE, SIRIUS_A, DELTA_CEPHEI) -- real astrometry values, no network |
| test_periodogram.py | 14 | OFFLINE | Injected sinusoids + pure noise (make_sinusoidal, make_noise) |
| test_lightcurve.py | 11 | MIXED | Offline: synthetic LCs (_make_lc_data, _make_trended_lc). Network (3 tests): real MAST search/download of KIC 6922244 |
| test_transit.py | 58 | MIXED | Offline: synthetic box and V-shaped transits, even/odd depth cases. Network (1 test): KIC 6922244 BLS recovery of Kepler-410A b |
| test_phase7a.py | 27 | OFFLINE | Synthetic RV sinusoids (200 obs); reads the real in-code TARGET_CATALOG; synthetic Gaia/Hipparcos PM dicts |
| test_phase7b.py | 28 | OFFLINE | Synthetic multi-instrument RV, injection-recovery grids, analytic sensitivity curves |
| test_phase7b2.py | 10 | OFFLINE (RadVel-dependent) | RadVel fits on SYNTHETIC injected planets (10 d/5 m/s single, 17 m/s two-instrument offset, HD 20794-like 3-planet). skipif RadVel missing. No mocking of network -- these modules take arrays as input |
| test_phase7c.py | 28 | OFFLINE (RadVel-dependent) | Synthetic nightly-binning scenarios, injected activity for decorrelation, planet + quasi-periodic activity for GP recovery, vary-flag regression tests |
| test_phase8.py | 43 | OFFLINE | Canned rv_result dicts for scorecards; fake SIMBAD rows for filter_shell_rows; monkeypatched survey_target for failure injection; READS REAL data/known_planets.json (so edits to that file can break assertions -- intentional, it pins the curated values) |
| test_82_eridani_integration.py | 4 | NETWORK | Real SIMBAD + NASA + Gaia/Hipparcos + full report for HD 20794. Slow; results depend on live archive state |

Testing philosophy in this repo:
- Analysis code takes arrays, never fetches -- so tests inject known
  signals and assert recovery (periods within tolerance, K within
  tolerance, offsets recovered). No HTTP mocking layer is needed.
- Network-touching functions are covered only by @pytest.mark.network
  integration tests; everything else runs offline.
- test_phase8.py deliberately reads the real known_planets.json as a
  regression pin on curated literature values.
- RadVel-based tests skip cleanly if RadVel is not installed.

---

## 7. Quick answers

- "Is X fetched live?" -- only via the six APIs in section 2; grep
  for the calling functions listed there.
- "Where do planet parameters come from?" -- data/known_planets.json
  for the curated targets; NASA archive live query otherwise;
  deep_dive._known_planets_fallback as last resort for HD 20794.
- "Where do stellar Teff/mass come from?" -- curated TARGET_CATALOG
  values (literature) for the 10 main targets; spectral-type
  interpolation (approximate) for the 176 shell-catalog stars; the
  Gaia-based pipeline (Modules 1-3) when run explicitly.
- "Can the whole test suite run offline?" -- yes: 253 of 261 tests;
  the 8 network-marked ones are excluded with -m "not network".
- "What must exist on disk for the survey to run?" --
  data/known_planets.json (committed) and, for --shell mode,
  data/shell_catalog.json (committed; regenerable in ~30 s).
