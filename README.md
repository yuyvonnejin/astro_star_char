# Astronomy Star Characterization Pipeline

Compute physical properties of stars -- distance, temperature, luminosity, radius, and mass -- from Gaia DR3 survey data. Optionally retrieve TESS/Kepler light curves for variability analysis and exoplanet transit detection.

## What This Does

Given photometric observations from the [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3) catalog, this pipeline runs five sequential modules:

1. **Distance** (Module 1): Bayesian parallax inversion with an exponentially decreasing space-density prior (Bailer-Jones 2015), or Cepheid period-luminosity relation (Leavitt Law) for classical Cepheids.
2. **Temperature + Luminosity + Radius** (Module 2): Effective temperature from the dereddened (BP-RP) color index using the Mucciarelli et al. (2021) calibration, then bolometric correction, luminosity via the distance modulus, and radius via the Stefan-Boltzmann law.
3. **Mass** (Module 3): Piecewise mass-luminosity relation for main-sequence stars.
4. **Light Curve + Variability** (Module 4, optional): Retrieve TESS/Kepler/K2 light curves from MAST via lightkurve, detect periodic signals with Lomb-Scargle periodogram, classify variability. Cepheid periods feed back into Module 1 for improved distance.
5. **Transit Detection + Planet Characterization** (Module 5, optional): BLS (Box Least Squares) transit search with multi-candidate extraction, then derive planet radius, orbital distance, equilibrium temperature, habitable zone status, and size classification.

```
Gaia DR3 observation
        |
        v
  [Parallax]  ---------->  Module 1: Distance (pc)
        |
        v
  [Color BP-RP]  -------->  Module 2: Teff (K) + Luminosity (Lsun) + Radius (Rsun)
  + distance
        |
        v
  [Luminosity]  ---------->  Module 3: Mass (Msun)
        |
        v  (optional, requires network)
  [MAST light curve]  --->  Module 4: Period + Variability classification
        |                     |
        |                     +---> Cepheid feedback -> recompute distance
        v
  [BLS periodogram]  ---->  Module 5: Transit detection + Planet properties
```

## Project Structure

```
astro_calib/
  src/
    pipeline.py        # Pipeline orchestration + CLI entry point
    distance.py        # Module 1: Bayesian parallax + Cepheid Leavitt law
    temperature.py     # Module 2: Teff, BC_G, luminosity, radius
    mass.py            # Module 3: Mass-luminosity relation
    data_access.py     # Gaia DR3 queries + SIMBAD name resolution
    lightcurve.py      # Module 4a: MAST light curve retrieval + cleaning
    periodogram.py     # Module 4b: Lomb-Scargle period detection + variability
    transit.py         # Module 5: BLS transit detection + planet characterization
  tests/
    test_distance.py
    test_temperature.py
    test_mass.py
    test_pipeline.py       # Integration tests against validation targets
    test_lightcurve.py     # Light curve retrieval + cleaning tests
    test_periodogram.py    # Period detection + variability classification tests
    test_transit.py        # BLS detection + planet property tests
  docs/
    detailed_pipeline_gaia_based.md   # Full technical spec (v3.0)
    phase2_lightcurve_plan.md         # Phase 2 planning document
    phase3_exoplanet_plan.md          # Phase 3 planning document
    good_questions.md                 # Research direction guidance
  logs/                # Runtime logs
  output/              # Pipeline output files
  run_stars.py         # Run pipeline on predefined or custom stars
  tutorial_stellar_properties.ipynb   # Educational walkthrough notebook
  requirements.txt
```

## Setup

```bash
# Create virtual environment
python -m venv venv

# Install dependencies
.\venv\Scripts\pip install -r requirements.txt

# For the tutorial notebook (optional)
.\venv\Scripts\pip install ipykernel matplotlib
.\venv\Scripts\python -m ipykernel install --user --name astro_calib --display-name "Python (astro_calib)"
```

## Usage

### Quick run (predefined stars)

```bash
# Run all predefined stars (Sun, Proxima Cen, Sirius A, Delta Cep, Alpha Cen A, Barnard's Star)
.\venv\Scripts\python run_stars.py

# Run specific stars
.\venv\Scripts\python run_stars.py sun alpha_cen_a
```

To add a new star, add its data dict to the `STARS` dictionary in `run_stars.py`.

### With light curve analysis

```bash
# Enable Module 4: download and analyze MAST light curves
.\venv\Scripts\python run_stars.py proxima_cen --lightcurve

# Enable Module 5: BLS transit detection (auto-enables light curve retrieval)
.\venv\Scripts\python run_stars.py --name "KIC 6922244" --transit

# Both
.\venv\Scripts\python run_stars.py --name "KIC 11497958" --lightcurve --transit
```

### SIMBAD name lookup

Look up any star by SIMBAD identifier -- HD numbers, Bayer names, common names, etc. The name is resolved to a Gaia DR3 source_id, then queried from the Gaia archive automatically.

```bash
# Single star
.\venv\Scripts\python run_stars.py --name "tau Cet"

# Multiple stars
.\venv\Scripts\python run_stars.py --name "HD 22049" "61 Cyg A" "eps Eri"

# Mix predefined + SIMBAD lookups
.\venv\Scripts\python run_stars.py sun --name "tau Cet" "HD 22049"
```

Note: Very bright stars (e.g. Vega, Sirius) may not be in Gaia DR3 due to detector saturation.

### CLI

Process a single star from a JSON file:

```bash
.\venv\Scripts\python -m src.pipeline --json-file input.json
```

Process inline JSON:

```bash
.\venv\Scripts\python -m src.pipeline --json-str '{"source_id":"5853498713190525696","parallax_mas":768.07,"parallax_error_mas":0.03,"phot_g_mean_mag":11.13,"bp_rp":3.50,"ag_gspphot":0.01,"ebpminrp_gspphot":0.005,"feh":0.0,"logg":4.5,"is_cepheid":false,"cepheid_period_days":null}'
```

Write results to a file:

```bash
.\venv\Scripts\python -m src.pipeline --json-file input.json --output output/results.json
```

### As a library

```python
from src.pipeline import process_star

star = {
    "source_id": "5853498713190525696",
    "parallax_mas": 768.07,
    "parallax_error_mas": 0.03,
    "phot_g_mean_mag": 11.13,
    "bp_rp": 3.50,
    "ag_gspphot": 0.01,
    "ebpminrp_gspphot": 0.005,
    "feh": 0.0,
    "logg": 4.5,
    "is_cepheid": False,
    "cepheid_period_days": None,
}

# Basic stellar properties only
result = process_star(star)

# With light curve analysis and transit detection
result = process_star(star, include_lightcurve=True, include_transit=True,
                      lc_target="Proxima Cen")
# result contains: distance_pc, teff_K, luminosity_Lsun, radius_Rsun, mass_Msun,
#   period_days, variability_class, transit_detected, planet_radius_Rearth, ...
```

### Querying Gaia DR3 directly

```python
from src.data_access import query_stars_by_id, query_cepheids

stars = query_stars_by_id([5853498713190525696])
cepheids = query_cepheids(limit=50)
```

## Example Output

Sun-like synthetic star at 10 pc (with reference comparison):

```
  Star                : sun
  Distance            : 10.0023 pc [9.99, 10.00]  (parallax_bayesian)
  Temperature         : 5684 +/- 61 K  (ok)
    (reference)       : 5772 K  (err: 1.5%)
  Luminosity          : 0.87657 Lsun
    (reference)       : 1.00000 Lsun  (err: 12.3%)
  Radius              : 0.9655 Rsun  (ref: 1.0000 Rsun, err: 3.5%)
  Mass                : 0.968 Msun  (ref: 1.000 Msun, err: 3.2%)
  Main sequence       : True
```

Reference values come from the Gaia FLAME module (`radius_flame`, `mass_flame`, `distance_gspphot`, `lum_flame`, `teff_gspphot`) for SIMBAD-resolved stars, or from literature for predefined stars. Percent error is shown when reference data is available.

## Tests

```bash
.\venv\Scripts\python -m pytest tests/ -v
```

96 tests covering all five modules: distance, temperature/luminosity/radius, mass, light curve retrieval, period detection, variability classification, BLS transit detection, planet property derivation, and integration tests against validation targets (Sun, Proxima Centauri, Sirius A, Delta Cephei). Tests requiring network access (MAST queries) are marked with `@pytest.mark.network`.

## Tutorial

The `tutorial_stellar_properties.ipynb` notebook walks through the physics and math behind each module step by step, with visualizations. Open it in Jupyter and select the **Python (astro_calib)** kernel.

## Style Guide

- Keep code simple and concise. No over-engineering.
- Flat is better than nested. Avoid unnecessary abstractions.
- No print statements -- use logging.
- Don't add error handling for scenarios that can't happen.
- Prefer adding a few lines inline over creating a helper for a one-time operation.

## References

- Andrae, R. et al. (2018). Gaia Data Release 2: first stellar parameters from Apsis. A&A 616, A8.
- Bailer-Jones, C. A. L. (2015). Estimating Distances from Parallaxes. PASP 127, 994.
- Fulton, B. J. et al. (2017). The California-Kepler Survey. III. A Gap in the Radius Distribution of Small Planets. AJ 154, 109.
- Kopparapu, R. K. et al. (2013). Habitable Zones around Main-sequence Stars: New Estimates. ApJ 765, 131.
- Kovacs, G., Zucker, S. & Mazeh, T. (2002). A box-fitting algorithm in the search for periodic transits. A&A 391, 369.
- Lightkurve Collaboration (2018). Lightkurve: Kepler and TESS time series analysis in Python. ascl:1812.013.
- Mucciarelli, A., Bellazzini, M. & Massari, D. (2021). Exploiting the Gaia EDR3 photometry to derive stellar temperatures. A&A 653, A90.
- Riess, A. G. et al. (2022). A Comprehensive Measurement of the Local Value of the Hubble Constant. ApJ 934, L7.
- VanderPlas, J. T. (2018). Understanding the Lomb-Scargle Periodogram. ApJS 236, 16.
- Winn, J. N. (2010). Exoplanet Transits and Occultations. In Exoplanets, ed. S. Seager.
