# Astronomy Object Property Pipeline — Spec v3.0

## Overall Goal

The objective is for me to be able to understand how to estimate the characteristics of a star based on astronomical observation. The obvious ones are distance, mass, luminosity, radius, temperature.

Below is the first stage implementation, where input is the Gaia DR3 catalog with pre-processed data.
In the future, I would like to be able to do so directly from light curve. 
Eventually, It would be cool to just specify an coordinate in the sky, resolve to object, fetch light curve data, and generate information about the object. 

For now, below is the stage 1 implmentation. 
Let's not overwrap or complicate the project, but still adhere to a modula, pipeline based coding framework. Keep code simple (don't overthink on exception handling etc.)

## Overview

Given observational data from the Gaia DR3 catalog, compute physical properties of stars: distance, effective temperature, luminosity, radius, and mass.

The pipeline consists of three modules executed sequentially. Distance (Module 1) feeds into temperature/luminosity/radius (Module 2), which feeds into mass (Module 3). Temperature is derived from photometric color; radius is derived from luminosity and temperature via the Stefan-Boltzmann law.

---

## Input Schema

Each star is represented as a single record, sourced from the Gaia DR3 `gaiadr3.gaia_source` table (with optional joins to `gaiadr3.vari_cepheid`).

```json
{
  "source_id": "5853498713190525696",
  "parallax_mas": 768.07,
  "parallax_error_mas": 0.03,
  "phot_g_mean_mag": 11.13,
  "phot_bp_mean_mag": 12.95,
  "phot_rp_mean_mag": 9.45,
  "bp_rp": 3.50,
  "ag_gspphot": 0.01,
  "ebpminrp_gspphot": 0.005,
  "feh": 0.0,
  "logg": 4.5,
  "is_cepheid": false,
  "cepheid_period_days": null,
  "teff_gspphot": null,
  "lum_gspphot": null,
  "ref_radius_Rsun": null,
  "ref_mass_Msun": null,
  "ref_distance_pc": null
}
```

**Field definitions:**

| Field | Source | Description |
|---|---|---|
| `source_id` | `gaia_source.source_id` | Unique Gaia identifier |
| `parallax_mas` | `gaia_source.parallax` | Parallax in milliarcseconds |
| `parallax_error_mas` | `gaia_source.parallax_error` | Parallax uncertainty (mas) |
| `phot_g_mean_mag` | `gaia_source.phot_g_mean_mag` | G-band mean magnitude |
| `phot_bp_mean_mag` | `gaia_source.phot_bp_mean_mag` | BP-band mean magnitude |
| `phot_rp_mean_mag` | `gaia_source.phot_rp_mean_mag` | RP-band mean magnitude |
| `bp_rp` | `gaia_source.bp_rp` | Pre-computed BP - RP color index |
| `ag_gspphot` | `gaia_source.ag_gspphot` | G-band extinction estimate (mag) |
| `ebpminrp_gspphot` | `gaia_source.ebpminrp_gspphot` | E(BP-RP) reddening estimate (mag) |
| `feh` | `astrophysical_parameters.mh_gspphot` or assumed | [Fe/H] metallicity (dex). Default: 0.0 |
| `logg` | `astrophysical_parameters.logg_gspphot` or assumed | Surface gravity log(g) (dex). Default: 4.0 |
| `is_cepheid` | Derived from join to `vari_cepheid` | Boolean flag |
| `cepheid_period_days` | `vari_cepheid.pf` (as 1/pf) | Pulsation period in days, null if not Cepheid |
| `teff_gspphot` | `gaia_source.teff_gspphot` | Gaia's own T_eff estimate (for comparison) |
| `lum_gspphot` | `astrophysical_parameters.lum_flame` | Gaia FLAME luminosity (for comparison) |
| `ref_radius_Rsun` | `astrophysical_parameters.radius_flame` | Gaia FLAME radius in R_sun (for comparison) |
| `ref_mass_Msun` | `astrophysical_parameters.mass_flame` | Gaia FLAME mass in M_sun (for comparison) |
| `ref_distance_pc` | `astrophysical_parameters.distance_gspphot` | Gaia photo-geometric distance in pc (for comparison) |

**Defaults and missing data:**
- If `ag_gspphot` or `ebpminrp_gspphot` is null, assume 0.0 (no extinction correction).
- If `feh` is null, assume 0.0 (solar metallicity).
- If `logg` is null, assume 4.0 (dwarf).

---

## Data Access

Data is queried from the Gaia DR3 archive via the `astroquery.gaia` Python library, which wraps the ESA TAP+ service. Stars can also be resolved by name using SIMBAD (see below).

**SIMBAD name resolution:**

Stars can be looked up by any SIMBAD identifier (HD numbers, Bayer names like "Alp Lyr", common names like "Vega"). The `resolve_simbad_name()` function in `data_access.py` first looks for a Gaia DR3 ID in SIMBAD's identifier list, then falls back to a coordinate-based cone search on Gaia if no direct ID is found.

```python
from src.data_access import resolve_simbad_name, query_stars_by_id

source_id = resolve_simbad_name("tau Cet")  # returns "2452378776434477184"
stars = query_stars_by_id([source_id])
```

**Example query for specific stars:**

```python
from astroquery.gaia import Gaia

query = """
SELECT g.source_id, g.parallax, g.parallax_error,
       g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
       g.bp_rp, g.teff_gspphot,
       g.ag_gspphot, g.ebpminrp_gspphot,
       ap.mh_gspphot, ap.logg_gspphot, ap.lum_flame,
       ap.radius_flame, ap.mass_flame, ap.distance_gspphot
FROM gaiadr3.gaia_source AS g
LEFT JOIN gaiadr3.astrophysical_parameters AS ap
  ON g.source_id = ap.source_id
WHERE g.source_id IN (
    5853498713190525696,
    5853498713160606720
)
"""

job = Gaia.launch_job(query)
result = job.get_results()
```

**Example query for Cepheids:**

Note: Gaia's TAP+ service uses ADQL, which requires `TOP N` instead of `LIMIT N`.

```python
query = """
SELECT TOP 100 g.source_id, g.parallax, g.parallax_error,
       g.phot_g_mean_mag, g.bp_rp,
       g.ag_gspphot, g.ebpminrp_gspphot,
       c.pf, c.type_best_classification
FROM gaiadr3.gaia_source AS g
JOIN gaiadr3.vari_cepheid AS c
  ON g.source_id = c.source_id
WHERE c.type_best_classification = 'DCEP'
"""
```

**Quality filters** (recommended): `parallax_over_error > 10` for high-confidence parallax measurements.

**Parallax zero-point:** Gaia DR3 has a known parallax zero-point bias of approximately -0.017 mas. For v1.0, apply a global additive correction of +0.017 mas to all parallax values before processing.

---

## Module 1: Distance

Computes distance in parsecs from either parallax (Method A) or the Cepheid period-luminosity relation (Method B).

### Method A: Bayesian Parallax Inversion

Used for all non-Cepheid stars with measured parallax.

**Why not naive inversion:** The naive estimate `d = 1000 / parallax_mas` is biased for noisy measurements due to the nonlinear transformation and the asymmetric error propagation. A Bayesian approach with an appropriate prior yields more accurate distances.

**Prior:** Exponentially decreasing space density (Bailer-Jones 2015):

```
P(r) proportional to r^2 * exp(-r / L)
```

where `r` is distance in parsecs and `L = 1350 pc` (the scale length used for the Gaia catalog).

**Likelihood:**

```
P(parallax_obs | r) = Normal(1000/r, parallax_error_mas)
```

where `parallax_obs` is the observed parallax in mas (after zero-point correction).

**Posterior:**

```
P(r | parallax_obs) proportional to r^2 * exp(-r/L) * exp(-0.5 * ((parallax_obs - 1000/r) / parallax_error_mas)^2)
```

**Implementation:** Evaluate the unnormalized posterior numerically on a grid of `r` values from 1 pc to 100,000 pc (log-spaced), normalize, and report the mode as the point estimate. Report the 16th and 84th percentiles as the credible interval.

**Output:** `distance_pc`, `distance_lower_pc`, `distance_upper_pc`, `distance_method = "parallax_bayesian"`

### Method B: Cepheid Period-Luminosity (Leavitt Law)

Used when `is_cepheid = true` and `cepheid_period_days` is not null.

**Calibration:** Riess et al. (2022), V-band Leavitt law:

```
M_V = -2.43 * log10(cepheid_period_days) - 4.05
```

**Distance modulus:**

```
G_0 = phot_g_mean_mag - ag_gspphot
distance_pc = 10^((G_0 - M_V + 5) / 5)
```

Note: This uses G_0 as a proxy for V-band apparent magnitude, which introduces a small systematic error dependent on the star's color. For v1.0 this is acceptable; a future improvement would apply a G-to-V color transformation.

**Output:** `distance_pc`, `distance_method = "cepheid_leavitt"`

### Module 1 Output

```json
{
  "distance_pc": 1.301,
  "distance_lower_pc": 1.298,
  "distance_upper_pc": 1.304,
  "distance_method": "parallax_bayesian"
}
```

---

## Module 2: Temperature, Luminosity, and Radius

Computes effective temperature from the dereddened (BP-RP) color index, then derives bolometric luminosity using the distance from Module 1, and stellar radius from luminosity and temperature via the Stefan-Boltzmann law.

### Step 1: Deredden

```
(BP-RP)_0 = bp_rp - ebpminrp_gspphot
G_0 = phot_g_mean_mag - ag_gspphot
```

### Step 2: Effective Temperature

**Source:** Mucciarelli, Bellazzini & Massari (2021), A&A 653, A90, Table 1.

**Formula:**

```
theta = b0 + b1*C + b2*C^2 + b3*[Fe/H] + b4*[Fe/H]^2 + b5*[Fe/H]*C
T_eff = 5040 / theta
```

where `C = (BP-RP)_0` and `[Fe/H] = feh`.

**Dwarf stars** (logg >= 3.0):

Valid color range: (BP-RP)_0 in [0.39, 1.50]. Dispersion: 61 K.

| b0 | b1 | b2 | b3 | b4 | b5 |
|---|---|---|---|---|---|
| 0.4929 | 0.5092 | -0.0353 | 0.0192 | -0.0020 | -0.0395 |

**Giant stars** (logg < 3.0):

Valid color range: (BP-RP)_0 in [0.33, 1.81]. Dispersion: 83 K.

| b0 | b1 | b2 | b3 | b4 | b5 |
|---|---|---|---|---|---|
| 0.5323 | 0.4775 | -0.0344 | -0.0110 | -0.0020 | -0.0009 |

**Flags:**
- If (BP-RP)_0 is outside the valid range, flag the star and report T_eff as unreliable.
- Report the calibration dispersion (61 K or 83 K) as the temperature uncertainty floor.

### Step 3: Absolute G Magnitude

```
M_G = G_0 - 5 * log10(distance_pc) + 5
```

### Step 4: Bolometric Correction BC_G

**Source:** Andrae et al. (2018), Gaia DR2 FLAME pipeline (Equation 8.9, Table 8.3).

**Formula:**

```
BC_G = a0 + a1*(T_eff - 5772) + a2*(T_eff - 5772)^2
     + a3*(T_eff - 5772)^3 + a4*(T_eff - 5772)^4
```

**For 4000 K <= T_eff <= 8000 K:**

| a0 | a1 | a2 | a3 | a4 |
|---|---|---|---|---|
| 6.000e-02 | 6.731e-05 | -6.647e-08 | 2.859e-11 | -7.197e-15 |

**For 3300 K <= T_eff < 4000 K:**

| a0 | a1 | a2 | a3 | a4 |
|---|---|---|---|---|
| 1.749e+00 | 1.977e-03 | 3.737e-07 | -8.966e-11 | -4.183e-14 |

**Flags:** Stars with T_eff outside 3300-8000 K should be excluded from luminosity computation.

### Step 5: Luminosity

```
M_bol = M_G + BC_G
L_over_Lsun = 10^((4.74 - M_bol) / 2.5)
```

where M_bol_sun = 4.74 (IAU 2015 B2 resolution).

### Step 6: Radius

Stellar radius is derived from luminosity and temperature using the Stefan-Boltzmann law:

```
R/R_sun = sqrt(L/L_sun) * (T_sun / T_eff)^2
```

where T_sun = 5772 K. This requires no additional data beyond what Module 2 already computes.

### Step 7: Validation Cross-Check

If `lum_gspphot` is available, compute:

```
lum_ratio = L_over_Lsun / lum_gspphot
```

Report this ratio. Values between 0.8 and 1.2 indicate good agreement.

### Module 2 Output

```json
{
  "bp_rp_0": 3.495,
  "teff_K": 3042,
  "teff_uncertainty_K": 83,
  "teff_flag": "outside_valid_range",
  "M_G": 11.85,
  "BC_G": 1.62,
  "M_bol": 13.47,
  "luminosity_Lsun": 0.00155,
  "radius_Rsun": 0.154,
  "luminosity_validation_ratio": 0.97
}
```

---

## Module 3: Mass

Estimates stellar mass using the mass-luminosity relation. Only valid for main-sequence stars.

### Main-Sequence Check

A star is considered main-sequence if its derived T_eff and luminosity place it in a region consistent with the main sequence. For v1.0, use a simplified bounding-box check:

```
is_main_sequence = (logg >= 3.5) AND (T_eff >= 3300 K) AND (T_eff <= 8000 K)
```

If `is_main_sequence` is false, flag the star and skip mass estimation.

### Mass-Luminosity Relation

**Formula:**

```
L/L_sun = k * (M/M_sun)^alpha
```

Inverted:

```
M/M_sun = (L_over_Lsun / k)^(1/alpha)
```

**Piecewise coefficients** (k = 1 for all segments in the standard formulation):

| Luminosity range | alpha | Approximate mass range |
|---|---|---|
| L < 0.033 L_sun | 2.3 | M < 0.43 M_sun |
| 0.033 <= L < 16 L_sun | 4.0 | 0.43 <= M < 2 M_sun |
| 16 <= L < 54000 L_sun | 3.5 | 2 <= M < 55 M_sun |
| L >= 54000 L_sun | 1.0 | M >= 55 M_sun |

The luminosity boundaries above are derived from the mass boundaries by applying each segment's power law: L(0.43) = 0.43^4.0 = 0.034; L(2) = 2^3.5 = 11.3; L(55) = 55^1.0 = 55. For implementation, use the luminosity directly to select the segment:

```python
if L < 0.033:
    alpha = 2.3
elif L < 16:
    alpha = 4.0
elif L < 54000:
    alpha = 3.5
else:
    alpha = 1.0

mass_Msun = L ** (1.0 / alpha)
```

**Flags:**
- If `is_main_sequence` is false, set `mass_Msun = null` and `mass_flag = "not_main_sequence"`.
- The mass-luminosity relation has intrinsic scatter of roughly 10-20%. Report this as a caveat.

### Module 3 Output

```json
{
  "mass_Msun": 0.12,
  "mass_flag": "ok",
  "is_main_sequence": true
}
```

---

## Pipeline Output

The complete output per star combines all three modules plus reference data:

```json
{
  "source_id": "5853498713190525696",

  "teff_gspphot": 3042.0,
  "lum_gspphot": 0.00155,
  "ref_radius_Rsun": 0.1542,
  "ref_mass_Msun": 0.1221,
  "ref_distance_pc": 1.301,

  "distance_pc": 1.301,
  "distance_lower_pc": 1.298,
  "distance_upper_pc": 1.304,
  "distance_method": "parallax_bayesian",

  "bp_rp_0": 3.495,
  "teff_K": 3042,
  "teff_uncertainty_K": 83,
  "teff_flag": "outside_valid_range",

  "M_G": 11.85,
  "BC_G": 1.62,
  "M_bol": 13.47,
  "luminosity_Lsun": 0.00155,
  "radius_Rsun": 0.154,
  "luminosity_validation_ratio": 0.97,

  "mass_Msun": 0.12,
  "mass_flag": "ok",
  "is_main_sequence": true
}
```

The `ref_*` and `*_gspphot` fields are official Gaia values for comparison against the pipeline's computed estimates. For SIMBAD-resolved stars these come from the Gaia FLAME module; for predefined stars they are literature values.

Output format: JSON (one object per star) or CSV for batch processing.

---

## Validation Targets

The pipeline should be tested against these well-characterized stars:

| Star | Type | Expected T_eff (K) | Expected L (L_sun) | Expected R (R_sun) | Expected M (M_sun) | Notes |
|---|---|---|---|---|---|---|
| Sun | G2V dwarf | 5772 | 1.0 | 1.0 | 1.0 | Fundamental reference. Use known values directly. |
| Alpha Centauri A | G2V dwarf | ~5790 | ~1.52 | ~1.22 | ~1.10 | Nearby, well-measured binary. |
| Proxima Centauri | M5.5V dwarf | ~3042 | ~0.0017 | ~0.15 | ~0.12 | Nearest star. Tests cool end of calibration. |
| Sirius A | A1V dwarf | ~9940 | ~25.4 | ~1.71 | ~2.06 | Tests hot boundary (outside BC_G range, should be flagged). |
| Delta Cephei | F5Ib Cepheid | ~5500-6800 (variable) | ~2000 | ~44 | ~4.5 | Tests Cepheid distance method. |

---

## Module 4: Light Curve Analysis (v2.0)

Retrieves space-based time-series photometry (TESS, Kepler, K2) from MAST via lightkurve, detects periodic signals, and classifies stellar variability. This module is optional and requires network access.

### Module 4a: Light Curve Retrieval

**Data source:** MAST (Mikulski Archive for Space Telescopes) via the `lightkurve` Python package.

| Mission | Sky coverage | Cadence | Catalog ID |
|---------|-------------|---------|------------|
| TESS | ~85% of sky | 2 min, 10 min, 30 min | TIC |
| Kepler | ~116 sq deg (Cygnus-Lyra) | 1 min, 30 min | KIC |
| K2 | ~20 fields along ecliptic | 1 min, 30 min | EPIC |

**Search priority:** SPOC (TESS), then Kepler, then K2.

**Steps:**
1. Search MAST for available light curves by star name or coordinates.
2. Download up to `max_sectors` sectors/quarters (default: 20).
3. Stitch multi-sector data into a single time series.
4. Remove NaN values and >5-sigma outliers.
5. Flatten with Savitzky-Golay filter (window_length=1001 cadences) to remove spacecraft systematics.

### Module 4b: Period Detection

**Method:** Lomb-Scargle periodogram (astropy implementation).

- Frequency range: 1/max_period to 1/min_period (default: 0.01-50 cycles/day).
- Oversample factor: 5.
- Maximum period capped at half the time baseline.
- Reports: best period, normalized power at peak, False Alarm Probability (FAP).
- Amplitude estimated from binned phase-folded light curve (peak-to-peak, parts per thousand).

**Note:** Lomb-Scargle is optimized for sinusoidal signals. Non-sinusoidal signals (eclipsing binaries, transits) may produce the strongest peak at a harmonic (P/2, P/4) rather than the fundamental period. This is expected behavior.

### Module 4c: Variability Classification

Simple rule-based classification:

| Class | Criterion |
|-------|----------|
| `periodic` | FAP < 0.001 |
| `possible_periodic` | 0.001 <= FAP < 0.01 |
| `non_variable` | FAP >= 0.01 |

Quality flags: `ok`, `few_points` (n < 100), `short_baseline` (< 1 day).

### Cepheid Period Feedback

When Module 4 detects a significant period on a star flagged as a Cepheid (`is_cepheid=True`), the detected period is fed back into Module 1 to recompute distance via the period-luminosity relation. The recomputed distance is reported as `distance_pc_lc`.

### Module 4 Output

```json
{
  "lightcurve_available": true,
  "lc_mission": "TESS Sector 11",
  "lc_author": "SPOC",
  "lc_n_sectors": 7,
  "lc_n_points": 314004,
  "lc_time_baseline_days": 1499.84,
  "lc_cadence_s": 20.0,
  "variability_class": "periodic",
  "variability_flag": "ok",
  "period_days": 3.5225,
  "period_fap": 1.2e-15,
  "amplitude_ppt": 2.5
}
```

### Usage

```bash
# Enable light curve analysis with --lightcurve flag
python run_stars.py proxima_cen --lightcurve

# SIMBAD stars with light curve
python run_stars.py --name "KIC 6922244" --lightcurve
```

```python
# Library API
from src.pipeline import process_star

result = process_star(star_dict, include_lightcurve=True, lc_target="Proxima Cen")
```

---

## Module 5: Transit Detection and Planet Characterization (v4.0)

Uses BLS (Box Least Squares) to detect transit signals in light curves, then combines transit observables with Phase 1 stellar properties to derive physical properties of the orbiting planet.

### Scientific Rationale: Focus on Quiet Sun-Like Stars

Transit photometry has a fundamental signal-to-noise problem: a Jupiter-sized planet transiting a Sun-like star produces a ~1% (10,000 ppm) dip, but an Earth-sized planet produces only ~0.01% (84 ppm). Detecting shallow transits requires the host star's intrinsic variability to be well below the transit depth.

The pipeline therefore adopts a **quiet-star-first strategy**:

1. **Gate on variability:** Module 4 measures stellar variability amplitude. Stars with amplitude exceeding ~10 ppt (parts per thousand) are flagged `too_variable` and skipped for transit analysis. This threshold roughly corresponds to a Jupiter-transit depth around a Sun-like star -- if the star's own variability exceeds this, shallow transit detection is unreliable.

2. **Focus on habitable zone:** For scientifically high-value targets, the BLS search can be narrowed to the habitable zone period range, improving sensitivity where it matters most.

3. **Validate against false positives:** Eclipsing binaries (EBs) are the dominant source of transit-like false positives. The pipeline uses even/odd depth comparison, transit shape classification, and planet-level sanity checks to flag likely EBs.

This approach prioritizes **reliable detections on amenable targets** over attempting to extract marginal signals from noisy data.

### Module 5a: BLS Transit Detection

**Method:** BLS periodogram via astropy `BoxLeastSquares` (Kovacs, Zucker & Mazeh 2002).

BLS is purpose-built for finding periodic box-shaped dips -- exactly what transits look like. This complements Lomb-Scargle (Module 4b) which is optimized for sinusoidal signals.

**Core parameters:**
- Period range: 0.5 to 100 days (capped at half the time baseline).
- Duration grid: 20 linearly-spaced values from 0.01 to 0.5 days (capped at 0.9x min period).
- Significance: Signal Detection Efficiency (SDE) = (max_power - mean) / std. SDE >= 6 for detection.
- Reports: period, depth, duration, epoch, SDE, number of observed transits.

#### Log-Uniform Period Grid

Standard BLS implementations use a frequency-uniform period grid, which packs most points at short periods (e.g., 100k of 200k points below 1 day). This biases detection toward spurious short-period signals from spacecraft systematics.

The pipeline uses a **log-uniform (geometric) period grid** of 50,000 points, giving equal density per decade of period. This provides balanced sensitivity across the full period range -- a transit at 50 days gets the same frequency resolution as one at 0.5 days.

#### Stratified Multi-Candidate Extraction

Rather than reporting only the single highest BLS peak (which is often a short-period systematic), the pipeline extracts **one candidate per period decade** using per-bin statistics:

1. Divide the period range into ~5 logarithmic bins (3 bins per decade, minimum 3 bins total).
2. Within each bin, find the highest-power peak.
3. Compute **local (per-bin) SDE**: compare the peak power to the mean and std of power within that bin only.
4. Skip bins where the local best peak has SDE < 3.0 (no significant signal in that period range).
5. Return all surviving candidates, sorted by descending local SDE.

**Why local SDE matters:** Global SDE compares every peak to the mean/std of the full power spectrum. When a dominant noise peak at P~0.5d inflates the global std, real transit signals at P~20d appear insignificant even though they are strong relative to their local noise floor. Local SDE removes this cross-contamination between period regimes.

#### Candidate Sanity Checks

Each candidate undergoes automated sanity checks before ranking:

| Flag | Trigger | Interpretation |
|------|---------|---------------|
| `duration_exceeds_15pct_of_period` | Transit duration > 15% of orbital period | Geometrically implausible for planets |
| `fewer_than_3_transits` | `floor(baseline / period) < 3` | Too few transits for reliable BLS fit |
| `very_deep_transit` | Depth > 5% (50,000 ppm) | Likely eclipsing binary, not planet |
| `negative_depth` | BLS returns negative depth | Non-physical; artifact |

A candidate with any flag gets `transit_flag` set to the flag string(s); candidates without flags get `transit_flag="ok"`.

#### Candidate Selection Logic

The best candidate is selected by preferring **clean (unflagged) candidates over flagged ones**, even if the flagged candidate has higher SDE. A physically plausible detection at lower SDE is more trustworthy than a high-SDE detection that fails basic sanity checks. Among clean candidates, the highest SDE wins.

### Module 5b: Planet Property Derivation

Given transit parameters + stellar properties from Modules 2-3:

| Property | Formula | Units |
|----------|---------|-------|
| Planet radius | Rp = R_star * sqrt(depth) | R_earth, R_jup |
| Orbital distance | a = (M_star)^(1/3) * (P_yr)^(2/3) | AU |
| Equilibrium temp | T_eq = T_star * sqrt(R_star/(2a)) * (1-A)^(1/4) | K |
| Insolation | S = L_star / a^2 | S_earth |
| Habitable zone | inner = sqrt(L/1.107), outer = sqrt(L/0.356) | AU |

Default Bond albedo: 0.3 (Earth-like).

**Planet size classification** (Fulton et al. 2017):

| Class | Radius range |
|-------|-------------|
| sub_earth | < 0.8 R_earth |
| earth_sized | 0.8 - 1.25 R_earth |
| super_earth | 1.25 - 2.0 R_earth |
| sub_neptune | 2.0 - 4.0 R_earth |
| neptune_sized | 4.0 - 6.0 R_earth |
| sub_jupiter | 6.0 - 10.0 R_earth |
| jupiter_sized | 10.0 - 15.0 R_earth |
| super_jupiter | >= 15.0 R_earth |

**Planet-level sanity checks:** After deriving properties, automatic flags are raised for implausible results:

| Flag | Trigger | Interpretation |
|------|---------|---------------|
| `orbit_inside_star` | Orbital distance (AU) < stellar radius (AU) | Non-physical orbit; period too short or stellar properties wrong |
| `planet_larger_than_half_star` | Planet radius > 50% of stellar radius | Likely eclipsing binary, not a planet |
| `extreme_temperature` | Equilibrium temperature > 4000 K | Planet would be destroyed at this temperature |

These flags are appended to `planet_flag` alongside any inherited transit-level flags.

**When stellar properties are missing:** Transit detection still works (reports period, depth, duration, Rp/R_star ratio). Planet radius in physical units and orbital properties require R_star and M_star from Phase 1.

**Important caveat:** All detections are transit *candidates*. Photometry alone cannot distinguish planet transits from eclipsing binaries or background eclipsing binaries. Confirmation requires follow-up observations.

### Module 5 Output

```json
{
  "transit_detected": true,
  "transit_period_days": 3.5225,
  "transit_depth": 0.0012,
  "transit_depth_ppm": 1200,
  "transit_depth_raw_ppm": 1350,
  "transit_duration_hours": 3.2,
  "transit_epoch": 2458325.123,
  "transit_sde": 15.3,
  "n_transits_observed": 10,
  "transit_flag": "ok",
  "detection_method": "bls",
  "variability_removed_period_days": null,
  "planet_radius_Rearth": 2.3,
  "planet_radius_Rjup": 0.205,
  "planet_radius_ratio": 0.034641,
  "planet_size_class": "sub_neptune",
  "orbital_semi_major_axis_AU": 0.048,
  "orbital_period_days": 3.5225,
  "equilibrium_temp_K": 1250,
  "insolation_Searth": 180.5,
  "hz_conservative_inner_AU": 0.9505,
  "hz_conservative_outer_AU": 1.6761,
  "hz_optimistic_inner_AU": 0.7503,
  "hz_optimistic_outer_AU": 1.7678,
  "in_habitable_zone": false,
  "planet_flag": "ok",
  "depth_even_ppm": 1180.5,
  "depth_odd_ppm": 1220.3,
  "depth_ratio_even_odd": 0.967,
  "even_odd_validation_pass": true,
  "even_odd_flag": "ok",
  "shape_class": "U_shape",
  "flat_bottom_fraction": 0.42,
  "shape_flag": "ok",
  "transit_candidates": [
    {"rank": 1, "transit_period_days": 3.5225, "transit_sde": 15.3, "transit_depth_ppm": 1200, "planet_size_class": "sub_neptune", "in_habitable_zone": false},
    {"rank": 2, "transit_period_days": 48.7, "transit_sde": 7.2, "transit_depth_ppm": 450, "planet_size_class": "super_earth", "in_habitable_zone": true}
  ]
}
```

### Module 5c: Transit Validation & HZ-Targeted Detection (v4.0)

Added in Phase 4 to improve detection reliability and false positive rejection. These features work together as a layered defense: filter unsuitable stars first, focus the search, then validate the detections.

#### Quiet-Star Filter (Gate)

Gates Module 5 based on Module 4's variability amplitude. Highly variable stars drown shallow transits, so transit analysis is skipped unless forced.

- **Threshold:** `max_variability_ppt=10.0` (configurable). Stars above this are flagged `too_variable`.
- **Override:** `force_transit=True` bypasses the check.
- **Missing data:** If amplitude is unavailable, transit analysis proceeds (conservative).
- **Implementation:** `_check_quiet_star()` in `pipeline.py`.

#### Pre-Whitening (Variability Removal)

When Module 4 detects a significant periodic signal (e.g., starspot rotation, pulsations), the pipeline removes it before running BLS to prevent stellar variability from masking or mimicking transit signals.

**Method:** Phase-folded template subtraction (`remove_variability()` in `transit.py`):

1. Phase-fold the light curve at the detected variability period.
2. Divide the phase into 50 bins and compute the median flux per bin.
3. Interpolate the binned template at each data point's phase (with wrap-around for smooth edges).
4. Subtract the template from the flux and re-center on the original median.

This produces a residual light curve where the periodic stellar signal is removed but transit dips are preserved. The method is model-free (no assumed waveform shape), making it robust for arbitrary variability patterns.

**Integration:** Only applied when `variability_class="periodic"` in Module 4 output and the detected period is available. The variability period is passed from `pipeline.py` to `analyze_transit()`.

#### HZ-Targeted BLS Mode

Narrows the BLS period search to the habitable zone for improved sensitivity.

- Computes optimistic HZ in AU from stellar luminosity, converts to period via Kepler's 3rd law.
- Broadens by configurable factor (default 2.0x) to account for stellar property uncertainties (~20-30%).
- Floor at 0.5 days minimum period.
- Falls back to full [0.5, 100] day range if stellar properties are missing.
- CLI: `--hz-only`, `--hz-broadening`
- **Implementation:** `compute_hz_period_range()` in `transit.py`.

#### Even/Odd Transit Validation

Compares transit depths in even vs odd numbered transits. Real planets produce equal depths at every epoch; eclipsing binaries often show different primary/secondary eclipse depths (the secondary eclipse is typically shallower when both stars differ in size/temperature).

**Algorithm** (`validate_even_odd()` in `transit.py`):

1. Assign epoch numbers to all time points: `epoch_num = round((time - epoch) / period)`.
2. Identify in-transit points (within half the transit duration of each transit center).
3. Compute out-of-transit baseline from median of non-transit points.
4. Split in-transit points by even/odd epoch number.
5. Measure median depth for even and odd groups independently.
6. Compute depth ratio = depth_even / depth_odd.

**Pass criterion:** Depth ratio within [0.5, 2.0].
- Generous tolerance because shallow transits are noisy and small-number statistics can produce scatter.
- EBs typically show ratios of 3-10x (very different primary vs secondary eclipse depths).
- Requires at least 5 in-transit data points per group (even and odd).

**Reports:** `depth_even_ppm`, `depth_odd_ppm`, `depth_ratio_even_odd`, `even_odd_validation_pass`, `n_even`, `n_odd`, `even_odd_flag`.

#### Transit Shape Classification (V vs U)

Classifies transit shape to distinguish planetary (U-shaped, flat bottom) from eclipsing binary (V-shaped, triangular) signals.

**Physics:** A planet fully crossing its star produces a U-shaped transit: flat ingress, flat bottom (full occultation), flat egress. A grazing eclipsing binary produces a V-shaped (triangular) eclipse because neither star fully occults the other.

**Algorithm** (`classify_transit_shape()` in `transit.py`):

1. Phase-fold the light curve at the transit period.
2. Select points within 1.5x transit duration of the transit center.
3. Bin finely (100 bins default) within the window.
4. Measure flat-bottom fraction: number of bins with flux within 20% of the minimum depth, divided by total valid bins.

**Classification thresholds:**

| Flat-bottom fraction | Classification | Interpretation |
|---------------------|---------------|---------------|
| > 0.3 | `U_shape` | Planet-like transit with well-defined flat bottom |
| < 0.15 | `V_shape` | Grazing eclipsing binary or blended EB |
| 0.15 - 0.3 | `ambiguous` | Insufficient evidence to classify |

Requires at least 20 in-transit data points and 10 valid bins.

#### Depth Refinement for Multi-Candidate Systems

In systems with multiple transit candidates, the BLS-measured depth for each candidate is contaminated by transits from other planets overlapping in the folded light curve.

**Method** (`refine_transit_depth()` in `transit.py`):

1. Start with the best-ranked candidate.
2. Mask in-transit points from all **clean** (unflagged) other candidates by replacing them with the median flux. Flagged candidates are skipped -- they are likely spurious, and masking them (especially short-period false positives) can corrupt >30% of data points.
3. Phase-fold the masked light curve at the best candidate's period.
4. Re-measure depth as the difference between out-of-transit median and in-transit median.
5. If the refined depth is <= 0, the "transit" was likely an artifact of overlapping signals. Demote this candidate and try the next one.

**Fallback logic:** If depth refinement fails for the top candidate (refined depth <= 0), the pipeline iterates through lower-ranked candidates. The first candidate with a positive refined depth is promoted to #1. This catches cases where the highest-SDE peak is actually an alias of another candidate's signal.

**Output fields:** `transit_depth_raw_ppm` (original BLS depth) and `transit_depth_ppm` (refined depth). Planet properties are recomputed from the refined depth.

#### Candidate Re-Ranking

After enriching all candidates with planet properties (Module 5b), the pipeline re-ranks them by scientific priority:

**Ranking tiers (highest to lowest priority):**
1. **Clean + HZ:** Candidates with no flags AND in the habitable zone (most scientifically valuable).
2. **Clean:** Candidates with no flags but outside the HZ.
3. **Flagged:** Candidates with any sanity check flags.

Within each tier, candidates are sorted by descending SDE. If re-ranking changes the #1 candidate, top-level transit/planet results are updated to reflect the new best.

**Implementation:** `_rerank_candidates()` in `transit.py`.

#### Extended Output Fields

```json
{
  "quiet_star_skip_reason": null,
  "variability_removed_period_days": 5.234,
  "transit_depth_raw_ppm": 1350.0,
  "depth_even_ppm": 1180.5,
  "depth_odd_ppm": 1220.3,
  "depth_ratio_even_odd": 0.967,
  "even_odd_validation_pass": true,
  "even_odd_flag": "ok",
  "shape_class": "U_shape",
  "flat_bottom_fraction": 0.42,
  "shape_flag": "ok",
  "transit_candidates": [
    {"rank": 1, "transit_period_days": 3.52, "transit_sde": 15.3, "in_habitable_zone": false},
    {"rank": 2, "transit_period_days": 48.7, "transit_sde": 7.2, "in_habitable_zone": true}
  ]
}
```

### Usage

```bash
# Enable transit detection (auto-enables light curve retrieval)
python run_stars.py --name "KIC 6922244" --transit

# Both light curve analysis and transit detection
python run_stars.py --name "KIC 6922244" --lightcurve --transit

# HZ-targeted transit detection
python run_stars.py --name "KIC 6922244" --transit --hz-only

# Force transit analysis on variable stars
python run_stars.py --name "KIC 11497958" --transit --force-transit

# Custom variability threshold
python run_stars.py --name "KIC 6922244" --transit --max-variability-ppt 20.0
```

```python
# Library API
from src.pipeline import process_star

result = process_star(star_dict, include_transit=True, lc_target="KIC 6922244")

# HZ-targeted mode
result = process_star(star_dict, include_transit=True, lc_target="KIC 6922244",
                      hz_targeted=True, hz_broadening=2.0)

# Force transit on variable star
result = process_star(star_dict, include_transit=True, lc_target="KIC 11497958",
                      force_transit=True)
```

---

## Dependencies

- **Python 3.9+**
- `astroquery` — Gaia TAP+ data access
- `numpy` — Numerical computation
- `scipy` — Posterior integration for Bayesian distance
- `lightkurve` — MAST light curve access, stitching, flattening (v2.0)

---

## Future Improvements (v5.0)

- Integrate `colte` library (Casagrande et al. 2021) for multi-band weighted-average T_eff with Monte Carlo uncertainties.
- Apply the full Lindegren et al. parallax zero-point correction (function of magnitude, color, and sky position) instead of the global +0.017 mas correction.
- Add G-to-V color transformation for more accurate Cepheid distance modulus computation.
- Implement comparison against theoretical isochrones for more robust main-sequence classification.
- Harmonic analysis to recover fundamental periods from non-sinusoidal signals.
- Gyrochronology: estimate stellar age from detected rotation periods.
- Batch light curve processing for multiple stars.
- Transit timing variations (TTV) for multi-planet system detection.
- Limb darkening modeling for improved transit depth measurement.
- ~~Multi-planet detection (search for secondary signals after removing primary transit).~~ (Done in v4.0: multi-candidate BLS with stratified extraction + depth refinement by masking other candidates)
- ~~Eclipsing binary vs. planet transit discrimination heuristics.~~ (Done in v4.0: even/odd validation + transit shape classification + planet-level sanity checks)

---

## References

- Andrae, R. et al. (2018). Gaia Data Release 2: first stellar parameters from Apsis. A&A 616, A8.
- Bailer-Jones, C. A. L. (2015). Estimating Distances from Parallaxes. PASP 127, 994.
- Casagrande, L. & VandenBerg, D. A. (2018). On the use of Gaia magnitudes and new tables of bolometric corrections. MNRAS 479, L102.
- Casagrande, L. et al. (2021). Effective temperature calibration from the InfraRed Flux Method in the Gaia system. MNRAS 507, 2684.
- Lightkurve Collaboration (2018). Lightkurve: Kepler and TESS time series analysis in Python. ascl:1812.013.
- Mucciarelli, A., Bellazzini, M. & Massari, D. (2021). Exploiting the Gaia EDR3 photometry to derive stellar temperatures. A&A 653, A90.
- Fulton, B. J. et al. (2017). The California-Kepler Survey. III. A Gap in the Radius Distribution of Small Planets. AJ 154, 109.
- Kopparapu, R. K. et al. (2013). Habitable Zones around Main-sequence Stars: New Estimates. ApJ 765, 131.
- Kovacs, G., Zucker, S. & Mazeh, T. (2002). A box-fitting algorithm in the search for periodic transits. A&A 391, 369.
- Riess, A. G. et al. (2022). A Comprehensive Measurement of the Local Value of the Hubble Constant. ApJ 934, L7.
- VanderPlas, J. T. (2018). Understanding the Lomb-Scargle Periodogram. ApJS 236, 16.
- Winn, J. N. (2010). Exoplanet Transits and Occultations. In Exoplanets, ed. S. Seager.