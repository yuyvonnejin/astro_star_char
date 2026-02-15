# Astronomy Object Property Pipeline — Spec v1.1

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
  "lum_gspphot": null
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
| `teff_gspphot` | `gaia_source.teff_gspphot` | Gaia's own T_eff estimate (for validation) |
| `lum_gspphot` | `astrophysical_parameters.lum_flame` | Gaia's own luminosity estimate (for validation) |

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
       ap.mh_gspphot, ap.logg_gspphot, ap.lum_flame
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

The complete output per star combines all three modules:

```json
{
  "source_id": "5853498713190525696",

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

## Dependencies

- **Python 3.9+**
- `astroquery` — Gaia TAP+ data access
- `numpy` — Numerical computation
- `scipy` — Posterior integration for Bayesian distance

---

## Future Improvements (v2.0)

- Integrate `colte` library (Casagrande et al. 2021) for multi-band weighted-average T_eff with Monte Carlo uncertainties.
- Apply the full Lindegren et al. parallax zero-point correction (function of magnitude, color, and sky position) instead of the global +0.017 mas correction.
- Add G-to-V color transformation for more accurate Cepheid distance modulus computation.
- Add period extraction from raw light curves via Lomb-Scargle periodogram.
- Implement comparison against theoretical isochrones for more robust main-sequence classification.

---

## References

- Andrae, R. et al. (2018). Gaia Data Release 2: first stellar parameters from Apsis. A&A 616, A8.
- Bailer-Jones, C. A. L. (2015). Estimating Distances from Parallaxes. PASP 127, 994.
- Casagrande, L. & VandenBerg, D. A. (2018). On the use of Gaia magnitudes and new tables of bolometric corrections. MNRAS 479, L102.
- Casagrande, L. et al. (2021). Effective temperature calibration from the InfraRed Flux Method in the Gaia system. MNRAS 507, 2684.
- Mucciarelli, A., Bellazzini, M. & Massari, D. (2021). Exploiting the Gaia EDR3 photometry to derive stellar temperatures. A&A 653, A90.
- Riess, A. G. et al. (2022). A Comprehensive Measurement of the Local Value of the Hubble Constant. ApJ 934, L7.