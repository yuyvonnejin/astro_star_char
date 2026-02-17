# Phase 3: Exoplanet Characterization from Transits -- Plan

## Goal

Add Module 5 (Transit Detection and Planet Characterization) that uses BLS to find transit signals in light curves, then combines the transit observables with Phase 1 stellar properties to derive physical properties of the orbiting planet: radius, orbital distance, equilibrium temperature, and habitable zone status.

This is the phase where Phase 1 (stellar characterization) and Phase 2 (light curves) **depend on each other** -- you need the star's radius to get the planet's radius, and the star's mass to get the orbit.

## The Science in Brief

When a planet crosses in front of its host star (a "transit"), it blocks a fraction of the star's light. From that dip:

```
What we measure from the light curve:
  - Period (P): time between transits
  - Transit depth (delta): fractional flux drop = (Rp / R_star)^2
  - Transit duration (T_dur): total time of the dip

What we need from Phase 1 (stellar properties):
  - R_star: stellar radius (from Module 2)
  - M_star: stellar mass (from Module 3)
  - T_star: effective temperature (from Module 2)
  - L_star: luminosity (from Module 2)

What we derive (planet properties):
  - Rp: planet radius = R_star * sqrt(delta)
  - a: orbital semi-major axis from Kepler's 3rd law
  - T_eq: equilibrium temperature
  - Habitable zone: is the planet between the inner and outer HZ boundaries?
```

## Module 5 Design

### 5a: BLS Transit Detection (`src/transit.py`)

BLS (Box Least Squares) is purpose-built for finding periodic box-shaped dips -- exactly what transits look like. This replaces Lomb-Scargle for transit-specific detection and directly fixes the harmonic problem we saw with KIC 6922244.

lightkurve has BLS built in:
```python
pg = lc.to_periodogram(method='bls',
                        minimum_period=0.5,
                        maximum_period=max_period)
best_period = pg.period_at_max_power
transit_depth = pg.depth_at_max_power      # fractional flux dip
transit_duration = pg.duration_at_max_power # hours
transit_time = pg.transit_time_at_max_power # epoch T0
```

**Steps**:
1. Run BLS periodogram on the flattened light curve.
   - Period range: 0.5 to 100 days (most known transiting planets).
   - Duration grid: 0.01 to 0.5 days (15 min to 12 hours).
2. Identify the strongest BLS peak.
3. Compute signal detection efficiency (SDE) or BLS power as significance metric.
   - SDE > 6-7 is a common threshold for candidate detection.
4. Phase-fold at the detected period and verify the dip shape.
5. Measure transit parameters: depth, duration, epoch.
6. Compute the number of observed transits: `n_transits = floor(baseline / period)`.

**Output**:
```python
{
    "transit_detected": True,
    "transit_period_days": 3.5225,
    "transit_depth": 0.0012,          # fractional (Rp/R*)^2
    "transit_depth_ppm": 1200,        # parts per million
    "transit_duration_hours": 3.2,
    "transit_epoch_btjd": 2458325.123,
    "transit_sde": 15.3,             # signal detection efficiency
    "n_transits_observed": 10,
    "detection_method": "bls",
}
```

### 5b: Planet Property Derivation (`src/transit.py`)

Given transit parameters + stellar properties, derive the planet's physical characteristics.

**Planet radius**:
```
Rp / R_star = sqrt(delta)
Rp = R_star * sqrt(delta)

Rp_earth = Rp * R_sun_earth    (R_sun = 109.076 R_earth)
```

Size classification (approximate, from Fulton et al. 2017):
| Category | Radius range |
|----------|-------------|
| Sub-Earth | Rp < 0.8 R_earth |
| Earth-sized | 0.8 <= Rp < 1.25 R_earth |
| Super-Earth | 1.25 <= Rp < 2.0 R_earth |
| Sub-Neptune | 2.0 <= Rp < 4.0 R_earth |
| Neptune-sized | 4.0 <= Rp < 6.0 R_earth |
| Sub-Jupiter | 6.0 <= Rp < 10.0 R_earth |
| Jupiter-sized | 10.0 <= Rp < 15.0 R_earth |
| Super-Jupiter | Rp >= 15.0 R_earth |

**Orbital semi-major axis** (Kepler's 3rd law):
```
a^3 = G * M_star * P^2 / (4 * pi^2)

In convenient units:
a_AU = (M_star / M_sun)^(1/3) * (P_years)^(2/3)
     = (M_star / M_sun)^(1/3) * (P_days / 365.25)^(2/3)
```

**Equilibrium temperature**:
```
T_eq = T_star * sqrt(R_star / (2 * a)) * (1 - A_bond)^(1/4)
```

Where A_bond is the Bond albedo. Default assumption: A_bond = 0.3 (Earth-like).
In solar units:
```
T_eq = T_star * sqrt(R_star_Rsun * R_sun_AU / (2 * a_AU)) * (1 - A_bond)^(1/4)

R_sun_AU = 0.00465047 AU
```

**Habitable zone boundaries** (Kopparapu et al. 2013, moist greenhouse / maximum greenhouse):
```
Conservative HZ:
  inner_AU = sqrt(L_star / 1.107)   (moist greenhouse limit)
  outer_AU = sqrt(L_star / 0.356)   (maximum greenhouse limit)

Optimistic HZ:
  inner_AU = sqrt(L_star / 1.776)   (recent Venus limit)
  outer_AU = sqrt(L_star / 0.320)   (early Mars limit)
```

Where L_star is in solar luminosities. A planet is "in the habitable zone" if `inner_AU <= a_AU <= outer_AU`.

**Insolation flux** (stellar energy received by the planet, relative to Earth):
```
S = L_star / a_AU^2    (in units of S_earth)
```

**Output**:
```python
{
    "planet_radius_Rearth": 2.3,
    "planet_radius_Rjup": 0.205,
    "planet_size_class": "sub_neptune",
    "orbital_semi_major_axis_AU": 0.048,
    "orbital_period_days": 3.5225,
    "equilibrium_temp_K": 1250,
    "insolation_Searth": 180.5,
    "hz_inner_AU": 0.95,
    "hz_outer_AU": 1.67,
    "in_habitable_zone": False,
    "planet_flag": "ok",
}
```

### When Phase 1 Data Is Missing

The planet characterization requires stellar radius and mass from Phase 1. If these are unavailable (e.g., star outside calibration range), we can still report:
- Transit detection: yes/no, period, depth, duration (from light curve alone)
- Planet radius: only as `Rp/R_star` ratio (dimensionless)
- Flag: `"stellar_params_missing"` -- planet properties are incomplete

This keeps the module useful even when Phase 1 can't fully characterize the host star.

---

## How It Integrates

```
Star identifier
    |
    +---> Phase 1: [Module 1-3: Distance, Teff, L, R, Mass]
    |                    |
    |                    +---> R_star, M_star, T_star, L_star
    |                                   |
    +---> Phase 2: [Module 4: Light Curve + Period Detection]
    |                    |
    |                    +---> flattened light curve
    |                                   |
    +---> Phase 3: [Module 5a: BLS Transit Detection]
                         |
                         +---> period, depth, duration
                                        |
                         [Module 5b: Planet Properties] <--- R_star, M_star, T_star, L_star
                                        |
                                        +---> Rp, a, T_eq, HZ status
```

### Pipeline changes

- `src/transit.py` (NEW): BLS detection + planet property derivation.
- `src/periodogram.py` (EDIT): Add `detect_transit_bls()` function alongside existing `detect_period()`.
  - Or: keep BLS entirely in `transit.py` since it serves a different purpose.
  - Decision: **put BLS in transit.py** -- conceptually it belongs to exoplanet detection, not general variability.
- `src/pipeline.py` (EDIT): Add Module 5 after Module 4, passing stellar properties.
  - New parameter: `include_transit=False` (separate from `include_lightcurve`).
  - When `--lightcurve --transit` are both set, run both analyses.
  - When `--transit` is set without `--lightcurve`, auto-enable light curve retrieval.
- `run_stars.py` (EDIT): Add `--transit` flag, display formatting for planet results.

### Relationship between Module 4 and Module 5

Module 4 (Lomb-Scargle) and Module 5a (BLS) serve different purposes:
- **Module 4**: General variability -- pulsations, rotation, eclipsing binaries. Reports: period, amplitude, variability class.
- **Module 5a**: Transit-specific detection. Reports: transit period, depth, duration, SDE.

Both run on the same flattened light curve. A star could be both a rotational variable AND have a transiting planet -- these are independent signals.

In practice: if `--transit` is set, we run BLS in addition to (not instead of) Lomb-Scargle.

---

## Implementation Steps

### Step 1: Implement BLS transit detection in `src/transit.py`

Functions:
- `detect_transit(time, flux, flux_err, min_period=0.5, max_period=100, duration_range=(0.01, 0.5))` -- BLS search.
- Returns transit parameters or `{"transit_detected": False}`.

Uses lightkurve's BLS:
```python
import lightkurve as lk

lc_obj = lk.LightCurve(time=time, flux=flux, flux_err=flux_err)
pg = lc_obj.to_periodogram(method='bls',
                            minimum_period=min_period,
                            maximum_period=max_period)
```

SDE (signal detection efficiency) computation:
```python
sde = (pg.max_power - np.mean(pg.power)) / np.std(pg.power)
```

### Step 2: Implement planet property derivation in `src/transit.py`

Functions:
- `compute_planet_properties(transit_result, stellar_props)` -- derives Rp, a, T_eq, HZ.
  - `stellar_props` is a dict with keys: `radius_Rsun`, `mass_Msun`, `teff_K`, `luminosity_Lsun`.
- `classify_planet_size(radius_Rearth)` -- returns size class string.
- `compute_habitable_zone(luminosity_Lsun)` -- returns inner/outer boundaries.
- `compute_equilibrium_temp(teff_K, radius_Rsun, a_AU, albedo=0.3)` -- returns T_eq.

Constants needed:
```python
R_SUN_REARTH = 109.076     # R_sun in Earth radii
R_SUN_AU = 0.00465047      # R_sun in AU
R_EARTH_RJUP = 0.08921     # R_earth in Jupiter radii
```

### Step 3: Integrate into pipeline

- Update `src/pipeline.py`:
  - Add `include_transit=False` parameter to `process_star()`.
  - Module 5 runs after Module 4 (needs the flattened light curve).
  - Pass stellar properties from Modules 2-3 into Module 5b.
- Update `run_stars.py`:
  - Add `--transit` flag.
  - Add display formatting for planet properties (radius, orbit, temperature, HZ).

### Step 4: Write unit tests

**`tests/test_transit.py`** (offline, synthetic signals):
- Test BLS detection on synthetic box-shaped transit signal.
  - Generate: constant flux with periodic box dips.
  - Verify: correct period, depth, duration recovered.
- Test planet radius calculation: known depth + known R_star -> correct Rp.
- Test Kepler's 3rd law: known period + known M_star -> correct semi-major axis.
- Test equilibrium temperature calculation against known values.
- Test habitable zone boundaries for Sun (should give ~0.95 - 1.67 AU).
- Test planet size classification across all categories.
- Test graceful handling when stellar properties are missing.

**Network tests** (marked `@pytest.mark.network`):
- Test BLS on KIC 6922244: should detect transit at P ~ 3.5225 d.

### Step 5: Validate with known exoplanet hosts

| Target | Planet | Expected | Why |
|--------|--------|----------|-----|
| KIC 6922244 | Kepler-410A b | P~3.52 d, Rp~2.8 Re | Our existing test case; BLS should find fundamental period |
| Kepler-10 | Kepler-10b | P~0.84 d, Rp~1.47 Re | First confirmed rocky exoplanet, well-characterized |
| HAT-P-7 | HAT-P-7b | P~2.20 d, Rp~1.36 Rj | Hot Jupiter, deep transit, easy detection |
| TOI-700 | TOI-700d | P~37.4 d, Rp~1.07 Re | TESS habitable zone planet |

### Step 6: Update documentation

- Update `docs/detailed_pipeline_gaia_based.md` with Module 5 spec.
- Update `docs/phase2_lightcurve_plan.md` noting that BLS moved from "future" to implemented.

---

## File Changes Summary

| File | Action | Description |
|------|--------|-------------|
| `src/transit.py` | NEW | BLS transit detection + planet property derivation |
| `src/pipeline.py` | EDIT | Add Module 5 integration, `include_transit` param |
| `src/periodogram.py` | MINOR | No changes needed (LS stays for variability) |
| `run_stars.py` | EDIT | Add `--transit` flag, planet display formatting |
| `tests/test_transit.py` | NEW | Transit detection + planet property tests |
| `docs/detailed_pipeline_gaia_based.md` | EDIT | Add Module 5 spec, update to v3.0 |

---

## Constants and References

**Physical constants used:**
```
R_sun = 6.957e8 m = 109.076 R_earth = 0.00465047 AU
M_sun = 1.989e30 kg
R_earth = 6.371e6 m
R_jupiter = 7.149e7 m = 11.209 R_earth
AU = 1.496e11 m
G = 6.674e-11 m^3 kg^-1 s^-2
M_bol_sun = 4.74
T_sun = 5772 K
```

**References:**
- Kovacs, G., Zucker, S. & Mazeh, T. (2002). BLS algorithm. A&A 391, 369.
- Kopparapu, R. K. et al. (2013). Habitable zone boundaries. ApJ 765, 131.
- Fulton, B. J. et al. (2017). The radius gap. AJ 154, 109.
- Winn, J. N. (2010). Exoplanet Transits and Occultations. In Exoplanets, ed. S. Seager.

---

## Scope Control

**In scope for Phase 3:**
- BLS transit detection
- Planet radius, orbital distance, equilibrium temperature, HZ check
- Size classification
- Insolation flux

**Not in scope:**
- Transit timing variations (TTV) for multi-planet systems
- Limb darkening modeling
- Planet mass estimation (requires radial velocity, not available from photometry alone)
- Atmospheric characterization
- Multi-planet detection in the same light curve (detect only the strongest signal)
- Eclipsing binary vs planet discrimination (flag as candidate, don't adjudicate)

**Known limitation:** We cannot distinguish a planet transit from an eclipsing binary or a background eclipsing binary (BEB) from photometry alone. All detections should be labeled as "candidates" not "confirmed planets." Confirmation requires follow-up (radial velocity, high-resolution imaging, etc.).
