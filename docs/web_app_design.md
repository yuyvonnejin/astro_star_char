# Stellar Property Pipeline -- Web App Design

## Overview

A single-page GitHub Pages web app that lets users enter a star identifier (SIMBAD name or Gaia DR3 source ID) and runs the stellar property pipeline in the browser, displaying intermediate steps for educational purposes.

## Architecture

Pure client-side JavaScript. No server required. Hosted on GitHub Pages.

```
web/
  index.html          -- Single-page app
  css/
    style.css         -- Styling (no framework)
  js/
    app.js            -- Main controller, UI logic
    api.js            -- Gaia TAP+ and SIMBAD REST API calls
    distance.js       -- Module 1: Bayesian parallax / Cepheid Leavitt law
    temperature.js    -- Module 2: Teff, luminosity, radius (Mucciarelli + Andrae)
    mass.js           -- Module 3: Piecewise mass-luminosity relation
    charts.js         -- Plotly.js chart rendering
    examples.js       -- Pre-cached Gaia data for 6 known stars
```

### External Dependencies (CDN)

- **Plotly.js** -- Interactive scientific charts (posterior distribution, HR diagram, mass-luminosity)
- No other frameworks. Vanilla JS + CSS.

## API Strategy

### Gaia TAP+ (ESA Archive)
- Endpoint: `https://gea.esac.esa.int/tap-server/tap/sync`
- Method: POST with `ADQL` query, `FORMAT=json`
- Used for: fetching star data by source_id, Cepheid lookup
- CORS: Public ESA service, expected to support cross-origin

### SIMBAD TAP (CDS Strasbourg)
- Endpoint: `https://simbad.cds.unistra.fr/simbad/sim-tap/sync`
- Used for: name resolution (star name -> coordinates -> Gaia cone search)
- CORS: Public CDS service, expected to support cross-origin

### Fallback
- 6 pre-cached stars always work offline (examples.js)
- If API calls fail due to CORS or network, user is informed and can use examples or manual input

## Page Layout (top to bottom)

### 1. Header
- Title: "Stellar Property Pipeline"
- Subtitle with brief description

### 2. Input Section
- Text field for star name or Gaia DR3 source ID
- "Run Pipeline" button
- Quick-pick buttons: Sun, Proxima Cen, Sirius A, Delta Cephei, Alpha Cen A, Barnard's Star
- Collapsible "Advanced: Manual Input" with fields for all star parameters

### 3. Pipeline Steps (vertical cards, revealed sequentially)

**Step 0: Data Retrieved**
- Table of raw Gaia fields with descriptions
- Highlights which fields feed which module

**Step 1: Distance (Module 1)**
- Method used (Bayesian or Cepheid)
- For Bayesian: formula display, posterior distribution chart (Plotly), MAP + credible interval
- For Cepheid: Leavitt law formula with numbers
- Comparison: naive 1/parallax vs Bayesian

**Step 2: Temperature, Luminosity, Radius (Module 2)**
- Dereddening calculation
- Color-Teff polynomial with coefficients + result
- Bolometric correction (BC_G)
- Luminosity from absolute bolometric magnitude
- Radius from Stefan-Boltzmann law
- Validation ratio vs Gaia estimate

**Step 3: Mass (Module 3)**
- Main-sequence check result
- Mass-luminosity relation chart with piecewise segments (Plotly)
- Which regime this star falls into

### 4. Results Summary
- Clean property card
- Comparison table: Pipeline result vs Gaia estimates vs literature (for known stars)

### 5. HR Diagram
- Teff (reversed x-axis) vs Luminosity (log y-axis)
- Main sequence track as reference
- This star plotted

## Math Porting Notes

All math is straightforward and translatable to vanilla JS:

| Python | JavaScript Equivalent |
|---|---|
| `numpy.logspace` | Manual loop: `Math.pow(10, ...)` |
| `scipy.stats.norm.pdf` | Gaussian formula: `exp(-0.5*((x-mu)/sigma)^2) / (sigma*sqrt(2*pi))` |
| `scipy.integrate.trapezoid` | Manual trapezoidal rule loop |
| `numpy.interp` | Linear interpolation function |
| `numpy.log10` | `Math.log10` |

## Pre-cached Example Stars

The 6 stars from `run_stars.py` are embedded in `examples.js` with their full input data. This guarantees the app works without any API calls.

## Literature Values for Comparison

| Star | Teff (K) | L (Lsun) | R (Rsun) | M (Msun) | Source |
|---|---|---|---|---|---|
| Sun | 5772 | 1.000 | 1.000 | 1.000 | IAU 2015 |
| Proxima Cen | 3042 | 0.0017 | 0.154 | 0.122 | Mann et al. 2015 |
| Sirius A | 9940 | 25.4 | 1.711 | 2.063 | Liebert et al. 2005 |
| Alpha Cen A | 5790 | 1.521 | 1.224 | 1.100 | Kervella et al. 2017 |
| Barnard's Star | 3278 | 0.0035 | 0.187 | 0.144 | Dawson & De Robertis 2004 |
| Delta Cephei | 5900 | 2000 | -- | -- | Turner 2012 |

## Testing Plan

1. Verify all 6 pre-cached stars produce results matching Python pipeline output (manual comparison)
2. Test live API with a known Gaia source ID
3. Test SIMBAD name resolution
4. Test manual input mode
5. Verify all charts render correctly
6. Test on GitHub Pages deployment
