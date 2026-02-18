# Phase 6: Future Directions

## Current State (2026-02-18)

Pipeline covers Modules 1-5: stellar properties (distance, Teff, luminosity, radius, mass)
from Gaia DR3, light curve retrieval + variability analysis, and BLS transit detection with
validation. Benchmarked on 13 TESS confirmed planets.

## Key Findings from Benchmark Testing

### What works
- Stellar property derivation: avg errors ~2-5% for distance, Teff, radius, mass
- Transit detection: 9/13 targets have period-matched candidates (within 5%, incl. aliases)
- Validation (even/odd, shape classification) catches eclipsing binaries

### Known limitations

**Period alias problem (unresolved):**
- 5 of 9 matched detections are at 2x or 0.5x the true period (BLS harmonic aliases)
- Causes ~38% average planet radius error due to distorted depth measurement
- Attempted fix: generate alias candidates at P/2 and 2P (commit 1063d29)
- Result: did not improve -- adding candidates without disambiguation strategy is insufficient
- Functions retained in transit.py but disabled (commit 8c7f533)
- Needs: chi-squared or depth-comparison based disambiguation to select true period

**Earth-analog detection gap:**
- Transit method is fundamentally biased toward short-period, large planets
- Earth-like transit: ~84 ppm depth, ~365d period, 0.5% geometric probability
- No confirmed Earth twin (1 R_earth, 1 AU, G-type star) exists by any method
- Pipeline benchmark targets are all short-period (P < 30d except TOI-2010)
- Current pipeline is best suited for hot/warm planet population

## Future Work Items

### 1. Period Alias Disambiguation (high priority)
The alias problem is the single largest source of error. Approaches to investigate:
- Compare chi-squared of phase-folded fits at P, P/2, 2P
- Use transit shape (ingress/egress duration ratio) to distinguish true period
- Check if odd/even transits at P/2 show alternating depths (indicates 2P is true)
- Literature: Huang et al. (2013), "A Gap in the Mass-Period Distribution"

### 2. Long-Period Transit Detection (medium priority)
People have detected long-period transits -- investigate how:
- Single-transit event detection (mono-transits): Osborn et al. (2016), Gill et al. (2020)
- BLS modifications for few-transit regimes
- Matched-filter approaches for individual transit events
- TESS extended mission enables longer baselines for some targets
- Kepler had 4-year baseline -- revisit Kepler data for long-period candidates

### 3. Literature Review (medium priority)
Go through recent literature to identify pipeline gaps:
- Detrending methods: Gaussian processes (Aigrain et al. 2016), pixel-level decorrelation
- False positive vetting: centroid analysis, nearby star contamination
- Transit parameter refinement: MCMC fitting of transit model (Mandel & Agol 2002)
- Multi-planet system detection: iterative BLS after signal removal
- Stellar contamination: third-light dilution corrections

### 4. Radial Velocity Integration (longer term)
RV is complementary to transit -- mass vs radius gives bulk density:
- Public RV archives: HARPS, HIRES, ESPRESSO
- Combined transit+RV gives composition constraints (rocky vs gaseous)
- Earth induces ~9 cm/s on Sun; current instruments at ~10-30 cm/s (almost feasible)

### 5. Alternative Detection Methods (exploratory)
| Method | What it gives | Data sources |
|--------|--------------|--------------|
| Transit Timing Variations | Hidden planets in known systems | Kepler/TESS timing data |
| Direct Imaging | Planet itself | Limited public data |
| Microlensing | ~1 AU sensitivity | Roman Space Telescope (future) |

## Pickup Point

Next session should start with:
1. Literature review on period alias disambiguation techniques
2. Literature review on long-period / mono-transit detection
3. Then decide which to implement first based on feasibility
