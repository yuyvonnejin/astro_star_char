# Exoplanet Detection Methods: Principles, Use Cases, and Limitations

## 1. Radial Velocity (RV)

### Principle
Measures the periodic Doppler shift in a star's spectrum caused by its wobble around the star-planet center of mass. As the star moves toward us, its light blueshifts; as it moves away, it redshifts.

The induced velocity amplitude is:
```
K = (2πG/P)^(1/3) · (M_p·sin(i))/(M_star^(2/3)) · 1/√(1-e²)
```

where P is orbital period, M_p is planet mass, i is orbital inclination, and e is eccentricity.

### Best Use Cases
- **Jupiter analogs** at 1-10 AU around nearby stars (requires multi-year baselines)
- **Mass determination** when combined with transits (breaks the sin(i) degeneracy)
- **Characterizing planetary systems** - can detect multiple planets and measure orbital architectures
- **Low-mass stars** - smaller stars have larger wobbles for the same planet mass

### Limitations
- **sin(i) degeneracy**: Measures M_p·sin(i), not true mass (unless orbit is edge-on)
- **Stellar activity**: Surface features (spots, plages) create velocity signals that mimic planets
- **Long periods require patience**: Jupiter would need 11+ years of observations
- **Massive stars problematic**: Fewer spectral lines, faster rotation blurs lines
- **Precision ceiling**: Currently ~0.1-1 m/s, limiting detection to super-Earths or above at habitable-zone distances

### Current State
Ground-based spectrographs (HARPS, ESPRESSO, NEID) achieve ~0.3 m/s precision. Upcoming instruments aim for 0.1 m/s to access Earth-mass planets in habitable zones around solar-type stars.

---

## 2. Transit Method

### Principle
Detects the periodic dimming when a planet crosses in front of its host star from our viewpoint. The depth of the transit gives the planet-to-star radius ratio:

```
ΔF/F = (R_p/R_star)²
```

Transit duration, ingress/egress times, and timing encode orbital parameters and stellar density.

### Best Use Cases
- **Statistical surveys** - can monitor thousands of stars simultaneously
- **Radius determination** - direct measurement of planetary size
- **Atmospheric characterization** - transmission spectroscopy during transit reveals atmospheric composition
- **Small planets** - Earth-sized planets detectable around solar-type stars, sub-Earths around M dwarfs
- **System architecture** - transit timing variations reveal additional planets

### Limitations
- **Geometric bias**: Only detects planets with orbits nearly edge-on (probability ~ R_star/a for random orientations)
- **False positives**: Eclipsing binaries, blended systems require careful vetting
- **Host star characterization critical**: Planet radius accuracy depends on knowing stellar radius
- **Limited to close-in planets** for space missions (periods < 100 days for Kepler/TESS due to mission lifetimes)
- **Cannot determine mass** without RV follow-up

### Current State
TESS (Transiting Exoplanet Survey Survey) is surveying nearby bright stars. Future missions (PLATO, Ariel) will extend to longer periods and characterize atmospheres.

---

## 3. Transit Timing Variations (TTVs)

### Principle
Measures deviations from strictly periodic transit times caused by gravitational interactions between planets. A planet in the system perturbs the transiting planet's orbit, causing transits to arrive early or late in a characteristic pattern.

The TTV signal encodes information about the perturbing planet's mass, orbital period, and proximity to resonance.

### Best Use Cases
- **Detecting non-transiting planets** - can find planets that don't transit but gravitationally affect those that do
- **Mass determination without RV** - particularly valuable for faint stars where RV is difficult
- **Near-resonant systems** - TTVs are amplified near mean-motion resonances (2:1, 3:2, etc.)
- **Low-mass planets** - can detect sub-Earth mass planets if they're near resonance with a transiting planet

### Limitations
- **Requires multiple planets** - at least one must transit
- **Long baseline needed** - years of observations to detect clear patterns
- **Degeneracies** - multiple orbital configurations can produce similar TTV patterns
- **Model complexity** - N-body simulations required for interpretation
- **Weak signals** for non-resonant configurations

### Current State
Kepler discovered numerous multi-planet systems with detectable TTVs. Ongoing analysis of TESS and future PLATO data will expand this catalog.

---

## 4. Astrometry

### Principle
Measures the positional wobble of a star on the sky caused by orbiting planets. The angular displacement is:

```
α = (M_p/M_star) · (a/d)
```

where a is the semi-major axis and d is the distance to the system.

### Best Use Cases
- **Wide-separation planets** - most sensitive to large semi-major axes (unlike RV which favors close-in planets)
- **Face-on orbits** - complementary to RV which is insensitive to face-on systems
- **Nearby stars** - angular displacement scales inversely with distance
- **Long-period giants** - Jupiter analogs at several AU
- **True mass determination** - no sin(i) degeneracy

### Limitations
- **Precision requirements extreme** - microarcsecond precision needed
- **Favors nearby systems** - signal drops as 1/distance
- **Long baselines required** - need to observe significant fraction of orbit
- **Crowded fields problematic** - need clean astrometric reference frame

### Current State
Gaia is revolutionizing the field, having detected thousands of giant planet candidates at several AU around nearby stars. Ground-based interferometry (GRAVITY) can achieve comparable precision for brightest targets.

---

## 5. Direct Imaging

### Principle
Directly photographs planets by blocking starlight with a coronagraph or using adaptive optics to resolve the planet spatially. Detects reflected light (optical) or thermal emission (infrared).

### Best Use Cases
- **Wide-separation giants** - young, self-luminous planets at > 10 AU
- **Atmospheric characterization** - can take spectra of the planet itself
- **Young systems** - planets are hotter and brighter when young (< 100 Myr)
- **Debris disk systems** - direct imaging can reveal planets sculpting disk structure

### Limitations
- **Extreme contrast** - star outshines planet by factors of 10⁶ to 10¹⁰
- **Small angular separations impossible** - habitable-zone planets around solar-type stars are too close to the star
- **Favors young, massive planets** - old Jupiter-mass planets too faint
- **Distance limitations** - need ~10 AU at ~10 pc for sufficient angular separation

### Current State
Ground-based instruments (SPHERE, GPI, SCExAO) have imaged dozens of planets. Future space missions (HWO, LIFE) aim to image rocky planets in habitable zones using starshades or advanced coronagraphs.

---

## 6. Microlensing

### Principle
Detects planets via gravitational lensing when a foreground star-planet system passes in front of a background star. The planet creates a brief spike or anomaly in the magnification curve.

Most sensitive at the Einstein ring radius:
```
R_E = √(4GM·D_L·(D_S - D_L)/(c²·D_S))
```

typically corresponding to 1-10 AU for galactic bulge events.

### Best Use Cases
- **Cold analogs to solar system planets** - sensitive to Jupiter/Saturn at 1-10 AU
- **Distant stellar populations** - works for stars thousands of light-years away
- **Low-mass planets** - can detect Earth-mass planets at several AU
- **Free-floating planets** - can detect planets without host stars

### Limitations
- **Non-repeating** - each event is unique, no follow-up observations possible
- **Degeneracies** - multiple parameter combinations can fit the same light curve
- **Host star characterization difficult** - lens is typically too distant/faint to study directly
- **Requires precise photometry** of ongoing events in real-time
- **Geometric coincidence** - relies on rare alignment events

### Current State
OGLE and MOA surveys have found several dozen planets. Roman Space Telescope will conduct a major microlensing survey, finding thousands of cold planets and providing statistical census at 1-10 AU.

---

## 7. Eclipse Timing Variations (ETVs)

### Principle
Similar to TTVs but applies to eclipsing binary stars. When a planet orbits one or both stars in a binary (circumbinary planet), it perturbs the binary orbit, causing eclipse times to vary periodically.

### Best Use Cases
- **Circumbinary planets** - planets orbiting both stars (like Tatooine)
- **Planets in binary systems** - where RV is complicated by the presence of two stars
- **Long-term dynamical studies** - can track orbital evolution over decades

### Limitations
- **Requires eclipsing binary** - only ~1% of binaries eclipse from our viewpoint
- **Complex dynamics** - three-body or four-body problem, difficult to model
- **Binary characterization required** - need accurate binary parameters
- **Stability constraints** - circumbinary planets must orbit beyond critical stability radius

### Current State
About a dozen circumbinary planets known, mostly from Kepler. These systems test planet formation theories in complex dynamical environments.

---

## 8. Disk Kinematics

### Principle
Observes kinematic disturbances in protoplanetary disks caused by embedded planets. Planets create gaps, spirals, and velocity perturbations in the gas disk that can be detected via high-resolution spectroscopy of molecular emission lines.

ALMA observations of CO, ¹³CO, and other molecules reveal:
- Kinks in velocity channels indicating local pressure perturbations
- Non-Keplerian velocity patterns from planet-disk gravitational interactions
- Gap edges with perturbed velocity structure

### Best Use Cases
- **Young, actively forming planets** - catches planets in the act of formation
- **Massive planets** - Jupiter-mass or larger create detectable disk perturbations
- **Nearby young systems** - requires high angular resolution (ALMA)
- **Understanding planet-disk interaction** - directly observes the physics of migration

### Limitations
- **Young systems only** - disks dissipate within ~10 Myr
- **Interpretation challenges** - disk features can arise from multiple mechanisms (planets, instabilities, shadows)
- **Angular resolution** - need to resolve disk structure (works best for nearby systems)
- **Ambiguous mass constraints** - difficult to uniquely determine planet mass from disk gaps

### Current State
ALMA has revealed stunning disk structures in HL Tau, PDS 70, and other systems. Some show clear planetary companions (PDS 70 b,c imaged directly) while others show suggestive gaps and spirals.

---

## Summary Table: Complementary Parameter Space Coverage

| Method | Optimal Period Range | Mass Sensitivity | Key Strength | Major Limitation |
|--------|---------------------|------------------|--------------|------------------|
| Radial Velocity | 1 day - 20 years | M·sin(i) ≥ 1 M⊕ (for HZ) | Mass measurement | sin(i) degeneracy |
| Transit | 0.5 - 500 days | Radius ≥ 0.5 R⊕ | Radius + atmosphere | Geometric bias |
| TTV | Variable (coupled to other planets) | Can reach sub-Earth | Detects non-transiting | Requires multi-planet system |
| Astrometry | 1 - 100+ years | Depends on distance | True mass, wide orbits | Precision requirements |
| Direct Imaging | > 10 AU (years - centuries) | > Jupiter mass (young systems) | Atmosphere spectra | Limited to young/wide/massive |
| Microlensing | Months - decades (at detection) | Earth mass at few AU | Cold planets, distant stars | Non-repeating events |
| ETV | Tied to binary period | Variable | Circumbinary systems | Requires eclipsing binary |
| Disk Kinematics | N/A (forming planets) | > Jupiter mass | Planet formation | Young systems only |

The complementary nature of these methods means no single technique gives a complete census. Understanding the generative process behind each detection method's biases is essential for interpreting the exoplanet population as a whole - much as understanding selection effects is crucial in any observational science.