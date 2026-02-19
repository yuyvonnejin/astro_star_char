# Phase 7: Exo-Earth Targeted Search -- Proof of Concept

**Date: 2026-02-19**
**Prerequisite: Phases 1-6 complete (stellar calibration + transit detection pipeline)**
**Nature: Research PoC -- defining the search space, not yet automating it**

---

## 1. Motivation and Scope

### 1.1 Why Pivot from Transit Pipeline Improvements

Phases 1-6 built a functional stellar characterization and transit detection pipeline. Benchmark testing on 13 TESS confirmed planets around quasi-solar stars demonstrated:
- Strong stellar property derivation (2-5% errors)
- Reliable transit detection for short-period planets (9/13 period-matched)
- A fundamental detection gap: transit surveys are systematically biased toward short-period, large-radius planets

The original goal -- detecting planetary systems truly like ours -- falls outside what transit photometry alone can achieve with current data. Rather than incrementally improving a method with diminishing returns toward this goal, Phase 7 pivots to a **multi-method, proximity-ordered, targeted search** of the nearest Sun-like stars.

### 1.2 What This Phase Is

This is a **Proof of Concept** phase. The primary deliverables are not planet detections (though any would be remarkable) but rather:

1. **Search hyperparameters**: What stellar selection criteria, distance limits, data quality thresholds, and method combinations define a tractable and scientifically productive search?
2. **Success criteria**: What constitutes a meaningful result? What detection limits can we achieve per target? What would a non-detection tell us?
3. **Methodological inventory**: For each target, what data exists, what methods can be applied, and what are the sensitivity limits of each?
4. **Tool proficiency**: Hands-on experience with RV archives, astrometric data, FFI photometry, and their respective analysis techniques -- understanding what each tool can and cannot do before attempting automation.

Research is fundamentally about clearing out the unknown. This phase maps the unknown before trying to traverse it.

---

## 2. Relationship to Existing and Upcoming Programs

A key question: where does this project sit relative to major funded efforts targeting the same scientific goal?

### 2.1 Landscape Assessment

| Program | Timeline | Method | Overlap with this project |
|---------|----------|--------|--------------------------|
| **Terra Hunting Experiment** | Starting 2026; 10-year campaign | RV (HARPS3, 50 Sun-like stars, nightly cadence) | **Complementary upstream**: will produce the definitive RV dataset. This project can analyze archival precursor data (HARPS, ESPRESSO) now, then incorporate Terra Hunting results as they become available. |
| **PLATO** | Launch late 2026 / early 2027 | Transit (200,000+ stars, HZ-focused) | **Complementary future data**: will extend transit sensitivity to Earth-analog periods. This project's transit pipeline could be adapted for PLATO targets. |
| **Gaia DR4** | Expected December 2026 | Astrometry (approximately 7,500 predicted planet detections) | **Complementary data source**: Gaia detects planets at wider separations (2-5 AU), not HZ. Useful for understanding full system architecture of our target stars. |
| **ESPRESSO GTO** | Ongoing | RV (approximately 10 cm/s) | **Direct precursor data**: ESPRESSO has already observed several of our target stars (tau Ceti, 82 G. Eridani, alpha Centauri). Public data releases inform our sensitivity analysis. |
| **Planet Hunters TESS** | Ongoing (Zooniverse) | Visual transit identification | **Parallel but different scope**: carpet search of all TESS light curves by citizen scientists. Finds unusual signals pipelines miss, but no target prioritization. |
| **Exoplanet Watch** (NASA) | Ongoing | Ground-based transit refinement | **Subsequent**: refines known transiting planets; not discovery-focused. |

### 2.2 This Project's Niche

This project is **not competing with** these programs. Instead, it occupies a specific niche:

- **Target-focused synthesis**: Major programs specialize in one method applied broadly. This project applies *all available methods* to a small, curated set of the very nearest and most Sun-like stars -- a depth-first rather than breadth-first approach.
- **Archival data integration**: Substantial public data already exists (HARPS RV, TESS photometry, Gaia astrometry) but has rarely been analyzed holistically per-target in a unified framework.
- **Preparatory work for future data**: By establishing per-target baselines and sensitivity maps now, this project is positioned to immediately incorporate Terra Hunting, PLATO, and Gaia DR4 data as they become available.
- **Methodological learning**: Understanding the limits and complementarity of each method on real targets is intrinsically valuable and transferable.

---

## 3. Target Selection

### 3.1 Why Exclude Hot Stars: The Case Against High-Temperature Hosts

The upper Teff boundary at approximately 6000 K (G0V / F9V) is driven primarily by a hard observational wall: above the **Kraft break** at approximately 6200 K, stars lose their deep convective envelopes, stop magnetically braking, and rotate rapidly (20-100+ km/s vs 2-5 km/s for the Sun). This destroys precision radial velocity capability -- the single most important detection method for HZ planets around nearby stars.

**Detection degradation with increasing Teff:**

| Factor | Effect as Teff increases | Critical threshold |
|--------|--------------------------|-------------------|
| **RV precision collapse** | Rapid stellar rotation broadens spectral lines, reducing RV information content by 100-1000x. ESPRESSO achieves approximately 10 cm/s on G dwarfs but only approximately 100-500 m/s on mid-F stars. Since RV is central to multi-method HZ planet detection, this is disqualifying. | **approximately 6200 K** (Kraft break). Above this, virtually all stars are fast rotators. The 6000 K upper bound stays safely below this wall. |
| **Fewer spectral lines** | Hotter stars have fewer and broader metal absorption lines (spectrum dominated by hydrogen Balmer and helium lines). This degrades both RV precision and the ability to measure stellar metallicity and detailed atmospheric parameters. | Gradual degradation above approximately 6500 K; severe by A-type (approximately 8000 K+). |
| **Shallower transits** | Larger stellar radius means smaller (Rp/R\*)^2. Earth transiting: G2V approximately 84 ppm; F5V approximately 43 ppm; A0V approximately 13 ppm. At 43 ppm, Earth transits are below TESS single-transit detection for most targets. | approximately 6500 K (F5V, R approximately 1.4 Rsun). Marginal for Earth-sized detection. |
| **Stellar pulsations** | Many F-type stars are delta Scuti pulsators (approximately 6300-8600 K, photometric variations 1-300 mmag) or gamma Doradus pulsators (approximately 6900-7500 K, multi-day variations). These add both photometric noise masking transits and RV noise masking planets. | approximately 6300 K. Pulsation incidence rises sharply into the F range. |
| **Shorter main-sequence lifetime** | F5V: approximately 4 Gyr; F0V: approximately 2 Gyr; A0V: approximately 0.5-1 Gyr. Complex life on Earth took approximately 4 Gyr to develop. Stars with MS lifetimes below this threshold may not host biologically interesting planets. | approximately 6500 K (F5V, approximately 4 Gyr). Fatal by approximately 7200 K (F0V, approximately 2 Gyr). |
| **HZ moves outward** | Higher luminosity pushes the HZ to wider orbits (F5V: HZ at approximately 1.5-2.5 AU). This reduces geometric transit probability (P_transit = R\*/a) and weakens the RV signal for a given planet mass (K scales as P^(-1/3)). | Gradual; no hard threshold. Each factor degrades by approximately 30-50% from G2V to F5V. |
| **Different space weather** | Above the Kraft break, the absence of a solar-type magnetic dynamo produces weaker, differently structured stellar winds. The planetary magnetosphere interaction and atmospheric erosion regime are fundamentally non-solar. | approximately 6200 K. |

**Summary of upper boundary temperatures:**

| Boundary | Teff | Spectral type | What changes |
|----------|------|--------------|-------------|
| **6000 K** (adopted upper bound) | G0V / F9V | Just below the Kraft break. All targets are slow rotators amenable to precision RV. Ensures multi-method approach is viable. |
| **6200 K** (Kraft break) | F8V | Above this: rapid rotation destroys RV precision. Most stars rotate at 20+ km/s. |
| **6500 K** | F5V | Above this: MS lifetime drops below approximately 4 Gyr. Transit depth for Earth drops below approximately 43 ppm. Delta Scuti pulsations become common. |
| **7200 K** | F0V | Above this: MS lifetime approximately 2 Gyr. Planet detection by any spectroscopic method becomes extremely difficult. |

The current target list includes Zeta Tucanae (F9.5V, 5956 K) and Eta Cassiopeiae A (G0V, 6012 K) right at this boundary -- both are slow rotators with usable RV data, confirming the 6000 K cut is appropriately placed.

### 3.2 Why Exclude Cool Stars: The Case Against Low-Temperature Hosts

The lower Teff boundary is not arbitrary. As the host star cools below solar temperatures, multiple compounding effects make an "Earth-like" outcome progressively less likely and harder to detect. The boundary at approximately 5200 K (G/K spectral type divide) reflects the balance between sample size and scientific relevance to our specific goal: finding systems where a planet could experience conditions similar to Earth's.

**Habitability degradation with decreasing Teff:**

| Factor | Effect as Teff decreases | Critical threshold |
|--------|--------------------------|-------------------|
| **Habitable zone moves inward** | HZ at approximately 0.4-0.8 AU for K5V (4400 K) vs approximately 0.95-1.7 AU for the Sun. Closer orbits mean stronger tidal forces, higher stellar wind flux, and a fundamentally different orbital environment. | Below approximately 4800 K (mid-K), HZ planets are likely **tidally locked** (synchronous rotation), producing permanent day/night sides -- a climate regime with no solar system analog. |
| **Stellar magnetic activity** | Cooler main-sequence stars have deeper convective envelopes, driving stronger magnetic dynamos. This produces more frequent and energetically significant flares, especially relative to the star's bolometric luminosity. | M dwarfs (Teff < 3700 K) are worst; late K (4000-4800 K) are moderate; early K and G are relatively quiet. The activity-luminosity ratio starts rising noticeably below approximately 5000 K. |
| **XUV radiation and atmospheric stripping** | Extreme UV and X-ray flux (relative to total luminosity) increases for cooler stars. This photoevaporates planetary atmospheres over Gyr timescales. Planets in the HZ of cooler stars receive higher XUV doses per unit of visible light. | Below approximately 4500 K, the cumulative XUV fluence over a planet's lifetime can strip Earth-like atmospheres unless the planet has strong magnetic field protection. |
| **Pre-main-sequence desiccation** | Lower-mass stars have longer pre-main-sequence phases where they are more luminous than their MS state. A planet that ends up in the MS-era HZ was inside the PMS-era HZ, potentially losing its volatiles. | PMS duration: approximately 30 Myr for 1 Msun (G2V), approximately 100 Myr for 0.7 Msun (K5V), approximately 1 Gyr for 0.1 Msun (M5V). The K/M boundary at approximately 3700 K is extreme, but even mid-K stars have 2-3x longer PMS volatile-loss windows. |
| **Spectral energy distribution** | Cooler stars shift emission toward infrared. This changes atmospheric photochemistry (less UV for ozone production), alters greenhouse gas behavior (different absorption band overlaps), and may inhibit prebiotic chemistry pathways that require UV photons. | Below approximately 5000 K, the UV flux drops sharply. Some origin-of-life models require UV for RNA synthesis (Ranjan & Sasselov 2016). |
| **Planet formation efficiency** | Protoplanetary disk mass correlates with stellar mass. Lower-mass stars form less massive disks, potentially producing fewer or smaller rocky planets in the HZ. | Kepler statistics show rocky planet occurrence *increases* for K dwarfs, but the giant planet occurrence decreases. The net effect on system architecture (i.e., whether Solar System-like configurations form) is unclear below approximately 0.8 Msun. |
| **RV detection challenges** | Cooler stars have stronger molecular bands (TiO, VO for late K/M) that complicate spectral fitting. Stellar activity induces larger RV jitter relative to the planet signal. | Below approximately 4500 K, molecular contamination becomes significant. However, K dwarfs in the 4800-5200 K range are often called the "Goldilocks zone" for RV because they have more spectral lines than G stars with less activity than M stars. |
| **Transit detection** | Actually *improves* for cooler stars -- smaller R_star means deeper transit for same-sized planet, and shorter HZ periods mean more transits per year. | This is the one factor that favors cool stars. Earth transiting a K5V star produces approximately 300 ppm (vs 84 ppm for the Sun). |
| **Luminosity evolution stability** | G-type stars have nearly constant luminosity over their main-sequence lifetime (approximately 10-30% variation over 10 Gyr). K dwarfs are even more stable. M dwarfs have longer lifetimes but their early high-luminosity phase is a concern. | Actually favors K dwarfs. The HZ of a K star is more stable over Gyr timescales than for a G star. |
| **Departure from "like ours"** | The Sun is G2V (5772 K). As Teff decreases, the stellar environment becomes progressively less solar: different UV flux, different HZ geometry, different tidal regime, different activity pattern. The goal of this project is specifically to find systems that resemble ours. | This is a goal-dependent criterion, not a physical one. For a different goal (e.g., "maximize habitable planet occurrence rate"), K dwarfs might be *preferred*. |

**Summary of boundary temperatures:**

| Boundary | Teff | Spectral type | What changes |
|----------|------|--------------|-------------|
| **5200 K** (current lower bound) | G9V / K0V | Below this: stellar activity rises, UV flux drops, protoplanetary disk mass decreases, departure from "solar" increases substantially. This is the G/K divide. |
| **4800 K** (relaxed bound) | K2V / K3V | Below this: tidal locking becomes probable for HZ planets. Molecular bands start complicating RV analysis. Still scientifically productive for planet detection, but "like ours" becomes a stretch. |
| **4400 K** (extended bound) | K5V | Below this: HZ planets are almost certainly tidally locked. PMS desiccation window exceeds 200 Myr. System architecture unlikely to resemble the Solar System. |
| **3700 K** (M dwarf boundary) | K7V / M0V | Below this: extreme flaring, very close-in HZ, strong atmospheric stripping. Excellent for detection (deep transits, large RV signal) but very different habitability regime. |

**This PoC adopts 5200 K as the baseline lower bound** (G stars only), consistent with the goal of finding "systems like ours." However, Section 3.2 below explores how relaxing this boundary affects the available target pool.

### 3.3 Selection Tiers: Strictness vs. Sample Size Trade-off

The tighter the match to solar properties, the fewer candidates exist nearby. This PoC needs approximately 10 targets to develop methodology, which constrains how strict we can be. The table below shows four tiers of increasing strictness and the distance radius needed to reach 10 candidates at each tier.

**Tier definitions:**

| Tier | Name | Teff range | Mass range | Additional constraints | G-dwarf density |
|------|------|-----------|------------|----------------------|-----------------|
| **A** | Strict Solar Twin | 5672-5872 K (within 100 K of Sun) | 0.95-1.05 Msun | [Fe/H] within 0.1 dex of solar; log g > 4.3 | Extremely rare |
| **B** | Solar Analog | 5500-6000 K (G0-G5) | 0.85-1.15 Msun | Main sequence only | approximately 0.003 per pc^3 |
| **C** | Broad G-type | 5200-6000 K (G0-G9) | 0.70-1.15 Msun | Main sequence only | approximately 0.006 per pc^3 |
| **D** | Extended (G + early K) | 4800-6000 K (G0-K2) | 0.65-1.15 Msun | Main sequence only | approximately 0.010 per pc^3 |

**Available targets at each distance:**

| Distance | Volume (pc^3) | Tier A (Solar Twin) | Tier B (Solar Analog) | Tier C (Broad G) | Tier D (G + early K) |
|----------|--------------|--------------------|-----------------------|------------------|---------------------|
| 5 pc | 524 | 0 | approximately 1 (alpha Cen A) | approximately 2 (+ tau Ceti) | approximately 3-4 |
| 10 pc | 4,189 | 0-1 | approximately 8-10 | approximately 19-21 | approximately 32-39 |
| 15 pc | 14,137 | approximately 2-3 | approximately 25-35 | approximately 50-60 | approximately 85-110 |
| 20 pc | 33,510 | approximately 5-8 | approximately 50-70 | approximately 120-150 | approximately 200-260 |
| 25 pc | 65,450 | approximately 10-15 | approximately 100-140 | approximately 230-280 | approximately 390-500 |
| 50 pc | 523,599 | approximately 10-20 confirmed | -- | -- | -- |

Sources: RECONS 10 pc census (2018.3); Reyle et al. (2021) Gaia 10 pc sample; Porto de Mello et al. (2014) solar twin survey within 50 pc; density extrapolations for larger volumes.

**Distance needed for 10 targets at each tier:**

| Tier | Distance for approximately 10 targets | Implication |
|------|---------------------------------------|-------------|
| **A** (Solar Twin) | approximately 25-50 pc | Very few exist; must go distant, reducing sensitivity of all methods. Solar twins are genuinely rare. |
| **B** (Solar Analog) | approximately 10 pc | Just barely achievable within 10 pc. Good balance: close enough for high sensitivity, strict enough to be scientifically focused. |
| **C** (Broad G) | approximately 7-8 pc | Easily achievable. approximately 19-21 G dwarfs within 10 pc provides a comfortable pool. |
| **D** (G + early K) | approximately 5-6 pc | Very tight distance cut possible. approximately 32-39 candidates within 10 pc. Includes stars where HZ planets may be tidally locked. |

**Recommendation for the PoC:**

Use **Tier C (Broad G-type, Teff 5200-6000 K) within 10 pc** as the primary selection, giving approximately 19-21 candidates from which to choose the best 10 based on data availability and binary status. This provides enough targets without needing to relax to K dwarfs or extend to large distances.

However, include **18 Scorpii** (14.1 pc, Tier A solar twin) and **HD 134060** (24.2 pc, Tier B solar analog with confirmed planets) as bonus targets worth the extra distance because of their high scientific value.

### 3.4 Known G-Type Main-Sequence Stars Within 10 Parsecs

The full inventory from RECONS and Gaia, from which the target list is drawn:

| # | Star | Sp. Type | Teff (K) | Dist (pc) | Mass (Msun) | Tier | Known Planets | Binary? |
|---|------|----------|----------|-----------|-------------|------|---------------|---------|
| 1 | Alpha Centauri A | G2V | 5790 | 1.34 | 1.10 | A/B | 1 unconfirmed | Yes (wide) |
| 2 | Tau Ceti | G8.5V | 5375 | 3.65 | 0.78 | C | 4 unconfirmed | No |
| 3 | Eta Cassiopeiae A | G0V | 6012 | 5.92 | 1.03 | B | 1 unconfirmed | Yes (wide, 70+ AU) |
| 4 | 36 Ophiuchi A | G7V | 5300 | 5.95 | 0.85 | C | None | Triple system |
| 5 | 82 G. Eridani | G6V | 5401 | 6.04 | 0.70 | C | 3 confirmed (d in HZ) | No |
| 6 | Delta Pavonis | G8IV-V | 5604 | 6.11 | 1.05 | B/C | None (deep limits) | No |
| 7 | Xi Bootis A | G8Ve | 5400 | 6.70 | 0.90 | C | None | Yes (close) |
| 8 | Mu Cassiopeiae A | G5Vp | 5308 | 7.55 | 0.74 | C | None | Yes (astrom. binary) |
| 9 | 107 Piscium | G5V | 5280 | 7.47 | 0.83 | C | None | No |
| 10 | Mu Herculis A | G5IV | 5560 | 8.40 | 1.10 | B | None | Yes (wide) |
| 11 | Beta CVn (Chara) | G0V | 5880 | 8.44 | 1.02 | B | None | No |
| 12 | 61 Virginis | G7V | 5531 | 8.53 | 0.94 | C | 3 confirmed | No |
| 13 | Zeta Tucanae | F9.5V | 5956 | 8.60 | 0.99 | B | None | No |
| 14 | Chi1 Orionis A | G0V | 5890 | 8.66 | 1.01 | B | None | Yes (wide) |
| 15 | Beta Comae Berenices | G0V | 5960 | 9.15 | 1.04 | B | None | No |
| 16 | Groombridge 1830 | G8Vp | 5200 | 9.16 | 0.66 | C | None | No (high proper motion, halo star) |
| 17 | Kappa1 Ceti | G5V | 5665 | 9.16 | 1.04 | B | None | No (young, active) |

Plus approximately 2-4 borderline G9V/K0V stars (e.g., Sigma Draconis at 5.77 pc, classified K0V by some sources). The exact count of 19-21 depends on where catalogs draw the G/K line.

### 3.5 Candidate Target List (Initial 10)

Selected from the 10 pc inventory above, prioritizing: (1) single stars or wide binaries, (2) low activity, (3) rich archival data. Two bonus targets beyond 10 pc included for their unique scientific value.

| Priority | Star | Tier | Dist (pc) | Teff (K) | M* (Msun) | [Fe/H] | Known Planets | Key Data Available |
|----------|------|------|-----------|----------|-----------|--------|---------------|-------------------|
| 1 | **Alpha Centauri A** | A/B | 1.34 | 5790 | 1.10 | +0.20 | 1 unconfirmed | HARPS, ESPRESSO, Gaia; TESS saturated |
| 2 | **Tau Ceti** (HD 10700) | C | 3.65 | 5375 | 0.78 | -0.50 | 4 unconfirmed (e,f,g,h) | 9000+ HARPS measurements, TESS multi-sector, Gaia |
| 3 | **82 G. Eridani** (HD 20794) | C | 6.04 | 5401 | 0.70 | -0.40 | 3 confirmed (d in HZ) | HARPS + ESPRESSO, TESS Sectors 3/4/30/31, Gaia |
| 4 | **Eta Cassiopeiae A** | B | 5.92 | 6012 | 1.03 | -0.30 | 1 unconfirmed (approximately 22 Mearth) | RV multi-campaign, TESS, Gaia; wide binary |
| 5 | **Delta Pavonis** | B/C | 6.11 | 5604 | 1.05 | +0.33 | None (RV limit K < 1 m/s) | CORALIE + HARPS since 1998, TESS Sectors 1-4+, Gaia |
| 6 | **61 Virginis** (HD 115617) | C | 8.53 | 5531 | 0.94 | -0.02 | 3 confirmed (b, c, d) | HARPS + Keck + AAT, TESS, Gaia |
| 7 | **Beta CVn / Chara** | B | 8.44 | 5880 | 1.02 | -0.21 | None | RV surveys, TESS, Gaia |
| 8 | **Zeta Tucanae** | B | 8.60 | 5956 | 0.99 | -0.15 | None | RV surveys, TESS, Gaia |
| 9* | **18 Scorpii** (HD 146233) | A | 14.13 | 5817 | 1.02 | +0.05 | 1 candidate (super-Earth) | HARPS, TESS, Gaia; nearest solar twin |
| 10* | **HD 134060** | B | 24.15 | 5890 | 1.07 | +0.10 | 2 confirmed (b: 3.3d, c: 1168d) | 8-year HARPS survey, TESS, Gaia |

*Targets 9-10 are beyond 10 pc but included for high scientific value.

**Excluded from the 10 pc pool and why:**
- **36 Ophiuchi A**: Triple star system; dynamical complexity makes HZ stability uncertain.
- **Xi Bootis A**: Close binary (approximately 4.9 AU separation); severe RV and astrometric contamination.
- **Mu Cassiopeiae A**: Astrometric binary companion; very metal-poor ([Fe/H] = -0.70).
- **Groombridge 1830**: Halo star with extreme kinematics; likely very old and metal-poor ([Fe/H] approximately -1.3); not representative of solar neighborhood.
- **Kappa1 Ceti**: Young (approximately 400 Myr) and magnetically active; high RV jitter.
- **107 Piscium, Mu Herculis A, Chi1 Orionis A, Beta Com**: Retained as alternates; can substitute if primary targets prove intractable.

### 3.6 Notes on Top Targets

- **Alpha Centauri A** is the nearest G2V star but the binary companion (alpha Cen B at 11-35 AU separation) complicates RV and astrometric analysis. It represents the hardest but highest-value target.
- **Tau Ceti** has the richest RV dataset of any single G star (9000+ HARPS measurements, approximately 30 cm/s RV dispersion). Its low metallicity ([Fe/H] = -0.50) may suppress rocky planet formation, but the massive debris disk suggests an active planetary system.
- **82 G. Eridani** is the strongest current result: a confirmed 5.8 Mearth super-Earth in the habitable zone (planet d, confirmed 2024-2025 via HARPS + ESPRESSO). This makes it both a validation target and a system where additional lower-mass planets may await detection.
- **18 Scorpii** at 14.1 pc is the nearest true "solar twin" (nearly identical Teff, log g, metallicity, and activity level to the Sun). If any star harbors a system like ours, this is a prime candidate.

---

## 4. Methodology Per Target

For each target star, the PoC will systematically inventory and apply available methods:

### 4.1 Method Matrix

| Method | Data Source | What It Detects | Sensitivity for G-star HZ | Applies to Targets |
|--------|------------|-----------------|---------------------------|-------------------|
| **Radial Velocity** | HARPS archive, ESPRESSO, Keck | Planet mass (M*sin(i)) | approximately 2-5 Mearth (current); approximately 1 Mearth (future ESPRESSO) | All 10 |
| **Transit (TESS)** | MAST / lightkurve | Planet radius | approximately 1-2 Rearth for short P; limited for HZ P | 8/10 (excluding alpha Cen A/B saturated) |
| **Transit (TESS FFI)** | MAST Full Frame Images | Same, with lower cadence | Lower sensitivity than 2-min cadence | Potentially all |
| **Astrometry (Gaia)** | Gaia DR3 (+ DR4 when available) | True planet mass, wide orbits | Super-Jupiter at 2-5 AU | All 10 |
| **Proper motion anomaly** | Gaia DR3 vs Hipparcos | Long-period massive companions | Jupiter-mass at 1-10 AU | All 10 |
| **Direct imaging (archival)** | VLT/SPHERE, Gemini/GPI archives | Wide-separation massive planets | Young Jupiters at > 10 AU | Select targets |
| **Literature synthesis** | NASA Exoplanet Archive, ADS | All known constraints | Comprehensive | All 10 |

### 4.2 Per-Target Analysis Workflow

For each of the 10 targets:

1. **Literature review**: What is known? What has been searched for and not found? What are the current detection limits?
2. **Existing data inventory**: Query MAST, ESO archive, Gaia archive. How many TESS sectors? How many RV measurements? What precision?
3. **Run existing pipeline**: Apply Modules 1-5 for stellar properties and transit search.
4. **RV analysis** (new capability): Download public RV time series. Compute periodogram. Assess detection limits via injection-recovery.
5. **Astrometric constraints**: Extract Gaia DR3 proper motion anomaly and RUWE. Compute sensitivity to companions.
6. **Sensitivity map**: For each method, compute the minimum detectable planet mass/radius as a function of orbital period. Overlay all methods to show the combined sensitivity.
7. **System architecture assessment**: Given known planets and non-detections, what can we say about the system's architecture? Where are the gaps?

### 4.3 New Code Required

| Component | Description | Builds on |
|-----------|-------------|-----------|
| RV data retrieval | Download public HARPS/ESPRESSO RV time series from ESO/DACE archives | data_access.py |
| RV periodogram | Lomb-Scargle on RV data; compute detection limits | periodogram.py |
| Gaia PMa analysis | Proper motion anomaly from Gaia-Hipparcos comparison | data_access.py (Gaia queries) |
| Sensitivity calculator | Injection-recovery for transit; analytical limits for RV/astrometry | New module |
| Per-target report generator | Automated summary of all methods per target | pipeline.py |

---

## 5. Success Criteria and Deliverables

### 5.1 What Constitutes Success for This PoC

This is a research phase. Success is not measured by planet detections but by the quality of the framework established:

| Deliverable | Success Criterion |
|------------|-------------------|
| **Target catalog** | Final prioritized list of approximately 10 stars with justified selection criteria. Criteria refined based on actual data availability. |
| **Per-target data inventory** | For each target: what data exists, its quality, and gaps. Documented systematically. |
| **Per-target sensitivity maps** | Combined detection limits across all methods. Clear visualization of "what we can and cannot detect" for each star. |
| **Method proficiency** | Demonstrated ability to retrieve and analyze RV, transit, and astrometric data. Documented gotchas and limitations. |
| **Gap analysis** | For each target: what planet parameter space remains unexplored? What new data would improve coverage? |
| **Architecture decisions** | Should the next phase automate this process? Scale to more targets? Focus on a specific method? Wait for specific data releases (Gaia DR4, PLATO)? |
| **Refined search criteria** | Based on the 10-target PoC, what Teff / distance / metallicity / variability thresholds actually matter? |

### 5.2 What a Non-Detection Tells Us

A non-detection is not a failure. For each target where no new planet is found, we produce:
- A quantified upper limit on planet mass/radius as a function of period
- Understanding of where the detection frontier lies for that specific star
- A clear statement of what observations would be needed to push deeper

This is standard practice in observational astronomy and has direct scientific value.

---

## 6. Implementation Roadmap

### Phase 7a: Foundation (approximately 2-3 sessions)

- [ ] Finalize target list (validate against Gaia DR3; confirm data availability)
- [ ] Build RV data retrieval module (DACE archive API or ESO archive download)
- [ ] Build Gaia proper motion anomaly analysis
- [ ] Write per-target report template

### Phase 7b: First Deep Dive -- 82 G. Eridani (approximately 2-3 sessions)

Start with the strongest target to develop methodology:
- [ ] Run full pipeline (Modules 1-5) on HD 20794
- [ ] Retrieve and analyze public HARPS + ESPRESSO RV time series
- [ ] Reproduce the confirmed planet d detection as validation
- [ ] Assess: can we detect additional lower-mass planets in the data?
- [ ] Compute combined transit + RV + astrometric sensitivity map
- [ ] Document the full workflow; identify pain points

**Why 82 G. Eridani first**: It has a confirmed HZ planet (d), rich archival data, and known RV detection limits. Starting here lets us validate our methods against a known result before applying them to targets with no confirmed planets.

### Phase 7c: Expand to Top 5 (approximately 3-5 sessions)

Apply the refined workflow to:
- [ ] Tau Ceti -- richest RV dataset; test planet candidate validation
- [ ] Alpha Centauri A -- hardest target; test binary decontamination
- [ ] Eta Cassiopeiae A -- test long-period RV signal assessment
- [ ] Delta Pavonis -- test deep non-detection characterization

### Phase 7d: Complete the 10 + Synthesis (approximately 2-3 sessions)

- [ ] Remaining 5 targets
- [ ] Synthesize: what have we learned about the search space?
- [ ] Refine hyperparameters and success criteria
- [ ] Decide: automate to larger sample? Wait for new data? Focus on specific method?

---

## 7. Technical Considerations

### 7.1 Data Access

| Archive | Access Method | Authentication | Notes |
|---------|-------------|----------------|-------|
| MAST (TESS/Kepler) | lightkurve (Python) | None | Already integrated in pipeline |
| Gaia DR3 | astroquery.gaia (TAP+) | None | Already integrated in pipeline |
| ESO Archive (HARPS, ESPRESSO) | HTTP download or astroquery.eso | ESO account (free registration) | New; need to assess API |
| DACE (RV database) | dace-query Python package or REST API | None for public data | Geneva Observatory RV repository |
| NASA Exoplanet Archive | TAP query or pyvo | None | Planet parameters and orbital solutions |

### 7.2 New Dependencies

| Package | Purpose | Security Notes |
|---------|---------|---------------|
| dace-query | Access DACE RV archive | Geneva Observatory package; open source |
| pyvo | Virtual Observatory access for catalog queries | Astropy-affiliated; well-maintained |

### 7.3 Existing Pipeline Reuse

The Phase 1-6 pipeline is directly reusable:
- Modules 1-3 for stellar characterization of targets
- Module 4 for light curve variability analysis
- Module 5 for transit detection (with RV-informed period constraints)
- Batch framework for running multiple targets

---

## 8. Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Public RV data is sparse for some targets | Cannot assess sensitivity | Focus on targets with rich archival data; document gaps |
| TESS saturates on very bright targets (V < 4) | Cannot do transit photometry for alpha Cen | Use FFI halo photometry or accept the limitation |
| Gaia DR3 astrometric solutions have limited planet sensitivity | Cannot constrain HZ planets | Acknowledge; plan for DR4 (Dec 2026) |
| PoC reveals the search is intractable with available data | Project stalls | This is itself a valid result -- documenting why is valuable |
| Scope creep into building a full multi-method pipeline | Over-engineering before understanding the problem | Keep Phase 7 manual and per-target; automation is Phase 8 |

---

## 9. References

- Anglada-Escude, G. et al. (2012). A planetary system around the nearby M dwarf GJ 667C. A&A 556, A126.
- Cretignier, M. et al. (2024). YARARA v2: improved RV extraction from HARPS and ESPRESSO. A&A.
- Eisner, N. et al. (2021). Planet Hunters TESS III: two transiting planets around the bright G dwarf HD 152843. MNRAS 501, 4669.
- Hall, R. D. et al. (2025). ESPRESSO observations of HD 10700. A&A.
- Kervella, P. et al. (2019). Stellar and substellar companions of nearby stars from Gaia DR2. A&A 623, A72.
- Pepe, F. et al. (2021). ESPRESSO at VLT: An Instrument for Exoplanet Research. A&A 645, A96.
- Rauer, H. et al. (2014). The PLATO 2.0 mission. Experimental Astronomy 38, 249.
- Thompson, S. J. et al. (2016). The Terra Hunting Experiment. Proc. SPIE 9908.
- Tuomi, M. et al. (2013). Signals embedded in the radial velocity noise. A&A 551, A79.
