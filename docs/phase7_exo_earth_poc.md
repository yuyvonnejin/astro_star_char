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

### 3.1 Selection Criteria

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Spectral type | G0V -- G9V (Teff 5200-6000 K) | Solar analog; habitable zone at approximately 0.7-1.5 AU |
| Luminosity class | V (main sequence) | Excludes subgiants/giants with different HZ geometry |
| Distance | < 10 pc preferred; up to approximately 25 pc for solar twins | Maximizes sensitivity across all methods |
| Stellar variability | Low (amplitude < 5 ppt if known) | Enables sensitive photometry and RV |
| Multiplicity | Singles or wide binaries (separation > 100 AU) preferred | Avoids dynamical disruption of HZ and RV contamination |
| Data availability | At least 2 of: TESS photometry, RV time series, Gaia astrometry | Enables multi-method analysis |

### 3.2 Candidate Target List (Initial 10)

Ranked by proximity and solar similarity. This list is a **starting point** to be refined during the PoC.

| Priority | Star | Sp. Type | Dist (pc) | M* (Msun) | [Fe/H] | Known Planets | Key Data Available |
|----------|------|----------|-----------|-----------|--------|---------------|-------------------|
| 1 | **Alpha Centauri A** | G2V | 1.33 | 1.10 | +0.20 | 1 unconfirmed candidate | HARPS, ESPRESSO, Gaia; TESS saturated |
| 2 | **Tau Ceti** (HD 10700) | G8.5V | 3.65 | 0.78 | -0.50 | 4 unconfirmed candidates (e,f,g,h) | 9000+ HARPS measurements, TESS multi-sector, Gaia |
| 3 | **82 G. Eridani** (HD 20794) | G6V | 6.04 | 0.70 | -0.40 | 3 confirmed (d in HZ, 5.8 Mearth) | HARPS + ESPRESSO, TESS Sectors 3/4/30/31, Gaia |
| 4 | **Eta Cassiopeiae A** | G0V | 5.92 | 1.03 | -0.30 | 1 unconfirmed (approximately 22 Mearth, P approximately 850d) | RV multi-campaign, TESS, Gaia; wide binary |
| 5 | **Delta Pavonis** | G8IV-V | 6.11 | 1.05 | +0.33 | None (RV limit K < 1 m/s) | CORALIE + HARPS since 1998, TESS Sectors 1-4+, Gaia |
| 6 | **61 Virginis** (HD 115617) | G7V | 8.53 | 0.94 | -0.02 | 3 confirmed (b, c, d) | HARPS + Keck + AAT, TESS, Gaia |
| 7 | **Beta CVn / Chara** | G0V | 8.44 | 1.02 | -0.21 | None | RV surveys, TESS, Gaia |
| 8 | **Zeta Tucanae** | F9.5V | 8.60 | 0.99 | -0.15 | None | RV surveys, TESS, Gaia |
| 9 | **18 Scorpii** (HD 146233) | G2Va | 14.13 | 1.02 | +0.05 | 1 candidate (super-Earth, P approximately 19.9d) | HARPS, TESS, Gaia; nearest solar twin |
| 10 | **HD 134060** | G0V | 24.15 | 1.07 | +0.10 | 2 confirmed (b: 3.3d, c: 1168d) | 8-year HARPS survey, TESS, Gaia |

### 3.3 Notes on Target Selection

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
