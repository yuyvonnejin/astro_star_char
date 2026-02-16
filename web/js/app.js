// Main application controller

let currentStarName = null;
let currentStarData = null;
let currentResult = null;

// -- Entry points ----------------------------------------------------------

function loadExample(name) {
  const star = EXAMPLE_STARS[name];
  if (!star) return;

  document.getElementById("star-input").value = name;
  currentStarName = name;
  currentStarData = Object.assign({}, star);
  runPipeline(name, currentStarData);
}

async function lookupStar() {
  const input = document.getElementById("star-input").value.trim();
  if (!input) return;

  // Check if it matches an example name (case-insensitive)
  for (const [name, data] of Object.entries(EXAMPLE_STARS)) {
    if (name.toLowerCase() === input.toLowerCase()) {
      loadExample(name);
      return;
    }
  }

  showStatus("Looking up star...", "loading");
  hideAllSteps();

  try {
    let sourceId;
    let starData;

    // Check if input looks like a Gaia source ID (all digits)
    if (/^\d+$/.test(input)) {
      sourceId = input;
    } else {
      // Try SIMBAD resolution
      showStatus("Resolving name via SIMBAD...", "loading");
      sourceId = await resolveSimbadName(input);
      if (!sourceId) {
        showStatus("Could not resolve star name via SIMBAD. Try a Gaia DR3 source ID or use an example star.", "error");
        return;
      }
    }

    showStatus(`Querying Gaia DR3 for source ${sourceId}...`, "loading");
    starData = await queryGaiaById(sourceId);
    if (!starData) {
      showStatus(`No Gaia DR3 data found for source_id ${sourceId}.`, "error");
      return;
    }

    currentStarName = input;
    currentStarData = starData;
    runPipeline(input, starData);
  } catch (e) {
    showStatus(`API error: ${e.message}. Try using an example star instead.`, "error");
  }
}

function runFromManualInput() {
  const fields = [
    "parallax_mas", "parallax_error_mas", "phot_g_mean_mag",
    "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "ag_gspphot", "ebpminrp_gspphot", "feh", "logg",
  ];

  const star = {
    source_id: "manual_input",
    is_cepheid: document.getElementById("manual-is-cepheid").checked,
    cepheid_period_days: null,
    teff_gspphot: null,
    lum_gspphot: null,
  };

  for (const f of fields) {
    const el = document.getElementById("manual-" + f);
    const val = el ? el.value.trim() : "";
    star[f] = val === "" ? null : parseFloat(val);
  }

  // Compute bp_rp if not provided but BP and RP are
  if (star.bp_rp == null && star.phot_bp_mean_mag != null && star.phot_rp_mean_mag != null) {
    star.bp_rp = star.phot_bp_mean_mag - star.phot_rp_mean_mag;
  }

  const periodEl = document.getElementById("manual-cepheid-period");
  if (periodEl && periodEl.value.trim()) {
    star.cepheid_period_days = parseFloat(periodEl.value);
  }

  const teffEl = document.getElementById("manual-teff-gspphot");
  if (teffEl && teffEl.value.trim()) star.teff_gspphot = parseFloat(teffEl.value);

  const lumEl = document.getElementById("manual-lum-gspphot");
  if (lumEl && lumEl.value.trim()) star.lum_gspphot = parseFloat(lumEl.value);

  if (star.parallax_mas == null || star.parallax_error_mas == null || star.phot_g_mean_mag == null || star.bp_rp == null) {
    showStatus("Parallax, parallax error, G magnitude, and BP-RP color are required.", "error");
    return;
  }

  currentStarName = "Manual Input";
  currentStarData = star;
  runPipeline("Manual Input", star);
}


// -- Pipeline execution ----------------------------------------------------

function runPipeline(name, star) {
  showStatus("Running pipeline...", "loading");
  hideAllSteps();

  // Apply defaults (matching pipeline.py)
  const s = Object.assign({}, star);
  if (s.ag_gspphot == null) s.ag_gspphot = 0.0;
  if (s.ebpminrp_gspphot == null) s.ebpminrp_gspphot = 0.0;
  if (s.feh == null) s.feh = 0.0;
  if (s.logg == null) s.logg = 4.0;

  const result = {source_id: s.source_id};

  // Show Step 0: input data
  renderStep0(name, s);

  // Module 1: Distance
  const distResult = computeDistance(s);
  Object.assign(result, distResult);
  renderStep1(distResult);

  if (distResult.distance_pc == null) {
    showStatus("Pipeline stopped: distance computation failed.", "error");
    currentResult = result;
    return;
  }

  // Module 2: Temperature / Luminosity / Radius
  const tempResult = computeTemperatureLuminosity(s, distResult.distance_pc);
  Object.assign(result, tempResult);
  renderStep2(tempResult, distResult.distance_pc);

  // Module 3: Mass
  const teff_K = tempResult.teff_K;
  const luminosity = tempResult.luminosity_Lsun;
  if (teff_K != null) {
    const massResult = computeMass(s, teff_K, luminosity);
    Object.assign(result, massResult);
    renderStep3(massResult, luminosity);
  } else {
    result.mass_Msun = null;
    result.mass_flag = "no_teff";
    result.is_main_sequence = null;
  }

  currentResult = result;

  // Results summary and HR diagram
  renderSummary(name, result);
  renderHRSection(name, result);

  showStatus("Pipeline complete.", "success");
}


// -- Rendering functions ---------------------------------------------------

function renderStep0(name, star) {
  const section = document.getElementById("step-0");
  section.style.display = "block";

  const fields = [
    ["source_id", "Source ID", star.source_id],
    ["parallax_mas", "Parallax (mas)", fmt(star.parallax_mas, 4)],
    ["parallax_error_mas", "Parallax Error (mas)", fmt(star.parallax_error_mas, 4)],
    ["phot_g_mean_mag", "G magnitude", fmt(star.phot_g_mean_mag, 3)],
    ["phot_bp_mean_mag", "BP magnitude", fmt(star.phot_bp_mean_mag, 3)],
    ["phot_rp_mean_mag", "RP magnitude", fmt(star.phot_rp_mean_mag, 3)],
    ["bp_rp", "BP - RP color", fmt(star.bp_rp, 4)],
    ["ag_gspphot", "G-band extinction A_G", fmt(star.ag_gspphot, 4)],
    ["ebpminrp_gspphot", "E(BP-RP) reddening", fmt(star.ebpminrp_gspphot, 4)],
    ["feh", "[Fe/H] metallicity", fmt(star.feh, 3)],
    ["logg", "log(g) surface gravity", fmt(star.logg, 3)],
    ["is_cepheid", "Cepheid variable?", star.is_cepheid ? "Yes" : "No"],
  ];
  if (star.is_cepheid) {
    fields.push(["cepheid_period_days", "Cepheid Period (days)", fmt(star.cepheid_period_days, 3)]);
  }
  if (star.teff_gspphot != null) {
    fields.push(["teff_gspphot", "Gaia Teff estimate (K)", fmt(star.teff_gspphot, 0)]);
  }
  if (star.lum_gspphot != null) {
    fields.push(["lum_gspphot", "Gaia Luminosity estimate (Lsun)", fmt(star.lum_gspphot, 4)]);
  }

  let html = '<table class="data-table"><thead><tr><th>Parameter</th><th>Value</th><th>Used in</th></tr></thead><tbody>';
  const usage = {
    parallax_mas: "Module 1",
    parallax_error_mas: "Module 1",
    phot_g_mean_mag: "Module 1, 2",
    bp_rp: "Module 2",
    ag_gspphot: "Module 1, 2",
    ebpminrp_gspphot: "Module 2",
    feh: "Module 2",
    logg: "Module 2, 3",
    is_cepheid: "Module 1",
    cepheid_period_days: "Module 1",
    lum_gspphot: "Validation",
    teff_gspphot: "Comparison",
  };

  for (const [key, label, val] of fields) {
    html += `<tr><td>${label}</td><td><code>${val}</code></td><td>${usage[key] || ""}</td></tr>`;
  }
  html += "</tbody></table>";

  document.getElementById("step-0-content").innerHTML = html;
}


function renderStep1(distResult) {
  const section = document.getElementById("step-1");
  section.style.display = "block";

  const inter = distResult._intermediates || {};
  let html = "";

  if (distResult.distance_method === "cepheid_leavitt") {
    html += '<h4>Method: Cepheid Period-Luminosity (Leavitt Law)</h4>';
    html += '<div class="formula">M_V = -2.43 * log10(P) - 4.05</div>';
    html += '<div class="formula">d = 10^((G_0 - M_V + 5) / 5)</div>';
    html += '<div class="calc-steps">';
    html += `<p>Period P = ${fmt(inter.period_days, 3)} days</p>`;
    html += `<p>M_V = -2.43 * log10(${fmt(inter.period_days, 3)}) - 4.05 = <strong>${fmt(inter.M_V, 3)}</strong></p>`;
    html += `<p>G_0 = G - A_G = ${fmt(inter.G_0, 3)}</p>`;
    html += `<p>Distance = 10^((${fmt(inter.G_0, 3)} - (${fmt(inter.M_V, 3)}) + 5) / 5) = <strong>${fmt(distResult.distance_pc, 2)} pc</strong></p>`;
    html += '</div>';
    html += '<div class="note">Note: Uses G magnitude as proxy for V (introduces ~1% systematic error).</div>';
  } else {
    html += '<h4>Method: Bayesian Parallax Inversion</h4>';
    html += '<div class="formula">Prior: P(r) = r^2 * exp(-r / L),  L = 1350 pc</div>';
    html += '<div class="formula">Likelihood: P(plx_obs | r) = Normal(1000/r, sigma_plx)</div>';
    html += '<div class="formula">Posterior = Prior * Likelihood  (normalized)</div>';
    html += '<div class="calc-steps">';
    html += `<p>Parallax (corrected for zero-point +0.017 mas): ${fmt(inter.parallax_corrected, 4)} mas</p>`;
    html += `<p>Naive distance (1000 / parallax): ${fmt(inter.naive_distance, 4)} pc</p>`;
    html += `<p>Bayesian MAP estimate: <strong>${fmt(distResult.distance_pc, 4)} pc</strong></p>`;
    if (distResult.distance_lower_pc != null) {
      html += `<p>68% credible interval: [${fmt(distResult.distance_lower_pc, 4)}, ${fmt(distResult.distance_upper_pc, 4)}] pc</p>`;
    }
    const diff = inter.naive_distance ? Math.abs(distResult.distance_pc - inter.naive_distance) : 0;
    const pctDiff = inter.naive_distance ? (diff / inter.naive_distance * 100) : 0;
    html += `<p>Bayesian vs Naive difference: ${fmt(diff, 4)} pc (${fmt(pctDiff, 2)}%)</p>`;
    html += '</div>';

    // Chart
    html += '<div id="posterior-chart"></div>';
  }

  document.getElementById("step-1-content").innerHTML = html;

  // Render chart after DOM update
  if (distResult.distance_method !== "cepheid_leavitt" && inter.r_grid) {
    setTimeout(() => {
      renderPosteriorChart("posterior-chart", inter, distResult.distance_pc, distResult.distance_lower_pc, distResult.distance_upper_pc);
    }, 50);
  }
}


function renderStep2(tempResult, distance_pc) {
  const section = document.getElementById("step-2");
  section.style.display = "block";

  const inter = tempResult._intermediates || {};
  let html = "";

  // Step 2a: Dereddening
  html += '<h4>Step 2a: Dereddening</h4>';
  html += '<div class="formula">(BP-RP)_0 = BP-RP - E(BP-RP)</div>';
  html += '<div class="formula">G_0 = G - A_G</div>';
  html += '<div class="calc-steps">';
  html += `<p>(BP-RP)_0 = ${fmt(tempResult.bp_rp_0, 4)}</p>`;
  html += `<p>G_0 = ${fmt(inter.g_0, 4)}</p>`;
  html += '</div>';

  // Step 2b: Temperature
  html += '<h4>Step 2b: Effective Temperature (Mucciarelli et al. 2021)</h4>';
  html += `<div class="note">Using ${inter.coeffs_used} coefficients (logg ${inter.is_dwarf ? ">=" : "<"} 3.0)</div>`;
  html += '<div class="formula">theta = b0 + b1*C + b2*C^2 + b3*[Fe/H] + b4*[Fe/H]^2 + b5*[Fe/H]*C</div>';
  html += '<div class="formula">Teff = 5040 / theta</div>';
  html += '<div class="calc-steps">';
  html += `<p>C = (BP-RP)_0 = ${fmt(tempResult.bp_rp_0, 4)}, [Fe/H] = ${fmt(inter.feh, 3)}</p>`;
  html += `<p>theta = ${fmt(inter.theta, 6)}</p>`;
  html += `<p>Teff = 5040 / ${fmt(inter.theta, 6)} = <strong>${fmt(tempResult.teff_K, 1)} +/- ${tempResult.teff_uncertainty_K} K</strong></p>`;
  if (tempResult.teff_flag !== "ok") {
    html += `<p class="warning">Warning: Color is outside valid calibration range (flag: ${tempResult.teff_flag})</p>`;
  }
  html += '</div>';

  // Step 2c: Absolute magnitude
  html += '<h4>Step 2c: Absolute Magnitude</h4>';
  html += '<div class="formula">M_G = G_0 - 5*log10(d) + 5</div>';
  html += '<div class="calc-steps">';
  html += `<p>M_G = ${fmt(inter.g_0, 4)} - 5*log10(${fmt(distance_pc, 4)}) + 5 = <strong>${fmt(tempResult.M_G, 4)}</strong></p>`;
  html += '</div>';

  if (tempResult.BC_G == null) {
    html += '<div class="warning">Bolometric correction unavailable: Teff outside valid range [3300, 8000] K. Luminosity/radius/mass cannot be computed.</div>';
    document.getElementById("step-2-content").innerHTML = html;
    return;
  }

  // Step 2d: Bolometric correction
  html += '<h4>Step 2d: Bolometric Correction (Andrae et al. 2018)</h4>';
  html += `<div class="note">Regime: ${inter.bc_g_regime}</div>`;
  html += '<div class="formula">BC_G = a0 + a1*dT + a2*dT^2 + a3*dT^3 + a4*dT^4, where dT = Teff - 5772</div>';
  html += '<div class="calc-steps">';
  html += `<p>dT = ${fmt(inter.dt, 1)}</p>`;
  html += `<p>BC_G = <strong>${fmt(tempResult.BC_G, 5)}</strong></p>`;
  html += '</div>';

  // Step 2e: Luminosity
  html += '<h4>Step 2e: Luminosity</h4>';
  html += '<div class="formula">M_bol = M_G + BC_G</div>';
  html += '<div class="formula">L/Lsun = 10^((M_bol_sun - M_bol) / 2.5),  M_bol_sun = 4.74</div>';
  html += '<div class="calc-steps">';
  html += `<p>M_bol = ${fmt(tempResult.M_G, 4)} + (${fmt(tempResult.BC_G, 5)}) = <strong>${fmt(tempResult.M_bol, 4)}</strong></p>`;
  html += `<p>L = 10^((4.74 - ${fmt(tempResult.M_bol, 4)}) / 2.5) = <strong>${fmt(tempResult.luminosity_Lsun, 5)} Lsun</strong></p>`;
  html += '</div>';

  // Step 2f: Radius
  html += '<h4>Step 2f: Radius (Stefan-Boltzmann)</h4>';
  html += '<div class="formula">R/Rsun = sqrt(L/Lsun) * (Tsun/Teff)^2,  Tsun = 5772 K</div>';
  html += '<div class="calc-steps">';
  html += `<p>R = sqrt(${fmt(tempResult.luminosity_Lsun, 5)}) * (5772 / ${fmt(tempResult.teff_K, 1)})^2 = <strong>${fmt(tempResult.radius_Rsun, 4)} Rsun</strong></p>`;
  html += '</div>';

  // Step 2g: Validation
  if (tempResult.luminosity_validation_ratio != null) {
    html += '<h4>Step 2g: Validation vs Gaia Estimate</h4>';
    html += '<div class="calc-steps">';
    html += `<p>Computed L / Gaia L = ${fmt(tempResult.luminosity_validation_ratio, 4)}</p>`;
    const good = tempResult.luminosity_validation_ratio >= 0.8 && tempResult.luminosity_validation_ratio <= 1.2;
    html += `<p class="${good ? 'good' : 'warning'}">${good ? "Good agreement (within 20%)" : "Significant discrepancy (>20%)"}</p>`;
    html += '</div>';
  }

  document.getElementById("step-2-content").innerHTML = html;
}


function renderStep3(massResult, luminosity) {
  const section = document.getElementById("step-3");
  section.style.display = "block";

  const inter = massResult._intermediates || {};
  let html = "";

  html += '<h4>Main-Sequence Check</h4>';
  if (!massResult.is_main_sequence) {
    html += `<div class="warning">Not classified as main-sequence: ${inter.reason}</div>`;
    html += '<div class="note">Mass estimation via the mass-luminosity relation is only valid for main-sequence stars.</div>';
    document.getElementById("step-3-content").innerHTML = html;
    return;
  }

  if (massResult.mass_flag === "no_luminosity") {
    html += '<div class="warning">Luminosity not available; cannot compute mass.</div>';
    document.getElementById("step-3-content").innerHTML = html;
    return;
  }

  html += '<div class="calc-steps"><p>logg >= 3.5 and 3300 K <= Teff <= 8000 K: <strong>Main-sequence star confirmed</strong></p></div>';

  html += '<h4>Piecewise Mass-Luminosity Relation</h4>';
  html += '<div class="formula">L = M^alpha  =>  M = L^(1/alpha)</div>';
  html += '<table class="data-table"><thead><tr><th>L range (Lsun)</th><th>alpha</th><th>M range (Msun)</th></tr></thead><tbody>';
  html += '<tr><td>L < 0.033</td><td>2.3</td><td>M < 0.43</td></tr>';
  html += '<tr><td>0.033 - 16</td><td>4.0</td><td>0.43 - 2</td></tr>';
  html += '<tr><td>16 - 54,000</td><td>3.5</td><td>2 - 55</td></tr>';
  html += '<tr><td>> 54,000</td><td>1.0</td><td>> 55</td></tr>';
  html += '</tbody></table>';

  html += '<div class="calc-steps">';
  html += `<p>L = ${fmt(inter.luminosity_Lsun, 5)} Lsun</p>`;
  html += `<p>Regime: ${inter.regime}</p>`;
  html += `<p>alpha = ${inter.alpha}</p>`;
  html += `<p>M = ${fmt(inter.luminosity_Lsun, 5)}^(1/${inter.alpha}) = <strong>${fmt(massResult.mass_Msun, 4)} Msun</strong></p>`;
  html += '</div>';

  html += '<div id="mass-lum-chart"></div>';

  document.getElementById("step-3-content").innerHTML = html;

  setTimeout(() => {
    renderMassLuminosityChart("mass-lum-chart", luminosity, massResult.mass_Msun);
  }, 50);
}


function renderSummary(name, result) {
  const section = document.getElementById("summary");
  section.style.display = "block";

  let html = `<h3>Results for: ${escapeHtml(name)}</h3>`;

  // Main results table
  html += '<table class="data-table results-table"><thead><tr><th>Property</th><th>Pipeline Result</th>';

  const lit = LITERATURE_VALUES[name];
  if (lit) html += '<th>Literature Value</th><th>Difference</th>';
  html += '</tr></thead><tbody>';

  const rows = [
    ["Distance", fmt(result.distance_pc, 4) + " pc", null, null],
    ["Teff", result.teff_K != null ? fmt(result.teff_K, 0) + " K" : "N/A", lit ? lit.teff_K + " K" : null, result.teff_K && lit ? pctDiff(result.teff_K, lit.teff_K) : null],
    ["Luminosity", result.luminosity_Lsun != null ? fmt(result.luminosity_Lsun, 5) + " Lsun" : "N/A", lit && lit.luminosity_Lsun != null ? fmt(lit.luminosity_Lsun, 5) + " Lsun" : null, result.luminosity_Lsun && lit && lit.luminosity_Lsun ? pctDiff(result.luminosity_Lsun, lit.luminosity_Lsun) : null],
    ["Radius", result.radius_Rsun != null ? fmt(result.radius_Rsun, 4) + " Rsun" : "N/A", lit && lit.radius_Rsun != null ? fmt(lit.radius_Rsun, 4) + " Rsun" : null, result.radius_Rsun && lit && lit.radius_Rsun ? pctDiff(result.radius_Rsun, lit.radius_Rsun) : null],
    ["Mass", result.mass_Msun != null ? fmt(result.mass_Msun, 4) + " Msun" : "N/A (" + (result.mass_flag || "?") + ")", lit && lit.mass_Msun != null ? fmt(lit.mass_Msun, 4) + " Msun" : null, result.mass_Msun && lit && lit.mass_Msun ? pctDiff(result.mass_Msun, lit.mass_Msun) : null],
  ];

  for (const [prop, pipeline, litVal, diff] of rows) {
    html += `<tr><td>${prop}</td><td><strong>${pipeline}</strong></td>`;
    if (lit) {
      html += `<td>${litVal || "--"}</td><td>${diff || "--"}</td>`;
    }
    html += '</tr>';
  }
  html += '</tbody></table>';

  if (lit) {
    html += `<div class="note">Literature source: ${lit.source}</div>`;
  }

  document.getElementById("summary-content").innerHTML = html;
}


function renderHRSection(name, result) {
  const section = document.getElementById("hr-diagram");
  section.style.display = "block";

  setTimeout(() => {
    renderHRDiagram("hr-chart", result.teff_K, result.luminosity_Lsun, name);
  }, 100);
}


// -- UI helpers ------------------------------------------------------------

function showStatus(msg, type) {
  const el = document.getElementById("status");
  el.textContent = msg;
  el.className = "status " + type;
  el.style.display = "block";
}

function hideAllSteps() {
  for (const id of ["step-0", "step-1", "step-2", "step-3", "summary", "hr-diagram"]) {
    document.getElementById(id).style.display = "none";
  }
}

function toggleManualInput() {
  const panel = document.getElementById("manual-panel");
  panel.style.display = panel.style.display === "none" ? "block" : "none";
}

function toggleStep(stepId) {
  const content = document.getElementById(stepId + "-content");
  const header = content.previousElementSibling;
  content.classList.toggle("collapsed");
  header.classList.toggle("collapsed");
}

function fmt(val, decimals) {
  if (val == null) return "N/A";
  return Number(val).toFixed(decimals);
}

function pctDiff(computed, reference) {
  const diff = Math.abs(computed - reference) / Math.abs(reference) * 100;
  return fmt(diff, 1) + "%";
}

function escapeHtml(str) {
  const div = document.createElement("div");
  div.textContent = str;
  return div.innerHTML;
}

// -- Init ------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", () => {
  document.getElementById("star-input").addEventListener("keydown", (e) => {
    if (e.key === "Enter") lookupStar();
  });
});
