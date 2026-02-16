// Chart rendering with Plotly.js

const DARK_LAYOUT = {
  paper_bgcolor: "rgba(0,0,0,0)",
  plot_bgcolor: "#0f172a",
  font: {color: "#e2e8f0"},
  xaxis: {gridcolor: "#334155", zerolinecolor: "#475569"},
  yaxis: {gridcolor: "#334155", zerolinecolor: "#475569"},
};

function renderPosteriorChart(containerId, intermediates, distance_pc, lower, upper) {
  if (!intermediates.r_grid) return;

  // Downsample for plotting (every 10th point)
  const step = 10;
  const x = [];
  const y = [];
  for (let i = 0; i < intermediates.r_grid.length; i += step) {
    x.push(intermediates.r_grid[i]);
    y.push(intermediates.posterior[i]);
  }

  // Focus plot around the peak: +/- 5x the distance
  const x_lo = Math.max(0.1, distance_pc * 0.5);
  const x_hi = distance_pc * 2.0;

  const traces = [{
    x: x,
    y: y,
    type: "scatter",
    mode: "lines",
    name: "Posterior P(d | parallax)",
    line: {color: "#2563eb", width: 2},
  }];

  // MAP line
  traces.push({
    x: [distance_pc, distance_pc],
    y: [0, Math.max(...y) * 1.05],
    type: "scatter",
    mode: "lines",
    name: `MAP = ${distance_pc.toFixed(3)} pc`,
    line: {color: "#dc2626", width: 2, dash: "dash"},
  });

  // Credible interval shading
  if (lower != null && upper != null) {
    const shade_x = [];
    const shade_y = [];
    for (let i = 0; i < intermediates.r_grid.length; i += step) {
      if (intermediates.r_grid[i] >= lower && intermediates.r_grid[i] <= upper) {
        shade_x.push(intermediates.r_grid[i]);
        shade_y.push(intermediates.posterior[i]);
      }
    }
    traces.push({
      x: shade_x,
      y: shade_y,
      type: "scatter",
      fill: "tozeroy",
      mode: "none",
      name: `68% CI [${lower.toFixed(3)}, ${upper.toFixed(3)}]`,
      fillcolor: "rgba(37, 99, 235, 0.2)",
    });
  }

  // Naive estimate line
  if (intermediates.naive_distance) {
    traces.push({
      x: [intermediates.naive_distance, intermediates.naive_distance],
      y: [0, Math.max(...y) * 1.05],
      type: "scatter",
      mode: "lines",
      name: `Naive 1/plx = ${intermediates.naive_distance.toFixed(3)} pc`,
      line: {color: "#9ca3af", width: 1.5, dash: "dot"},
    });
  }

  const layout = {
    ...DARK_LAYOUT,
    title: "Bayesian Distance Posterior Distribution",
    xaxis: {...DARK_LAYOUT.xaxis, title: "Distance (pc)", range: [x_lo, x_hi], type: "log"},
    yaxis: {...DARK_LAYOUT.yaxis, title: "Probability Density", showticklabels: false},
    showlegend: true,
    legend: {x: 0.6, y: 0.95, bgcolor: "rgba(0,0,0,0)"},
    margin: {t: 40, b: 50, l: 50, r: 20},
    height: 350,
  };

  Plotly.newPlot(containerId, traces, layout, {responsive: true});
}


function renderMassLuminosityChart(containerId, luminosity, mass) {
  // Piecewise mass-luminosity relation: L = M^alpha
  // Generate reference curves for each regime
  const regimes = [
    {alpha: 2.3, L_min: 0.001, L_max: 0.033, label: "alpha=2.3 (M dwarfs)", color: "#ef4444"},
    {alpha: 4.0, L_min: 0.033, L_max: 16, label: "alpha=4.0 (K/G/F dwarfs)", color: "#f59e0b"},
    {alpha: 3.5, L_min: 16, L_max: 54000, label: "alpha=3.5 (A+ stars)", color: "#10b981"},
    {alpha: 1.0, L_min: 54000, L_max: 500000, label: "alpha=1.0 (massive)", color: "#6366f1"},
  ];

  const traces = [];
  for (const r of regimes) {
    const masses = [];
    const lums = [];
    const n = 50;
    const logLmin = Math.log10(r.L_min);
    const logLmax = Math.log10(r.L_max);
    for (let i = 0; i <= n; i++) {
      const L = Math.pow(10, logLmin + (logLmax - logLmin) * i / n);
      const M = Math.pow(L, 1.0 / r.alpha);
      masses.push(M);
      lums.push(L);
    }
    traces.push({
      x: masses,
      y: lums,
      type: "scatter",
      mode: "lines",
      name: r.label,
      line: {color: r.color, width: 2},
    });
  }

  // Plot the star
  if (luminosity != null && mass != null) {
    traces.push({
      x: [mass],
      y: [luminosity],
      type: "scatter",
      mode: "markers",
      name: `This star (M=${mass.toFixed(3)})`,
      marker: {color: "#dc2626", size: 12, symbol: "star"},
    });
  }

  const layout = {
    ...DARK_LAYOUT,
    title: "Mass-Luminosity Relation",
    xaxis: {...DARK_LAYOUT.xaxis, title: "Mass (Msun)", type: "log"},
    yaxis: {...DARK_LAYOUT.yaxis, title: "Luminosity (Lsun)", type: "log"},
    showlegend: true,
    legend: {x: 0.02, y: 0.98, bgcolor: "rgba(0,0,0,0)"},
    margin: {t: 40, b: 50, l: 60, r: 20},
    height: 350,
  };

  Plotly.newPlot(containerId, traces, layout, {responsive: true});
}


function renderHRDiagram(containerId, teff_K, luminosity_Lsun, starName) {
  // Main sequence reference track (approximate)
  const msTrack = [
    {teff: 30000, L: 100000},
    {teff: 20000, L: 10000},
    {teff: 10000, L: 50},
    {teff: 7500, L: 10},
    {teff: 6000, L: 2},
    {teff: 5772, L: 1},
    {teff: 5000, L: 0.4},
    {teff: 4000, L: 0.08},
    {teff: 3500, L: 0.02},
    {teff: 3000, L: 0.003},
    {teff: 2500, L: 0.0005},
  ];

  const traces = [];

  // Main sequence
  traces.push({
    x: msTrack.map(p => p.teff),
    y: msTrack.map(p => p.L),
    type: "scatter",
    mode: "lines",
    name: "Main Sequence (approx.)",
    line: {color: "#d1d5db", width: 3},
  });

  // Plot all example stars
  const exampleColors = {
    "Sun": "#f59e0b",
    "Proxima Cen": "#ef4444",
    "Sirius A": "#3b82f6",
    "Alpha Cen A": "#f97316",
    "Barnard's Star": "#ef4444",
    "Delta Cephei": "#8b5cf6",
  };

  for (const [name, lit] of Object.entries(LITERATURE_VALUES)) {
    if (lit.teff_K && lit.luminosity_Lsun) {
      traces.push({
        x: [lit.teff_K],
        y: [lit.luminosity_Lsun],
        type: "scatter",
        mode: "markers+text",
        name: name,
        text: [name],
        textposition: "top right",
        textfont: {size: 10, color: "#94a3b8"},
        marker: {
          color: exampleColors[name] || "#9ca3af",
          size: 8,
          symbol: "circle",
        },
        showlegend: false,
      });
    }
  }

  // Plot the current star prominently
  if (teff_K != null && luminosity_Lsun != null) {
    traces.push({
      x: [teff_K],
      y: [luminosity_Lsun],
      type: "scatter",
      mode: "markers",
      name: starName || "This star",
      marker: {color: "#dc2626", size: 14, symbol: "star", line: {width: 2, color: "#000"}},
    });
  }

  const layout = {
    ...DARK_LAYOUT,
    title: "HR Diagram",
    xaxis: {...DARK_LAYOUT.xaxis, title: "Effective Temperature (K)", autorange: "reversed", type: "log"},
    yaxis: {...DARK_LAYOUT.yaxis, title: "Luminosity (Lsun)", type: "log"},
    showlegend: true,
    legend: {x: 0.02, y: 0.02, yanchor: "bottom", bgcolor: "rgba(0,0,0,0)"},
    margin: {t: 40, b: 50, l: 60, r: 20},
    height: 400,
  };

  Plotly.newPlot(containerId, traces, layout, {responsive: true});
}
