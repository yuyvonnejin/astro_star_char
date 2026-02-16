// Quick validation: run JS pipeline on Sun and compare with Python output.
// Run with: node test_pipeline.js

const fs = require("fs");
const vm = require("vm");

// Create a shared context for all scripts
const context = vm.createContext({
  console: console,
  Math: Math,
  Number: Number,
  Array: Array,
  Object: Object,
  String: String,
  isNaN: isNaN,
  parseFloat: parseFloat,
  parseInt: parseInt,
});

// Load modules into shared context
const files = [
  "js/examples.js",
  "js/distance.js",
  "js/temperature.js",
  "js/mass.js",
];
for (const f of files) {
  vm.runInContext(fs.readFileSync(f, "utf8"), context, {filename: f});
}

// Run test in the same context
const testCode = `
const PYTHON_RESULTS = {
  sun: {
    distance_pc: 10.002303080542816,
    teff_K: 5683.94376558658,
    M_G: 4.829499949995,
    BC_G: 0.0535375797363369,
    luminosity_Lsun: 0.8765667402514375,
    radius_Rsun: 0.9654851398535276,
    mass_Msun: 0.9676008629800456,
  },
};

function testStar(name, starData, expected) {
  const s = Object.assign({}, starData);
  if (s.ag_gspphot == null) s.ag_gspphot = 0.0;
  if (s.ebpminrp_gspphot == null) s.ebpminrp_gspphot = 0.0;
  if (s.feh == null) s.feh = 0.0;
  if (s.logg == null) s.logg = 4.0;

  const distResult = computeDistance(s);
  const tempResult = computeTemperatureLuminosity(s, distResult.distance_pc);
  const massResult = computeMass(s, tempResult.teff_K, tempResult.luminosity_Lsun);

  const checks = [
    ["distance_pc", distResult.distance_pc, expected.distance_pc],
    ["teff_K", tempResult.teff_K, expected.teff_K],
    ["M_G", tempResult.M_G, expected.M_G],
    ["BC_G", tempResult.BC_G, expected.BC_G],
    ["luminosity_Lsun", tempResult.luminosity_Lsun, expected.luminosity_Lsun],
    ["radius_Rsun", tempResult.radius_Rsun, expected.radius_Rsun],
    ["mass_Msun", massResult.mass_Msun, expected.mass_Msun],
  ];

  console.log("\\n=== " + name + " ===");
  let allPass = true;
  for (const [field, jsVal, pyVal] of checks) {
    const diff = Math.abs(jsVal - pyVal);
    const relDiff = pyVal !== 0 ? (diff / Math.abs(pyVal) * 100) : diff;
    const pass = relDiff < 0.1; // 0.1% tolerance (Bayesian grid may differ slightly)
    const status = pass ? "PASS" : "FAIL";
    if (!pass) allPass = false;
    console.log("  " + status + "  " + field + ": JS=" + jsVal.toFixed(8) + ", PY=" + pyVal.toFixed(8) + ", diff=" + relDiff.toFixed(6) + "%");
  }
  return allPass;
}

let ok = true;
ok = testStar("Sun", EXAMPLE_STARS["Sun"], PYTHON_RESULTS.sun) && ok;

// Also test all example stars run without errors
for (const [name, data] of Object.entries(EXAMPLE_STARS)) {
  if (name === "Sun") continue;
  const s = Object.assign({}, data);
  if (s.ag_gspphot == null) s.ag_gspphot = 0.0;
  if (s.ebpminrp_gspphot == null) s.ebpminrp_gspphot = 0.0;
  if (s.feh == null) s.feh = 0.0;
  if (s.logg == null) s.logg = 4.0;

  try {
    const dr = computeDistance(s);
    const tr = computeTemperatureLuminosity(s, dr.distance_pc);
    if (tr.teff_K != null) {
      computeMass(s, tr.teff_K, tr.luminosity_Lsun);
    }
    console.log("  PASS  " + name + " completed without errors (d=" + (dr.distance_pc ? dr.distance_pc.toFixed(2) : "null") + " pc)");
  } catch (e) {
    console.log("  FAIL  " + name + ": " + e.message);
    ok = false;
  }
}

console.log(ok ? "\\nAll tests passed." : "\\nSome tests FAILED.");
ok;
`;

const result = vm.runInContext(testCode, context);
process.exit(result ? 0 : 1);
