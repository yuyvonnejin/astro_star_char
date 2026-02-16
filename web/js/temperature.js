// Module 2: Effective temperature, bolometric correction, luminosity, and radius.
// Ported from src/temperature.py

// Mucciarelli, Bellazzini & Massari (2021) -- Table 1
const DWARF_COEFFS = {
  b: [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395],
  color_range: [0.39, 1.50],
  dispersion_K: 61,
};
const GIANT_COEFFS = {
  b: [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009],
  color_range: [0.33, 1.81],
  dispersion_K: 83,
};

// Andrae et al. (2018) -- BC_G polynomial coefficients
const BC_COEFFS_HOT = {
  a: [6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15],
  teff_range: [4000, 8000],
};
const BC_COEFFS_COOL = {
  a: [1.749e+00, 1.977e-03, 3.737e-07, -8.966e-11, -4.183e-14],
  teff_range: [3300, 4000],
};

const M_BOL_SUN = 4.74;   // IAU 2015 B2
const T_EFF_SUN = 5772.0; // K

function computeTemperatureLuminosity(star, distance_pc) {
  const result = {};

  // Step 1: Deredden
  const ebpminrp = star.ebpminrp_gspphot || 0.0;
  const ag = star.ag_gspphot || 0.0;
  const bp_rp_0 = star.bp_rp - ebpminrp;
  const g_0 = star.phot_g_mean_mag - ag;
  result.bp_rp_0 = bp_rp_0;

  // Step 2: Effective temperature
  const logg = star.logg != null ? star.logg : 4.0;
  const feh = star.feh != null ? star.feh : 0.0;
  const is_dwarf = logg >= 3.0;

  const coeffs = is_dwarf ? DWARF_COEFFS : GIANT_COEFFS;
  const b = coeffs.b;
  const color_lo = coeffs.color_range[0];
  const color_hi = coeffs.color_range[1];

  let teff_flag = "ok";
  if (bp_rp_0 < color_lo || bp_rp_0 > color_hi) {
    teff_flag = "outside_valid_range";
  }

  const c = bp_rp_0;
  const theta = b[0] + b[1] * c + b[2] * c * c + b[3] * feh + b[4] * feh * feh + b[5] * feh * c;
  const teff_K = 5040.0 / theta;

  result.teff_K = teff_K;
  result.teff_uncertainty_K = coeffs.dispersion_K;
  result.teff_flag = teff_flag;

  // Step 3: Absolute G magnitude
  const m_g = g_0 - 5.0 * Math.log10(distance_pc) + 5.0;
  result.M_G = m_g;

  // Step 4: Bolometric correction BC_G
  const bc_g = computeBCG(teff_K);
  if (bc_g === null) {
    result.BC_G = null;
    result.M_bol = null;
    result.luminosity_Lsun = null;
    result.radius_Rsun = null;
    result.luminosity_validation_ratio = null;
    result._intermediates = {
      g_0: g_0,
      bp_rp_0: bp_rp_0,
      theta: theta,
      is_dwarf: is_dwarf,
      coeffs_used: is_dwarf ? "dwarf" : "giant",
      bc_g_failed: true,
    };
    return result;
  }
  result.BC_G = bc_g;

  // Step 5: Luminosity
  const m_bol = m_g + bc_g;
  const luminosity = Math.pow(10.0, (M_BOL_SUN - m_bol) / 2.5);
  result.M_bol = m_bol;
  result.luminosity_Lsun = luminosity;

  // Step 6: Radius via Stefan-Boltzmann
  const radius_Rsun = Math.sqrt(luminosity) * Math.pow(T_EFF_SUN / teff_K, 2);
  result.radius_Rsun = radius_Rsun;

  // Step 7: Validation cross-check
  const lum_gspphot = star.lum_gspphot;
  if (lum_gspphot != null && lum_gspphot > 0) {
    result.luminosity_validation_ratio = luminosity / lum_gspphot;
  } else {
    result.luminosity_validation_ratio = null;
  }

  // Store intermediates for display
  result._intermediates = {
    g_0: g_0,
    bp_rp_0: bp_rp_0,
    theta: theta,
    is_dwarf: is_dwarf,
    coeffs_used: is_dwarf ? "dwarf" : "giant",
    feh: feh,
    logg: logg,
    bc_g_regime: teff_K >= 4000 ? "hot (4000-8000K)" : "cool (3300-4000K)",
    dt: teff_K - T_EFF_SUN,
  };

  return result;
}

function computeBCG(teff_K) {
  let a;
  if (teff_K >= 4000 && teff_K <= 8000) {
    a = BC_COEFFS_HOT.a;
  } else if (teff_K >= 3300 && teff_K < 4000) {
    a = BC_COEFFS_COOL.a;
  } else {
    return null;
  }

  const dt = teff_K - 5772.0;
  return a[0] + a[1] * dt + a[2] * dt * dt + a[3] * dt * dt * dt + a[4] * dt * dt * dt * dt;
}
