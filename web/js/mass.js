// Module 3: Stellar mass estimation from the mass-luminosity relation.
// Ported from src/mass.py

function computeMass(star, teff_K, luminosity_Lsun) {
  const logg = star.logg != null ? star.logg : 4.0;
  const is_ms = (logg >= 3.5) && (teff_K >= 3300) && (teff_K <= 8000);

  if (!is_ms) {
    return {
      mass_Msun: null,
      mass_flag: "not_main_sequence",
      is_main_sequence: false,
      _intermediates: {
        logg: logg,
        teff_K: teff_K,
        reason: logg < 3.5 ? "logg < 3.5 (evolved star)" :
                teff_K < 3300 ? "Teff < 3300 K (too cool)" :
                "Teff > 8000 K (too hot for calibration)",
      },
    };
  }

  if (luminosity_Lsun == null) {
    return {
      mass_Msun: null,
      mass_flag: "no_luminosity",
      is_main_sequence: true,
      _intermediates: {reason: "No luminosity available"},
    };
  }

  const L = luminosity_Lsun;
  let alpha;
  let regime;
  if (L < 0.033) {
    alpha = 2.3;
    regime = "Very low-mass M dwarfs (L < 0.033 Lsun, M < 0.43 Msun)";
  } else if (L < 16) {
    alpha = 4.0;
    regime = "K/G/F dwarfs (0.033 < L < 16 Lsun, 0.43 < M < 2 Msun)";
  } else if (L < 54000) {
    alpha = 3.5;
    regime = "A-type and hotter (16 < L < 54000 Lsun, 2 < M < 55 Msun)";
  } else {
    alpha = 1.0;
    regime = "Massive stars (L > 54000 Lsun, M > 55 Msun)";
  }

  const mass_Msun = Math.pow(L, 1.0 / alpha);

  return {
    mass_Msun: mass_Msun,
    mass_flag: "ok",
    is_main_sequence: true,
    _intermediates: {
      luminosity_Lsun: L,
      alpha: alpha,
      regime: regime,
    },
  };
}
