"""Module 2: Effective temperature, bolometric correction, and luminosity."""

import logging
import numpy as np

logger = logging.getLogger(__name__)

# Mucciarelli, Bellazzini & Massari (2021) -- Table 1
DWARF_COEFFS = {
    "b": [0.4929, 0.5092, -0.0353, 0.0192, -0.0020, -0.0395],
    "color_range": (0.39, 1.50),
    "dispersion_K": 61,
}
GIANT_COEFFS = {
    "b": [0.5323, 0.4775, -0.0344, -0.0110, -0.0020, -0.0009],
    "color_range": (0.33, 1.81),
    "dispersion_K": 83,
}

# Andrae et al. (2018) -- BC_G polynomial coefficients
BC_COEFFS_HOT = {
    "a": [6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15],
    "teff_range": (4000, 8000),
}
BC_COEFFS_COOL = {
    "a": [1.749e+00, 1.977e-03, 3.737e-07, -8.966e-11, -4.183e-14],
    "teff_range": (3300, 4000),
}

M_BOL_SUN = 4.74  # IAU 2015 B2
T_EFF_SUN = 5772.0  # K


def compute_temperature_luminosity(star_dict, distance_pc):
    """Compute Teff, BC_G, and luminosity for a star.

    Parameters
    ----------
    star_dict : dict
        Star record with photometric and spectroscopic fields.
    distance_pc : float
        Distance in parsecs from Module 1.

    Returns
    -------
    dict
        Temperature, luminosity, and diagnostic fields.
    """
    result = {}

    # Step 1: Deredden
    bp_rp_0 = star_dict["bp_rp"] - star_dict.get("ebpminrp_gspphot", 0.0)
    g_0 = star_dict["phot_g_mean_mag"] - star_dict.get("ag_gspphot", 0.0)
    result["bp_rp_0"] = bp_rp_0

    # Step 2: Effective temperature
    logg = star_dict.get("logg", 4.0)
    feh = star_dict.get("feh", 0.0)
    is_dwarf = logg >= 3.0

    coeffs = DWARF_COEFFS if is_dwarf else GIANT_COEFFS
    b = coeffs["b"]
    color_lo, color_hi = coeffs["color_range"]

    teff_flag = "ok"
    if bp_rp_0 < color_lo or bp_rp_0 > color_hi:
        teff_flag = "outside_valid_range"
        logger.warning(
            "Color (BP-RP)_0=%.3f outside valid range [%.2f, %.2f]",
            bp_rp_0, color_lo, color_hi,
        )

    c = bp_rp_0
    theta = b[0] + b[1] * c + b[2] * c**2 + b[3] * feh + b[4] * feh**2 + b[5] * feh * c
    teff_K = 5040.0 / theta

    result["teff_K"] = teff_K
    result["teff_uncertainty_K"] = coeffs["dispersion_K"]
    result["teff_flag"] = teff_flag

    # Step 3: Absolute G magnitude
    m_g = g_0 - 5.0 * np.log10(distance_pc) + 5.0
    result["M_G"] = m_g

    # Step 4: Bolometric correction BC_G
    bc_g = _compute_bc_g(teff_K)
    if bc_g is None:
        logger.warning("Teff=%.0f K outside BC_G range [3300, 8000]; skipping luminosity", teff_K)
        result["BC_G"] = None
        result["M_bol"] = None
        result["luminosity_Lsun"] = None
        result["radius_Rsun"] = None
        result["luminosity_validation_ratio"] = None
        return result

    result["BC_G"] = bc_g

    # Step 5: Luminosity
    m_bol = m_g + bc_g
    luminosity = 10.0 ** ((M_BOL_SUN - m_bol) / 2.5)
    result["M_bol"] = m_bol
    result["luminosity_Lsun"] = luminosity

    # Step 6: Radius via Stefan-Boltzmann law
    # R/R_sun = sqrt(L/L_sun) * (T_sun / T_eff)^2
    radius_Rsun = luminosity**0.5 * (T_EFF_SUN / teff_K)**2
    result["radius_Rsun"] = radius_Rsun

    # Step 7: Validation cross-check
    lum_gspphot = star_dict.get("lum_gspphot")
    if lum_gspphot is not None and lum_gspphot > 0:
        ratio = luminosity / lum_gspphot
        result["luminosity_validation_ratio"] = ratio
        logger.info("Luminosity validation ratio: %.3f", ratio)
    else:
        result["luminosity_validation_ratio"] = None

    logger.info(
        "Teff=%.0f K, M_G=%.3f, BC_G=%.4f, L=%.5f Lsun, R=%.4f Rsun",
        teff_K, m_g, bc_g, luminosity, radius_Rsun,
    )
    return result


def _compute_bc_g(teff_K):
    """Compute bolometric correction BC_G from Teff.

    Returns None if Teff is outside the valid range [3300, 8000] K.
    """
    if 4000 <= teff_K <= 8000:
        a = BC_COEFFS_HOT["a"]
    elif 3300 <= teff_K < 4000:
        a = BC_COEFFS_COOL["a"]
    else:
        return None

    dt = teff_K - 5772.0
    bc_g = a[0] + a[1] * dt + a[2] * dt**2 + a[3] * dt**3 + a[4] * dt**4
    return bc_g
