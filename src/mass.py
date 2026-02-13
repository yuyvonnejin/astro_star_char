"""Module 3: Stellar mass estimation from the mass-luminosity relation."""

import logging

logger = logging.getLogger(__name__)


def compute_mass(star_dict, teff_K, luminosity_Lsun):
    """Compute stellar mass using the piecewise mass-luminosity relation.

    Only valid for main-sequence stars.

    Parameters
    ----------
    star_dict : dict
        Star record (needs logg).
    teff_K : float
        Effective temperature in Kelvin from Module 2.
    luminosity_Lsun : float or None
        Luminosity in solar units from Module 2.

    Returns
    -------
    dict
        Mass estimate and flags.
    """
    logg = star_dict.get("logg", 4.0)

    is_ms = (logg >= 3.5) and (teff_K >= 3300) and (teff_K <= 8000)

    if not is_ms:
        logger.info(
            "Not main-sequence (logg=%.2f, Teff=%.0f K); skipping mass",
            logg, teff_K,
        )
        return {
            "mass_Msun": None,
            "mass_flag": "not_main_sequence",
            "is_main_sequence": False,
        }

    if luminosity_Lsun is None:
        logger.warning("Luminosity unavailable; cannot compute mass")
        return {
            "mass_Msun": None,
            "mass_flag": "no_luminosity",
            "is_main_sequence": True,
        }

    L = luminosity_Lsun
    if L < 0.033:
        alpha = 2.3
    elif L < 16:
        alpha = 4.0
    elif L < 54000:
        alpha = 3.5
    else:
        alpha = 1.0

    mass_Msun = L ** (1.0 / alpha)

    logger.info("Mass estimate: L=%.4f Lsun, alpha=%.1f, M=%.3f Msun", L, alpha, mass_Msun)
    return {
        "mass_Msun": mass_Msun,
        "mass_flag": "ok",
        "is_main_sequence": True,
    }
