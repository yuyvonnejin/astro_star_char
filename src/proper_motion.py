"""Phase 7: Gaia-Hipparcos proper motion anomaly (PMa) analysis.

The proper motion anomaly is the difference between Gaia DR3 proper motion
(quasi-instantaneous, epoch ~2016) and Hipparcos proper motion (epoch ~1991.25).
A significant PMa indicates an unseen companion pulling the star over the
~25-year baseline (Kervella et al. 2019, 2022).

This module queries both catalogs and computes PMa significance and companion
mass sensitivity.
"""

import logging
import numpy as np

logger = logging.getLogger(__name__)

# Hipparcos-Gaia epoch difference (years)
EPOCH_DIFF_YR = 24.75  # Gaia DR3 epoch 2016.0 - Hipparcos epoch 1991.25

# Constants for mass sensitivity
AU_M = 1.496e11  # 1 AU in meters
YR_S = 3.156e7  # 1 year in seconds
MSUN_KG = 1.989e30
MEARTH_KG = 5.972e24
MAS_TO_RAD = np.pi / (180.0 * 3600.0 * 1000.0)


def query_gaia_proper_motion(source_id):
    """Query Gaia DR3 for proper motion and astrometric parameters.

    Parameters
    ----------
    source_id : str or int
        Gaia DR3 source_id.

    Returns
    -------
    dict or None
        Astrometric parameters: pmra, pmdec, pmra_error, pmdec_error,
        parallax, ruwe, astrometric_excess_noise.
    """
    from astroquery.gaia import Gaia

    query = f"""
    SELECT source_id, ra, dec, parallax, parallax_error,
           pmra, pmdec, pmra_error, pmdec_error,
           ruwe, astrometric_excess_noise, astrometric_excess_noise_sig,
           ipd_gof_harmonic_amplitude, visibility_periods_used
    FROM gaiadr3.gaia_source
    WHERE source_id = {source_id}
    """

    logger.info("Querying Gaia DR3 astrometry for source_id=%s", source_id)
    try:
        job = Gaia.launch_job(query)
        result = job.get_results()
        if len(result) == 0:
            logger.warning("No Gaia DR3 data for source_id=%s", source_id)
            return None

        row = result[0]
        data = {}
        for col in result.colnames:
            val = row[col]
            try:
                if hasattr(val, "mask") and val.mask:
                    data[col] = None
                elif np.isnan(float(val)):
                    data[col] = None
                else:
                    data[col] = float(val)
            except (TypeError, ValueError):
                data[col] = val

        logger.info("Gaia DR3 PM for source_id=%s: pmra=%.3f, pmdec=%.3f mas/yr",
                     source_id, data.get("pmra", 0), data.get("pmdec", 0))
        return data

    except Exception as e:
        logger.warning("Gaia DR3 query failed for source_id=%s: %s", source_id, e)
        return None


def query_hipparcos_proper_motion(hip_id):
    """Query Hipparcos-2 catalog (van Leeuwen 2007) for proper motion.

    Parameters
    ----------
    hip_id : int
        Hipparcos catalog number.

    Returns
    -------
    dict or None
        Hipparcos PM: pmra, pmdec, pmra_error, pmdec_error, parallax.
    """
    from astroquery.vizier import Vizier

    logger.info("Querying Hipparcos-2 for HIP %d", hip_id)
    try:
        # Hipparcos-2 catalog: I/311/hip2
        vizier = Vizier(columns=["HIP", "Plx", "pmRA", "pmDE", "e_pmRA", "e_pmDE"])
        result = vizier.query_constraints(catalog="I/311/hip2", HIP=hip_id)

        if not result or len(result) == 0 or len(result[0]) == 0:
            logger.warning("No Hipparcos-2 data for HIP %d", hip_id)
            return None

        row = result[0][0]
        data = {
            "hip_id": hip_id,
            "pmra": float(row["pmRA"]),
            "pmdec": float(row["pmDE"]),
            "pmra_error": float(row["e_pmRA"]),
            "pmdec_error": float(row["e_pmDE"]),
            "parallax": float(row["Plx"]),
        }

        logger.info("Hipparcos-2 PM for HIP %d: pmra=%.3f, pmdec=%.3f mas/yr",
                     hip_id, data["pmra"], data["pmdec"])
        return data

    except Exception as e:
        logger.warning("Hipparcos-2 query failed for HIP %d: %s", hip_id, e)
        return None


def compute_pma(gaia_pm, hip_pm):
    """Compute the proper motion anomaly between Gaia and Hipparcos.

    The PMa is the vector difference between the two proper motions. A
    significant PMa indicates the star's velocity changed between epochs,
    likely due to an unseen companion.

    Parameters
    ----------
    gaia_pm : dict
        Gaia proper motion with keys: pmra, pmdec, pmra_error, pmdec_error.
    hip_pm : dict
        Hipparcos proper motion with same keys.

    Returns
    -------
    dict
        PMa results: delta_pmra (mas/yr), delta_pmdec (mas/yr),
        pma_total (mas/yr), pma_error (mas/yr), pma_snr,
        significant (bool, SNR > 3).
    """
    dpmra = gaia_pm["pmra"] - hip_pm["pmra"]
    dpmdec = gaia_pm["pmdec"] - hip_pm["pmdec"]

    # Error propagation (assuming independent)
    err_pmra = np.sqrt(gaia_pm["pmra_error"]**2 + hip_pm["pmra_error"]**2)
    err_pmdec = np.sqrt(gaia_pm["pmdec_error"]**2 + hip_pm["pmdec_error"]**2)

    # Total PMa magnitude
    pma_total = np.sqrt(dpmra**2 + dpmdec**2)

    # Error on total PMa (first-order propagation)
    if pma_total > 0:
        pma_error = np.sqrt(
            (dpmra * err_pmra)**2 + (dpmdec * err_pmdec)**2
        ) / pma_total
    else:
        pma_error = np.sqrt(err_pmra**2 + err_pmdec**2)

    snr = pma_total / pma_error if pma_error > 0 else 0.0

    result = {
        "delta_pmra_mas_yr": round(float(dpmra), 4),
        "delta_pmdec_mas_yr": round(float(dpmdec), 4),
        "pma_total_mas_yr": round(float(pma_total), 4),
        "pma_error_mas_yr": round(float(pma_error), 4),
        "pma_snr": round(float(snr), 2),
        "significant": snr > 3.0,
    }

    logger.info("PMa: dpmra=%.4f, dpmdec=%.4f, total=%.4f mas/yr, SNR=%.1f (%s)",
                dpmra, dpmdec, pma_total, snr,
                "SIGNIFICANT" if result["significant"] else "not significant")

    return result


def pma_companion_mass(pma_mas_yr, distance_pc, stellar_mass_msun, separation_au=None):
    """Estimate companion mass from PMa magnitude.

    For long-period companions (P >> Delta_t), the PMa reflects the
    gravitational acceleration of the star by the companion. The change
    in tangential velocity over the Gaia-Hipparcos baseline is:

        Delta_v = PMa * d

    And the gravitational acceleration at separation a is:

        a_star = G * M_c / a^2

    Setting Delta_v = a_star * Delta_t gives:

        M_c = PMa * d * a^2 / (G * Delta_t)

    In natural units (AU, yr, Msun): G = 4*pi^2.

    Parameters
    ----------
    pma_mas_yr : float
        PMa magnitude in mas/yr.
    distance_pc : float
        Distance to the star in parsecs.
    stellar_mass_msun : float
        Stellar mass in solar masses (used for context, not in this formula).
    separation_au : float or array, optional
        Assumed orbital separation(s) in AU. Default: grid from 1-50 AU.

    Returns
    -------
    dict
        Companion mass estimates: separation_au (array), mass_mjup (array),
        mass_mearth (array). For each assumed separation, gives the companion
        mass that would produce the observed PMa.
    """
    if separation_au is None:
        separation_au = np.array([1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0])

    separation_au = np.asarray(separation_au, dtype=float)
    pma_rad_yr = pma_mas_yr * MAS_TO_RAD

    # Distance in AU (1 pc = 206265 AU)
    d_au = distance_pc * 206265.0

    # G in AU^3 / (Msun * yr^2) = 4 * pi^2
    G_au = 4.0 * np.pi**2

    # M_c [Msun] = PMa [rad/yr] * d [AU] * a^2 [AU^2] / (G [AU^3/(Msun*yr^2)] * Delta_t [yr])
    m_companion_msun = (pma_rad_yr * d_au * separation_au**2) / (G_au * EPOCH_DIFF_YR)

    # Convert to Jupiter and Earth masses
    m_jupiter = m_companion_msun * (MSUN_KG / 1.898e27)  # Msun to Mjup
    m_earth = m_companion_msun * (MSUN_KG / MEARTH_KG)

    logger.info("PMa mass estimate at %.1f AU: %.1f Mjup for PMa=%.4f mas/yr at %.1f pc",
                float(separation_au[0]) if len(separation_au) > 0 else 0,
                float(m_jupiter[0]) if len(m_jupiter) > 0 else 0,
                pma_mas_yr, distance_pc)

    return {
        "separation_au": separation_au,
        "mass_msun": m_companion_msun,
        "mass_mjup": m_jupiter,
        "mass_mearth": m_earth,
    }


def analyze_pma(gaia_source_id, hip_id, distance_pc=None, stellar_mass_msun=None):
    """Full PMa analysis pipeline for a target.

    Queries Gaia DR3 and Hipparcos, computes PMa, and estimates companion
    mass sensitivity.

    Parameters
    ----------
    gaia_source_id : str or int
        Gaia DR3 source_id.
    hip_id : int
        Hipparcos catalog number.
    distance_pc : float, optional
        Distance for mass estimation. If None, uses Gaia parallax.
    stellar_mass_msun : float, optional
        Stellar mass for mass estimation. Default: 1.0.

    Returns
    -------
    dict
        Full PMa analysis: gaia_pm, hipparcos_pm, pma (from compute_pma),
        companion_mass (if PMa significant), ruwe, astrometric_excess_noise.
    """
    if stellar_mass_msun is None:
        stellar_mass_msun = 1.0

    result = {
        "gaia_source_id": str(gaia_source_id),
        "hip_id": hip_id,
    }

    # Query both catalogs
    gaia_pm = query_gaia_proper_motion(gaia_source_id)
    hip_pm = query_hipparcos_proper_motion(hip_id)

    if gaia_pm is None:
        result["status"] = "gaia_query_failed"
        return result
    if hip_pm is None:
        result["status"] = "hipparcos_query_failed"
        return result

    result["gaia_pm"] = gaia_pm
    result["hipparcos_pm"] = hip_pm

    # Gaia astrometric quality indicators
    result["ruwe"] = gaia_pm.get("ruwe")
    result["astrometric_excess_noise"] = gaia_pm.get("astrometric_excess_noise")

    # Compute PMa
    pma = compute_pma(gaia_pm, hip_pm)
    result["pma"] = pma

    # Distance from parallax if not provided
    if distance_pc is None:
        plx = gaia_pm.get("parallax")
        if plx is not None and plx > 0:
            distance_pc = 1000.0 / plx

    # Companion mass estimate if PMa is significant
    if pma["significant"] and distance_pc is not None:
        mass_est = pma_companion_mass(
            pma["pma_total_mas_yr"], distance_pc, stellar_mass_msun
        )
        result["companion_mass_estimate"] = {
            "separation_au": mass_est["separation_au"].tolist(),
            "mass_mjup": mass_est["mass_mjup"].tolist(),
            "mass_mearth": mass_est["mass_mearth"].tolist(),
        }

    result["status"] = "ok"
    return result


def assess_ruwe(ruwe_value):
    """Interpret the Gaia RUWE (Renormalized Unit Weight Error).

    RUWE > 1.4 suggests the single-star astrometric model is a poor fit,
    which can indicate an unresolved companion.

    Parameters
    ----------
    ruwe_value : float
        RUWE from Gaia DR3.

    Returns
    -------
    dict
        Assessment: value, interpretation, companion_hint (bool).
    """
    if ruwe_value is None:
        return {"value": None, "interpretation": "not available", "companion_hint": False}

    if ruwe_value > 1.4:
        return {
            "value": round(ruwe_value, 3),
            "interpretation": "elevated - possible unresolved companion or extended source",
            "companion_hint": True,
        }
    else:
        return {
            "value": round(ruwe_value, 3),
            "interpretation": "normal - consistent with single star",
            "companion_hint": False,
        }
