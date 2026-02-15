"""Data access layer: query Gaia DR3 via astroquery TAP+ service."""

import logging

logger = logging.getLogger(__name__)


def resolve_simbad_name(name):
    """Resolve a SIMBAD name (e.g. 'Vega', 'Alp Lyr', 'HD 172167') to a Gaia DR3 source_id.

    Queries SIMBAD for all known identifiers and looks for a Gaia DR3 match.
    Falls back to coordinate-based Gaia cone search if no direct ID match.

    Returns
    -------
    str or None
        Gaia DR3 source_id, or None if resolution failed.
    """
    from astroquery.simbad import Simbad
    import re

    logger.info("Resolving SIMBAD name: %s", name)

    # Try to find Gaia DR3 ID directly from SIMBAD identifiers
    try:
        id_table = Simbad.query_objectids(name)
    except Exception as e:
        logger.warning("SIMBAD query_objectids failed for '%s': %s", name, e)
        id_table = None

    if id_table is not None:
        for row in id_table:
            id_str = str(row["id"]).strip()
            match = re.match(r"Gaia DR3\s+(\d+)", id_str)
            if match:
                source_id = match.group(1)
                logger.info("Resolved '%s' -> Gaia DR3 %s", name, source_id)
                return source_id

    # Fallback: get coordinates from SIMBAD, cone search Gaia
    logger.info("No Gaia DR3 ID in SIMBAD identifiers for '%s'; trying coordinate cone search", name)
    try:
        from astroquery.gaia import Gaia

        obj_table = Simbad.query_object(name)
        if obj_table is None or len(obj_table) == 0:
            logger.warning("SIMBAD cannot resolve '%s'", name)
            return None

        ra = float(obj_table["ra"][0])
        dec = float(obj_table["dec"][0])

        query = f"""
        SELECT TOP 1 source_id, parallax
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT(ra, dec),
            CIRCLE({ra}, {dec}, 0.005)
        ) = 1
        ORDER BY phot_g_mean_mag ASC
        """
        job = Gaia.launch_job(query)
        result = job.get_results()
        if len(result) > 0:
            source_id = str(result["source_id"][0])
            logger.info("Resolved '%s' via cone search -> Gaia DR3 %s", name, source_id)
            return source_id
    except Exception as e:
        logger.warning("Coordinate fallback failed for '%s': %s", name, e)

    logger.warning("Could not resolve '%s' to a Gaia DR3 source_id", name)
    return None


def query_stars_by_id(source_ids):
    """Query Gaia DR3 for specific stars by source_id.

    Parameters
    ----------
    source_ids : list[int or str]
        Gaia source IDs to retrieve.

    Returns
    -------
    list[dict]
        Star records in pipeline input schema format.
    """
    from astroquery.gaia import Gaia

    id_list = ", ".join(str(sid) for sid in source_ids)
    query = f"""
    SELECT g.source_id, g.parallax, g.parallax_error,
           g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
           g.bp_rp, g.teff_gspphot,
           g.ag_gspphot, g.ebpminrp_gspphot,
           ap.mh_gspphot, ap.logg_gspphot, ap.lum_flame
    FROM gaiadr3.gaia_source AS g
    LEFT JOIN gaiadr3.astrophysical_parameters AS ap
      ON g.source_id = ap.source_id
    WHERE g.source_id IN ({id_list})
    """

    logger.info("Querying Gaia DR3 for %d source(s)", len(source_ids))
    job = Gaia.launch_job(query)
    table = job.get_results()

    stars = []
    for row in table:
        star = {
            "source_id": str(row["source_id"]),
            "parallax_mas": _val(row, "parallax"),
            "parallax_error_mas": _val(row, "parallax_error"),
            "phot_g_mean_mag": _val(row, "phot_g_mean_mag"),
            "phot_bp_mean_mag": _val(row, "phot_bp_mean_mag"),
            "phot_rp_mean_mag": _val(row, "phot_rp_mean_mag"),
            "bp_rp": _val(row, "bp_rp"),
            "ag_gspphot": _val(row, "ag_gspphot"),
            "ebpminrp_gspphot": _val(row, "ebpminrp_gspphot"),
            "feh": _val(row, "mh_gspphot"),
            "logg": _val(row, "logg_gspphot"),
            "teff_gspphot": _val(row, "teff_gspphot"),
            "lum_gspphot": _val(row, "lum_flame"),
            "is_cepheid": False,
            "cepheid_period_days": None,
        }
        stars.append(star)

    logger.info("Retrieved %d star(s)", len(stars))
    return stars


def query_cepheids(limit=100):
    """Query Gaia DR3 for classical Cepheids (DCEP type).

    Parameters
    ----------
    limit : int
        Maximum number of Cepheids to return.

    Returns
    -------
    list[dict]
        Cepheid star records in pipeline input schema format.
    """
    from astroquery.gaia import Gaia

    query = f"""
    SELECT TOP {limit} g.source_id, g.parallax, g.parallax_error,
           g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
           g.bp_rp, g.ag_gspphot, g.ebpminrp_gspphot,
           g.teff_gspphot,
           c.pf, c.type_best_classification
    FROM gaiadr3.gaia_source AS g
    JOIN gaiadr3.vari_cepheid AS c
      ON g.source_id = c.source_id
    WHERE c.type_best_classification = 'DCEP'
      AND g.parallax_over_error > 10
    """

    logger.info("Querying Gaia DR3 for up to %d Cepheids", limit)
    job = Gaia.launch_job(query)
    table = job.get_results()

    stars = []
    for row in table:
        pf = _val(row, "pf")
        period = (1.0 / pf) if pf is not None and pf > 0 else None
        star = {
            "source_id": str(row["source_id"]),
            "parallax_mas": _val(row, "parallax"),
            "parallax_error_mas": _val(row, "parallax_error"),
            "phot_g_mean_mag": _val(row, "phot_g_mean_mag"),
            "phot_bp_mean_mag": _val(row, "phot_bp_mean_mag"),
            "phot_rp_mean_mag": _val(row, "phot_rp_mean_mag"),
            "bp_rp": _val(row, "bp_rp"),
            "ag_gspphot": _val(row, "ag_gspphot"),
            "ebpminrp_gspphot": _val(row, "ebpminrp_gspphot"),
            "feh": None,
            "logg": None,
            "teff_gspphot": _val(row, "teff_gspphot"),
            "lum_gspphot": None,
            "is_cepheid": True,
            "cepheid_period_days": period,
        }
        stars.append(star)

    logger.info("Retrieved %d Cepheid(s)", len(stars))
    return stars


def _val(row, col):
    """Extract a column value, converting masked/NaN to None."""
    import numpy as np

    try:
        v = row[col]
    except KeyError:
        return None
    if hasattr(v, "mask") and v.mask:
        return None
    try:
        if np.isnan(v):
            return None
    except (TypeError, ValueError):
        pass
    return float(v)
