"""Phase 7 target catalog: 10 nearest Sun-like stars for exo-Earth search.

Provides a curated target list with cross-matched identifiers (HD, HIP, TIC,
Gaia DR3) and functions to validate targets against Gaia DR3 and check data
availability in TESS/MAST archives.
"""

import logging

logger = logging.getLogger(__name__)

# Phase 7 target catalog -- 10 prioritized stars from phase7_exo_earth_poc.md
# Section 3.5. Ordered by priority (proximity + solar similarity + data richness).
#
# Identifiers will be populated/validated at runtime via SIMBAD resolution
# for any missing Gaia DR3 source_id or TIC values.
TARGET_CATALOG = [
    {
        "name": "Alpha Centauri A",
        "hd": "HD 128620",
        "hip": 71683,
        "tic": None,  # resolved at runtime; TESS saturated for this target
        "gaia_dr3_id": None,  # resolved at runtime via SIMBAD
        "spectral_type": "G2V",
        "teff_k": 5790,
        "distance_pc": 1.34,
        "mass_msun": 1.10,
        "feh": 0.20,
        "tier": "A/B",
        "known_planets": 1,
        "planet_status": "1 unconfirmed candidate",
        "notes": "Nearest G2V; binary with alpha Cen B; TESS saturated",
    },
    {
        "name": "Tau Ceti",
        "hd": "HD 10700",
        "hip": 8102,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G8.5V",
        "teff_k": 5375,
        "distance_pc": 3.65,
        "mass_msun": 0.78,
        "feh": -0.50,
        "tier": "C",
        "known_planets": 4,
        "planet_status": "4 unconfirmed candidates (e,f,g,h)",
        "notes": "9000+ HARPS measurements; low metallicity; debris disk",
    },
    {
        "name": "82 G. Eridani",
        "hd": "HD 20794",
        "hip": 15510,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G6V",
        "teff_k": 5401,
        "distance_pc": 6.04,
        "mass_msun": 0.70,
        "feh": -0.40,
        "tier": "C",
        "known_planets": 3,
        "planet_status": "3 confirmed (d in HZ, 5.8 Mearth)",
        "notes": "Primary validation target; HARPS + ESPRESSO data",
    },
    {
        "name": "Eta Cassiopeiae A",
        "hd": "HD 4614",
        "hip": 3821,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G0V",
        "teff_k": 6012,
        "distance_pc": 5.92,
        "mass_msun": 1.03,
        "feh": -0.30,
        "tier": "B",
        "known_planets": 1,
        "planet_status": "1 unconfirmed (~22 Mearth)",
        "notes": "Wide binary (70+ AU separation)",
    },
    {
        "name": "Delta Pavonis",
        "hd": "HD 190248",
        "hip": 99240,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G8IV-V",
        "teff_k": 5604,
        "distance_pc": 6.11,
        "mass_msun": 1.05,
        "feh": 0.33,
        "tier": "B/C",
        "known_planets": 0,
        "planet_status": "None (deep RV limits K < 1 m/s)",
        "notes": "CORALIE + HARPS since 1998; subgiant boundary",
    },
    {
        "name": "61 Virginis",
        "hd": "HD 115617",
        "hip": 64924,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G7V",
        "teff_k": 5531,
        "distance_pc": 8.53,
        "mass_msun": 0.94,
        "feh": -0.02,
        "tier": "C",
        "known_planets": 3,
        "planet_status": "3 confirmed (b, c, d)",
        "notes": "Near-solar metallicity; multi-planet system",
    },
    {
        "name": "Beta CVn",
        "hd": "HD 109358",
        "hip": 61317,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G0V",
        "teff_k": 5880,
        "distance_pc": 8.44,
        "mass_msun": 1.02,
        "feh": -0.21,
        "tier": "B",
        "known_planets": 0,
        "planet_status": "None",
        "notes": "Also known as Chara; quiet single star",
    },
    {
        "name": "Zeta Tucanae",
        "hd": "HD 1581",
        "hip": 1599,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "F9.5V",
        "teff_k": 5956,
        "distance_pc": 8.60,
        "mass_msun": 0.99,
        "feh": -0.15,
        "tier": "B",
        "known_planets": 0,
        "planet_status": "None",
        "notes": "At upper Teff boundary; slow rotator confirmed",
    },
    {
        "name": "18 Scorpii",
        "hd": "HD 146233",
        "hip": 79672,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G2Va",
        "teff_k": 5817,
        "distance_pc": 14.13,
        "mass_msun": 1.02,
        "feh": 0.05,
        "tier": "A",
        "known_planets": 1,
        "planet_status": "1 candidate (super-Earth, P~19.9d)",
        "notes": "Nearest solar twin; beyond 10 pc but high value",
    },
    {
        "name": "HD 134060",
        "hd": "HD 134060",
        "hip": 74273,
        "tic": None,
        "gaia_dr3_id": None,
        "spectral_type": "G0V",
        "teff_k": 5890,
        "distance_pc": 24.15,
        "mass_msun": 1.07,
        "feh": 0.10,
        "tier": "B",
        "known_planets": 2,
        "planet_status": "2 confirmed (b: 3.3d, c: 1168d)",
        "notes": "Beyond 10 pc; long-period planet c in wide orbit",
    },
]


def get_targets():
    """Return the full Phase 7 target catalog.

    Returns
    -------
    list[dict]
        List of target dictionaries.
    """
    return list(TARGET_CATALOG)


def get_target(name):
    """Look up a target by common name, HD number, or HIP number.

    Parameters
    ----------
    name : str
        Target identifier (case-insensitive). Matches against 'name', 'hd',
        or 'HIP <number>'.

    Returns
    -------
    dict or None
        Target dictionary, or None if not found.
    """
    name_lower = name.strip().lower()
    for target in TARGET_CATALOG:
        if target["name"].lower() == name_lower:
            return dict(target)
        if target["hd"].lower() == name_lower:
            return dict(target)
        if f"hip {target['hip']}" == name_lower:
            return dict(target)
    logger.warning("Target '%s' not found in catalog", name)
    return None


def resolve_target_ids(target):
    """Resolve missing identifiers (Gaia DR3, TIC) via SIMBAD.

    Queries SIMBAD for the target and extracts Gaia DR3 source_id and
    TIC ID from the identifier list.

    Parameters
    ----------
    target : dict
        Target dictionary from the catalog.

    Returns
    -------
    dict
        Updated target dictionary with resolved IDs. Unresolved IDs remain None.
    """
    from astroquery.simbad import Simbad
    import re

    result = dict(target)
    name = target["name"]
    logger.info("Resolving identifiers for %s", name)

    # Try common name first, then HD number as fallback
    search_names = [name]
    if target.get("hd") and target["hd"] != name:
        search_names.append(target["hd"])

    id_table = None
    for search_name in search_names:
        try:
            id_table = Simbad.query_objectids(search_name)
            if id_table is not None and len(id_table) > 0:
                logger.info("SIMBAD resolved '%s' via query name '%s'", name, search_name)
                break
        except Exception as e:
            logger.warning("SIMBAD query_objectids failed for '%s': %s", search_name, e)

    if id_table is None or len(id_table) == 0:
        logger.warning("No SIMBAD identifiers found for '%s'", name)
        return result

    for row in id_table:
        id_str = str(row["id"]).strip()

        # Gaia DR3
        if result["gaia_dr3_id"] is None:
            match = re.match(r"Gaia DR3\s+(\d+)", id_str)
            if match:
                result["gaia_dr3_id"] = match.group(1)
                logger.info("  Gaia DR3: %s", result["gaia_dr3_id"])

        # TIC
        if result["tic"] is None:
            match = re.match(r"TIC\s+(\d+)", id_str)
            if match:
                result["tic"] = int(match.group(1))
                logger.info("  TIC: %d", result["tic"])

    if result["gaia_dr3_id"] is None:
        logger.warning("Could not resolve Gaia DR3 ID for '%s'", name)
    if result["tic"] is None:
        logger.info("No TIC ID found for '%s' (may be too bright for TESS)", name)

    return result


def validate_against_gaia(target):
    """Validate catalog values against Gaia DR3.

    Queries Gaia DR3 for the target and compares Teff and distance
    against the catalog values.

    Parameters
    ----------
    target : dict
        Target dictionary (must have 'gaia_dr3_id' resolved, or 'name' for
        SIMBAD fallback).

    Returns
    -------
    dict
        Validation report with keys: gaia_data (raw Gaia result),
        teff_match (bool), teff_catalog, teff_gaia, teff_diff_pct,
        distance_match (bool), distance_catalog, distance_gaia, distance_diff_pct.
        Returns None if Gaia query fails.
    """
    from src.data_access import query_stars_by_id, resolve_simbad_name

    gaia_id = target.get("gaia_dr3_id")
    if gaia_id is None:
        # Try resolving via SIMBAD
        gaia_id = resolve_simbad_name(target["name"])
        if gaia_id is None:
            logger.warning("Cannot validate %s: no Gaia DR3 ID", target["name"])
            return None

    logger.info("Validating %s against Gaia DR3 (source_id=%s)", target["name"], gaia_id)
    stars = query_stars_by_id([gaia_id])

    if not stars:
        logger.warning("Gaia DR3 returned no data for source_id=%s", gaia_id)
        return None

    gaia = stars[0]
    report = {"gaia_data": gaia}

    # Teff comparison
    teff_gaia = gaia.get("teff_gspphot")
    teff_catalog = target.get("teff_k")
    if teff_gaia is not None and teff_catalog is not None:
        diff_pct = abs(teff_gaia - teff_catalog) / teff_catalog * 100
        report["teff_catalog"] = teff_catalog
        report["teff_gaia"] = teff_gaia
        report["teff_diff_pct"] = round(diff_pct, 1)
        report["teff_match"] = diff_pct < 5.0  # within 5%
    else:
        report["teff_match"] = None

    # Distance comparison
    dist_gaia = gaia.get("ref_distance_pc")
    dist_catalog = target.get("distance_pc")
    if dist_gaia is not None and dist_catalog is not None:
        diff_pct = abs(dist_gaia - dist_catalog) / dist_catalog * 100
        report["distance_catalog"] = dist_catalog
        report["distance_gaia"] = dist_gaia
        report["distance_diff_pct"] = round(diff_pct, 1)
        report["distance_match"] = diff_pct < 10.0  # within 10%
    else:
        report["distance_match"] = None

    logger.info("Validation for %s: Teff match=%s, distance match=%s",
                target["name"], report.get("teff_match"), report.get("distance_match"))
    return report


def check_tess_availability(target):
    """Check TESS light curve availability for a target.

    Parameters
    ----------
    target : dict
        Target dictionary.

    Returns
    -------
    dict
        Availability report with keys: available (bool), n_sectors (int),
        mission, author. Returns None if query fails.
    """
    from src.lightcurve import search_lightcurve

    # Try TIC ID first, then HD name, then common name
    search_names = []
    if target.get("tic"):
        search_names.append(f"TIC {target['tic']}")
    search_names.append(target["hd"])
    search_names.append(target["name"])

    for search_name in search_names:
        logger.info("Checking TESS availability for %s (searching as '%s')",
                     target["name"], search_name)
        try:
            info = search_lightcurve(search_name, mission="TESS")
            if info is not None:
                return {
                    "available": True,
                    "n_sectors": info["n_available"],
                    "mission": info["mission"],
                    "author": info["author"],
                    "search_name": search_name,
                }
        except Exception as e:
            logger.warning("TESS search failed for '%s': %s", search_name, e)

    logger.info("No TESS data found for %s", target["name"])
    return {"available": False, "n_sectors": 0, "mission": None, "author": None}


def resolve_all_targets():
    """Resolve identifiers for all targets in the catalog.

    Returns
    -------
    list[dict]
        List of target dictionaries with resolved Gaia DR3 and TIC IDs.
    """
    resolved = []
    for target in TARGET_CATALOG:
        resolved_target = resolve_target_ids(target)
        resolved.append(resolved_target)
    return resolved
