"""Module 4a: Light curve retrieval from MAST via lightkurve.

Searches for TESS, Kepler, or K2 light curves, downloads and stitches
multi-sector/quarter data, cleans, and flattens for period analysis.
"""

import logging
import numpy as np

logger = logging.getLogger(__name__)

# Mission search priority: TESS has broadest sky coverage, then Kepler, then K2.
DEFAULT_AUTHOR_PRIORITY = ["SPOC", "Kepler", "K2"]


def search_lightcurve(target, author_priority=None, mission=None):
    """Search MAST for available light curves for a target.

    Parameters
    ----------
    target : str
        Star identifier -- name ("Proxima Cen"), KIC/TIC ID ("KIC 6922244"),
        or coordinates ("19:12:34 +42:15:00").
    author_priority : list[str], optional
        Ordered list of authors/pipelines to try. Default: SPOC, Kepler, K2.
    mission : str, optional
        Force a specific mission ("TESS", "Kepler", "K2"). If set, overrides
        author_priority.

    Returns
    -------
    dict
        Search result info with keys: search_result (lightkurve object),
        mission, author, n_available. Returns None if nothing found.
    """
    import lightkurve as lk

    if author_priority is None:
        author_priority = DEFAULT_AUTHOR_PRIORITY

    if mission is not None:
        logger.info("Searching MAST for '%s' (mission=%s)", target, mission)
        search_res = lk.search_lightcurve(target, mission=mission)
        if len(search_res) > 0:
            author = search_res.author[0] if hasattr(search_res, "author") else mission
            logger.info("Found %d light curve(s) for '%s' (%s)", len(search_res), target, mission)
            return {
                "search_result": search_res,
                "mission": mission,
                "author": str(author),
                "n_available": len(search_res),
            }
        logger.info("No light curves found for '%s' (mission=%s)", target, mission)
        return None

    # Try each author in priority order
    for author in author_priority:
        logger.info("Searching MAST for '%s' (author=%s)", target, author)
        search_res = lk.search_lightcurve(target, author=author)
        if len(search_res) > 0:
            mission_name = search_res.mission[0] if hasattr(search_res, "mission") else author
            logger.info("Found %d light curve(s) for '%s' (%s)", len(search_res), target, author)
            return {
                "search_result": search_res,
                "mission": str(mission_name),
                "author": author,
                "n_available": len(search_res),
            }
        logger.info("No results for '%s' with author=%s", target, author)

    logger.warning("No light curves found for '%s' in any mission", target)
    return None


def download_and_stitch(search_info, max_sectors=20):
    """Download and stitch light curves into a single time series.

    Parameters
    ----------
    search_info : dict
        Output from search_lightcurve().
    max_sectors : int
        Maximum number of sectors/quarters to download. Set to 0 for all.

    Returns
    -------
    dict
        Light curve data with keys: time (array), flux (array), flux_err (array),
        mission, author, n_sectors, n_points_raw, cadence_s, time_baseline_days.
        Returns None if download fails.
    """
    search_result = search_info["search_result"]

    # Limit number of sectors if requested
    if max_sectors > 0 and len(search_result) > max_sectors:
        logger.info("Limiting download to %d of %d available sectors", max_sectors, len(search_result))
        search_result = search_result[:max_sectors]

    n_sectors = len(search_result)
    logger.info("Downloading %d light curve(s) for stitching", n_sectors)

    try:
        lc_collection = search_result.download_all()
    except Exception as e:
        logger.error("Download failed: %s", e)
        return None

    if lc_collection is None or len(lc_collection) == 0:
        logger.error("Download returned empty collection")
        return None

    logger.info("Downloaded %d light curve(s), stitching", len(lc_collection))
    lc = lc_collection.stitch()

    # Extract arrays
    time = lc.time.value
    flux = lc.flux.value
    flux_err = lc.flux_err.value if lc.flux_err is not None else np.zeros_like(flux)

    # Estimate cadence from median time differences
    dt = np.diff(time)
    cadence_days = float(np.median(dt[dt > 0]))
    cadence_s = cadence_days * 86400.0

    time_baseline = float(time[-1] - time[0]) if len(time) > 1 else 0.0

    logger.info("Stitched light curve: %d points, %.1f day baseline, %.0f s cadence",
                len(time), time_baseline, cadence_s)

    return {
        "time": time,
        "flux": flux,
        "flux_err": flux_err,
        "mission": search_info["mission"],
        "author": search_info["author"],
        "n_sectors": n_sectors,
        "n_points_raw": len(time),
        "cadence_s": round(cadence_s, 1),
        "time_baseline_days": round(time_baseline, 2),
    }


def clean_lightcurve(lc_data, sigma_clip=5.0):
    """Remove NaN values and outliers from light curve data.

    Parameters
    ----------
    lc_data : dict
        Output from download_and_stitch().
    sigma_clip : float
        Remove flux points deviating more than this many sigma from median.

    Returns
    -------
    dict
        Cleaned light curve data (same keys as input, updated arrays and n_points_clean).
    """
    time = lc_data["time"]
    flux = lc_data["flux"]
    flux_err = lc_data["flux_err"]

    # Remove NaNs
    mask = np.isfinite(flux) & np.isfinite(time)
    n_nan = int(np.sum(~mask))

    # Sigma clip outliers
    if np.any(mask):
        median_flux = np.median(flux[mask])
        std_flux = np.std(flux[mask])
        if std_flux > 0:
            outlier_mask = np.abs(flux - median_flux) > sigma_clip * std_flux
            n_outlier = int(np.sum(mask & outlier_mask))
            mask = mask & ~outlier_mask
        else:
            n_outlier = 0
    else:
        n_outlier = 0

    logger.info("Cleaning: removed %d NaN, %d outlier (%.1f-sigma) of %d points",
                n_nan, n_outlier, sigma_clip, len(flux))

    out = dict(lc_data)
    out["time"] = time[mask]
    out["flux"] = flux[mask]
    out["flux_err"] = flux_err[mask]
    out["n_points_clean"] = int(np.sum(mask))
    return out


def flatten_lightcurve(lc_data, window_length=1001):
    """Flatten light curve by removing long-term trends.

    Uses lightkurve's Savitzky-Golay based flatten method.

    Parameters
    ----------
    lc_data : dict
        Output from clean_lightcurve().
    window_length : int
        Window length for the SG filter in cadences. Must be odd.
        Larger = removes slower trends, preserves longer-period signals.
        Default 1001 is good for periods up to ~10 days at 2-min cadence.

    Returns
    -------
    dict
        Flattened light curve data (same keys, updated flux array, added flux_flat).
    """
    import lightkurve as lk

    # Reconstruct a lightkurve object for flattening
    lc_obj = lk.LightCurve(
        time=lc_data["time"],
        flux=lc_data["flux"],
        flux_err=lc_data["flux_err"],
    )

    # Ensure window_length is odd
    if window_length % 2 == 0:
        window_length += 1

    # Flatten -- window_length must not exceed data length
    n = len(lc_data["flux"])
    if window_length >= n:
        window_length = n - 1 if n % 2 == 0 else n - 2
        if window_length < 3:
            logger.warning("Too few points (%d) to flatten", n)
            out = dict(lc_data)
            out["flux_flat"] = lc_data["flux"].copy()
            out["window_length"] = 0
            return out

    lc_flat = lc_obj.flatten(window_length=window_length)

    out = dict(lc_data)
    out["flux_flat"] = lc_flat.flux.value
    out["window_length"] = window_length

    logger.info("Flattened light curve with window_length=%d", window_length)
    return out


def bin_lightcurve(lc_data, bin_cadence_s=120.0):
    """Bin a light curve to a coarser cadence by median-averaging.

    Reduces data volume for high-cadence data (e.g., 20s TESS) before
    expensive operations like Lomb-Scargle and BLS. Skipped if the
    native cadence is already >= bin_cadence_s.

    Parameters
    ----------
    lc_data : dict
        Light curve data with time, flux, flux_err arrays.
    bin_cadence_s : float
        Target cadence in seconds for binning.

    Returns
    -------
    dict
        Binned light curve data.
    """
    native_cadence = lc_data.get("cadence_s", 0)
    if native_cadence >= bin_cadence_s:
        logger.info("Cadence %.0fs already >= %.0fs, skipping binning",
                     native_cadence, bin_cadence_s)
        return lc_data

    time = lc_data["time"]
    flux = lc_data["flux"]
    flux_err = lc_data["flux_err"]

    # Ensure time-sorted order (lightkurve stitch can produce slight overlaps
    # at sector boundaries, breaking searchsorted)
    n_unsorted = int(np.sum(np.diff(time) < 0))
    if n_unsorted > 0:
        sort_idx = np.argsort(time)
        time = time[sort_idx]
        flux = flux[sort_idx]
        flux_err = flux_err[sort_idx]
        logger.info("Sorted %d out-of-order point(s) before binning", n_unsorted)

    bin_width_days = bin_cadence_s / 86400.0
    t_start = time[0]
    t_end = time[-1]
    n_bins = int((t_end - t_start) / bin_width_days) + 1

    bin_time = []
    bin_flux = []
    bin_err = []

    bin_edges = t_start + np.arange(n_bins + 1) * bin_width_days
    # Use searchsorted for fast binning (requires sorted time)
    indices = np.searchsorted(time, bin_edges)
    for i in range(n_bins):
        lo, hi = indices[i], indices[i + 1]
        if hi > lo:
            bin_time.append(np.median(time[lo:hi]))
            bin_flux.append(np.median(flux[lo:hi]))
            # Error of median ~ 1.253 * mean(err) / sqrt(n)
            n = hi - lo
            bin_err.append(1.253 * np.mean(flux_err[lo:hi]) / np.sqrt(n))

    n_before = len(time)
    n_after = len(bin_time)
    logger.info("Binned light curve: %d -> %d points (%.0fs -> %.0fs cadence)",
                n_before, n_after, native_cadence, bin_cadence_s)

    out = dict(lc_data)
    out["time"] = np.array(bin_time)
    out["flux"] = np.array(bin_flux)
    out["flux_err"] = np.array(bin_err)
    out["n_points_clean"] = n_after
    out["cadence_s"] = bin_cadence_s
    return out


def retrieve_lightcurve(target, mission=None, max_sectors=20, sigma_clip=5.0,
                        window_length=1001, bin_cadence_s=120.0):
    """Full pipeline: search, download, clean, bin, flatten.

    Convenience function that chains all steps.

    Parameters
    ----------
    target : str
        Star identifier for MAST search.
    mission : str, optional
        Force specific mission.
    max_sectors : int
        Max sectors to download.
    sigma_clip : float
        Outlier rejection threshold.
    window_length : int
        Flatten window length.
    bin_cadence_s : float
        Target cadence for binning (seconds). High-cadence data (e.g., 20s TESS)
        is binned to this cadence before flattening and analysis. Set to 0 to
        disable binning.

    Returns
    -------
    dict or None
        Complete light curve data ready for period analysis, or None if
        no data available.
    """
    search_info = search_lightcurve(target, mission=mission)
    if search_info is None:
        return None

    lc_data = download_and_stitch(search_info, max_sectors=max_sectors)
    if lc_data is None:
        return None

    lc_data = clean_lightcurve(lc_data, sigma_clip=sigma_clip)

    # Bin high-cadence data before expensive operations (flatten, LS, BLS)
    if bin_cadence_s > 0:
        lc_data = bin_lightcurve(lc_data, bin_cadence_s=bin_cadence_s)

    lc_data = flatten_lightcurve(lc_data, window_length=window_length)

    return lc_data
