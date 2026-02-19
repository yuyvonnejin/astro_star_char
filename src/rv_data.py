"""Phase 7: Radial velocity data retrieval and analysis.

Retrieves known planet parameters from the NASA Exoplanet Archive and
RV time series from the DACE (Data and Analysis Center for Exoplanets)
REST API. Provides RV periodogram analysis and detection limit estimation.
"""

import logging
import numpy as np

logger = logging.getLogger(__name__)

# NASA Exoplanet Archive TAP endpoint
NASA_EXOPLANET_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

# DACE REST API base URL (Geneva Observatory)
DACE_API_BASE = "https://dace-api.unige.ch/api"


def query_known_planets(star_name):
    """Query the NASA Exoplanet Archive for confirmed planets around a star.

    Parameters
    ----------
    star_name : str
        Host star name (e.g., 'HD 20794', 'tau Cet'). The query searches
        the hostname field in the Planetary Systems table.

    Returns
    -------
    list[dict]
        List of planet records, each with keys: pl_name, hostname,
        pl_orbper (days), pl_bmassj (Jupiter masses), pl_radj (Jupiter radii),
        discoverymethod, disc_year, pl_orbsmax (AU), pl_eqt (K).
        Returns empty list if no planets found or query fails.
    """
    import requests

    # Clean the name for ADQL query -- try both original and HD format
    search_names = _generate_search_names(star_name)

    for search_name in search_names:
        # ADQL query against the Planetary Systems composite table
        query = (
            "SELECT pl_name, hostname, pl_orbper, pl_bmassj, pl_bmasse, "
            "pl_radj, pl_rade, discoverymethod, disc_year, "
            "pl_orbsmax, pl_eqt, pl_orbeccen, sy_dist "
            "FROM ps "
            f"WHERE hostname = '{search_name}' "
            "AND default_flag = 1"
        )

        params = {
            "query": query,
            "format": "json",
        }

        logger.info("Querying NASA Exoplanet Archive for planets around '%s'", search_name)
        try:
            response = requests.get(NASA_EXOPLANET_TAP, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            if data:
                logger.info("Found %d confirmed planet(s) around '%s'", len(data), search_name)
                return _clean_planet_records(data)
        except requests.exceptions.Timeout:
            logger.warning("NASA Exoplanet Archive query timed out for '%s'", search_name)
        except requests.exceptions.RequestException as e:
            logger.warning("NASA Exoplanet Archive query failed for '%s': %s", search_name, e)
        except ValueError as e:
            logger.warning("Failed to parse NASA Exoplanet Archive response: %s", e)

    logger.info("No confirmed planets found for '%s'", star_name)
    return []


def query_dace_rv(star_name):
    """Retrieve public RV time series from the DACE archive via dace-query.

    DACE (Data and Analysis Center for Exoplanets) hosts reduced RV data
    from major spectrographs (HARPS, ESPRESSO, CORALIE, CORAVEL, etc.).
    Uses the dace-query Python package (v2.0.0) for access.

    Parameters
    ----------
    star_name : str
        Target name (e.g., 'HD 20794', 'tau Cet'). DACE typically expects
        HD numbers without spaces (e.g., 'HD20794').

    Returns
    -------
    dict or None
        RV data with keys:
        - time: array of RJD timestamps (reduced JD)
        - rv: array of radial velocities (m/s)
        - rv_err: array of RV uncertainties (m/s)
        - instrument: list of instrument names per measurement
        - n_measurements: total count
        - time_baseline_days: total time span
        - instruments: list of unique instruments
        - drs_qc: array of DRS quality control flags
        Returns None if no data found or query fails.
    """
    search_names = _generate_dace_names(star_name)

    for search_name in search_names:
        logger.info("Querying DACE (dace-query) for '%s'", search_name)
        try:
            from dace_query.spectroscopy import Spectroscopy

            data = Spectroscopy.get_timeseries(
                target=search_name,
                sorted_by_instrument=False,
                output_format='dict',
            )

            if data and 'rjd' in data and len(data['rjd']) > 0:
                return _parse_dace_query_rv(data, star_name)

        except ImportError:
            logger.warning("dace-query not installed; pip install dace-query")
            return None
        except Exception as e:
            logger.warning("DACE query failed for '%s': %s", search_name, e)

    logger.info("No DACE RV data found for '%s'", star_name)
    return None


def query_nasa_rv_data(star_name):
    """Query NASA Exoplanet Archive for RV orbital solution parameters.

    Note: The NASA Exoplanet Archive TAP does not host raw RV time series
    measurements. Individual RV time series must be obtained from instrument-
    specific archives (ESO/DACE for HARPS/ESPRESSO, Keck Archive for HIRES).

    This function queries the Planetary Systems table for orbital solution
    parameters (period, K amplitude, eccentricity) which indirectly characterize
    the RV signal.

    Parameters
    ----------
    star_name : str
        Host star name.

    Returns
    -------
    dict or None
        RV orbital solution summary, or None if not available.
    """
    import requests

    search_names = _generate_search_names(star_name)

    for search_name in search_names:
        query = (
            "SELECT pl_name, hostname, pl_orbper, pl_rvamp, pl_orbeccen, "
            "pl_bmasse, discoverymethod "
            "FROM ps "
            f"WHERE hostname = '{search_name}' "
            "AND discoverymethod = 'Radial Velocity' "
            "AND default_flag = 1"
        )

        params = {"query": query, "format": "json"}

        logger.info("Querying NASA for RV orbital solutions for '%s'", search_name)
        try:
            response = requests.get(NASA_EXOPLANET_TAP, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            if data and len(data) > 0:
                logger.info("Found %d RV orbital solution(s) for '%s'", len(data), search_name)
                return {
                    "source": "NASA_orbital_solutions",
                    "n_solutions": len(data),
                    "solutions": data,
                    "note": "These are orbital solutions, not raw RV time series. "
                            "Raw RV data requires DACE or ESO archive access.",
                }
        except Exception as e:
            logger.info("NASA RV orbital query for '%s': %s", search_name, e)

    return None


def rv_periodogram(time, rv, rv_err=None, min_period=1.0, max_period=None,
                   oversample_factor=5):
    """Compute Lomb-Scargle periodogram on RV time series.

    Parameters
    ----------
    time : array-like
        BJD timestamps.
    rv : array-like
        Radial velocities (m/s).
    rv_err : array-like, optional
        RV uncertainties (m/s).
    min_period : float
        Minimum period to search (days). Default 1.0.
    max_period : float, optional
        Maximum period to search (days). Default: half the time baseline.
    oversample_factor : int
        Frequency grid oversampling.

    Returns
    -------
    dict
        Periodogram results with keys: periods (array), power (array),
        best_period, best_power, fap, peaks (list of top 5 periods).
    """
    from astropy.timeseries import LombScargle

    time = np.asarray(time, dtype=float)
    rv = np.asarray(rv, dtype=float)

    baseline = time[-1] - time[0]
    if baseline <= 0:
        logger.error("RV time baseline is zero or negative")
        return _empty_periodogram()

    if max_period is None:
        max_period = baseline / 2.0

    max_period = min(max_period, baseline / 2.0)
    if max_period <= min_period:
        logger.warning("max_period (%.1f) <= min_period (%.1f)", max_period, min_period)
        return _empty_periodogram()

    logger.info("RV periodogram: %.1f-%.1f days, %d measurements, baseline %.1f days",
                min_period, max_period, len(time), baseline)

    if rv_err is not None:
        rv_err = np.asarray(rv_err, dtype=float)
        ls = LombScargle(time, rv, rv_err)
    else:
        ls = LombScargle(time, rv)

    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    frequency, power = ls.autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq,
        samples_per_peak=oversample_factor,
    )

    periods = 1.0 / frequency

    if len(power) == 0:
        return _empty_periodogram()

    # Find top peaks
    best_idx = np.argmax(power)
    best_period = float(periods[best_idx])
    best_power = float(power[best_idx])
    fap = float(ls.false_alarm_probability(best_power))

    # Find top 5 peaks (separated by at least 10% in period)
    peaks = _find_peaks(periods, power, n_peaks=5, min_separation_frac=0.1)

    logger.info("RV best period: %.2f days (power=%.4f, FAP=%.2e)", best_period, best_power, fap)

    return {
        "periods": periods,
        "power": power,
        "best_period": best_period,
        "best_power": best_power,
        "fap": fap,
        "peaks": peaks,
        "n_measurements": len(time),
        "time_baseline_days": float(baseline),
    }


def rv_detection_limit(time, rv_err, period_grid=None):
    """Estimate the RV semi-amplitude detection limit at each period.

    Uses the analytical approximation: K_min = SNR_threshold * sigma_eff / sqrt(N/2),
    where sigma_eff accounts for time sampling and period.

    Parameters
    ----------
    time : array-like
        BJD timestamps of observations.
    rv_err : array-like
        RV uncertainties (m/s) at each epoch.
    period_grid : array-like, optional
        Periods (days) at which to compute limits. Default: log-spaced 1-10000 days.

    Returns
    -------
    dict
        Detection limits with keys: periods (array), k_min_ms (array of minimum
        detectable semi-amplitudes in m/s), mass_min_mearth (array, if stellar
        mass provided), snr_threshold (float used).
    """
    time = np.asarray(time, dtype=float)
    rv_err = np.asarray(rv_err, dtype=float)

    if period_grid is None:
        period_grid = np.logspace(0, 4, 200)  # 1 to 10000 days

    period_grid = np.asarray(period_grid, dtype=float)
    n_obs = len(time)
    baseline = time[-1] - time[0]

    # Effective RV precision (weighted mean error)
    sigma_eff = np.sqrt(np.mean(rv_err**2))

    # SNR threshold for detection (typically 5-6 for significance)
    snr_threshold = 5.0

    # For each period, estimate the effective number of independent phase samples
    k_min = np.zeros_like(period_grid)
    for i, period in enumerate(period_grid):
        # Number of full cycles covered
        n_cycles = baseline / period
        if n_cycles < 1.5:
            # Fewer than 1.5 cycles: very poor sensitivity
            k_min[i] = np.inf
        else:
            # Effective N scales with total observations but with diminishing
            # returns for sparse sampling within each cycle
            n_eff = min(n_obs, n_cycles * 10)  # cap at ~10 points per cycle utility
            k_min[i] = snr_threshold * sigma_eff / np.sqrt(n_eff / 2.0)

    # Convert K to minimum planet mass (Earth masses) for a Sun-like star
    # K = (2*pi*G/P)^(1/3) * M_p*sin(i) / M_star^(2/3) / sqrt(1-e^2)
    # For circular orbit around 1 Msun: M_p [Mearth] ~ K [m/s] * P[days]^(1/3) * 3.16
    # This is approximate; exact conversion depends on stellar mass
    mass_min = k_min * (period_grid ** (1.0 / 3.0)) * 3.16  # very rough

    logger.info("RV detection limits: %.2f m/s median precision, %d obs, %.0f day baseline",
                sigma_eff, n_obs, baseline)

    return {
        "periods": period_grid,
        "k_min_ms": k_min,
        "mass_min_mearth_approx": mass_min,
        "sigma_eff_ms": float(sigma_eff),
        "n_obs": n_obs,
        "baseline_days": float(baseline),
        "snr_threshold": snr_threshold,
    }


def rv_to_planet_mass(k_ms, period_days, stellar_mass_msun=1.0, eccentricity=0.0):
    """Convert RV semi-amplitude to minimum planet mass.

    Parameters
    ----------
    k_ms : float or array
        RV semi-amplitude in m/s.
    period_days : float or array
        Orbital period in days.
    stellar_mass_msun : float
        Host star mass in solar masses.
    eccentricity : float
        Orbital eccentricity.

    Returns
    -------
    float or array
        Minimum planet mass (M*sin(i)) in Earth masses.
    """
    # Constants
    G = 6.674e-11  # m^3 kg^-1 s^-2
    M_sun = 1.989e30  # kg
    M_earth = 5.972e24  # kg

    period_s = np.asarray(period_days) * 86400.0
    k = np.asarray(k_ms)
    M_star = stellar_mass_msun * M_sun

    # M_p * sin(i) = K * (M_star)^(2/3) * (P/(2*pi*G))^(1/3) * sqrt(1-e^2)
    factor = (period_s / (2.0 * np.pi * G)) ** (1.0 / 3.0)
    m_planet_kg = k * (M_star ** (2.0 / 3.0)) * factor * np.sqrt(1.0 - eccentricity**2)

    return m_planet_kg / M_earth


# --- Internal helpers ---

def _generate_search_names(star_name):
    """Generate variant names for archive queries."""
    names = [star_name]

    # If it looks like "HD NNNNN", also try without space
    name_stripped = star_name.strip()
    if name_stripped.upper().startswith("HD "):
        names.append(name_stripped.upper())
        names.append(name_stripped.lower())

    # Common aliases
    aliases = {
        "82 G. Eridani": ["HD 20794", "e Eri", "82 Eridani"],
        "Tau Ceti": ["HD 10700", "tau Cet", "HIP 8102"],
        "Alpha Centauri A": ["HD 128620", "alf Cen A", "HIP 71683"],
        "Eta Cassiopeiae A": ["HD 4614", "eta Cas A", "HIP 3821"],
        "Delta Pavonis": ["HD 190248", "del Pav", "HIP 99240"],
        "61 Virginis": ["HD 115617", "61 Vir", "HIP 64924"],
        "Beta CVn": ["HD 109358", "bet CVn", "Chara", "HIP 61317"],
        "Zeta Tucanae": ["HD 1581", "zet Tuc", "HIP 1599"],
        "18 Scorpii": ["HD 146233", "18 Sco", "HIP 79672"],
        "HD 134060": ["HIP 74273"],
    }

    for key, alias_list in aliases.items():
        if star_name.lower() in key.lower() or key.lower() in star_name.lower():
            for alias in alias_list:
                if alias not in names:
                    names.append(alias)

    return names


def _generate_dace_names(star_name):
    """Generate DACE-specific name variants.

    DACE typically uses HD numbers without spaces (e.g., 'HD20794').
    """
    names = []
    name_stripped = star_name.strip()

    # If already HD format, try with and without space
    if name_stripped.upper().startswith("HD "):
        hd_nospace = name_stripped.upper().replace("HD ", "HD")
        names.append(hd_nospace)
        names.append(name_stripped)
    elif name_stripped.upper().startswith("HD"):
        names.append(name_stripped)
        # Try with space
        hd_space = "HD " + name_stripped[2:]
        names.append(hd_space)

    # Look up HD number from aliases
    hd_map = {
        "82 G. Eridani": "HD20794",
        "Tau Ceti": "HD10700",
        "Alpha Centauri A": "HD128620",
        "Eta Cassiopeiae A": "HD4614",
        "Delta Pavonis": "HD190248",
        "61 Virginis": "HD115617",
        "Beta CVn": "HD109358",
        "Zeta Tucanae": "HD1581",
        "18 Scorpii": "HD146233",
        "HD 134060": "HD134060",
    }

    for key, hd_name in hd_map.items():
        if name_stripped.lower() in key.lower() or key.lower() in name_stripped.lower():
            if hd_name not in names:
                names.insert(0, hd_name)  # HD without space first

    # If nothing matched, just try the name as-is
    if not names:
        names.append(name_stripped)

    return names


def _parse_dace_query_rv(data, star_name):
    """Parse dace-query get_timeseries() output into standard RV data format.

    The dace-query package returns a flat dict with parallel arrays when
    sorted_by_instrument=False.
    """
    times = np.array(data.get('rjd', []), dtype=float)
    rvs = np.array(data.get('rv', []), dtype=float)
    rv_errs = np.array(data.get('rv_err', []), dtype=float)
    instruments = list(data.get('ins_name', ['unknown'] * len(times)))
    drs_qc = data.get('drs_qc', [True] * len(times))

    if len(times) == 0:
        return None

    # Filter: valid RV, valid error, positive error
    valid = np.isfinite(rvs) & np.isfinite(rv_errs) & (rv_errs > 0)

    # Apply DRS quality control if available
    if drs_qc is not None:
        qc = np.array(drs_qc, dtype=bool)
        valid = valid & qc

    times = times[valid]
    rvs = rvs[valid]
    rv_errs = rv_errs[valid]
    instruments = [str(instruments[i]) for i in range(len(valid)) if valid[i]]

    if len(times) == 0:
        logger.info("DACE returned data but all measurements filtered for '%s'", star_name)
        return None

    # Sort by time
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    rvs = rvs[sort_idx]
    rv_errs = rv_errs[sort_idx]
    instruments = [instruments[i] for i in sort_idx]

    unique_instruments = sorted(set(instruments))
    baseline = float(times[-1] - times[0])

    # Per-instrument summary
    inst_summary = {}
    for inst in unique_instruments:
        mask = np.array([i == inst for i in instruments])
        inst_summary[inst] = {
            "n_measurements": int(np.sum(mask)),
            "median_err_ms": round(float(np.median(rv_errs[mask])), 3),
            "time_span_days": round(float(times[mask][-1] - times[mask][0]), 1),
        }

    logger.info("DACE RV for '%s': %d measurements (QC-filtered), %.0f day baseline, "
                "instruments: %s",
                star_name, len(times), baseline, ", ".join(unique_instruments))
    for inst, summary in inst_summary.items():
        logger.info("  %s: %d meas, median err %.3f m/s, %.0f day span",
                     inst, summary["n_measurements"], summary["median_err_ms"],
                     summary["time_span_days"])

    return {
        "time": times,
        "rv": rvs,
        "rv_err": rv_errs,
        "instrument": instruments,
        "n_measurements": len(times),
        "time_baseline_days": baseline,
        "instruments": unique_instruments,
        "instrument_summary": inst_summary,
    }


def _clean_planet_records(data):
    """Clean planet records from NASA Exoplanet Archive response."""
    planets = []
    for row in data:
        planet = {
            "pl_name": row.get("pl_name"),
            "hostname": row.get("hostname"),
            "period_days": row.get("pl_orbper"),
            "mass_jupiter": row.get("pl_bmassj"),
            "mass_earth": row.get("pl_bmasse"),
            "radius_jupiter": row.get("pl_radj"),
            "radius_earth": row.get("pl_rade"),
            "discovery_method": row.get("discoverymethod"),
            "discovery_year": row.get("disc_year"),
            "semi_major_axis_au": row.get("pl_orbsmax"),
            "eq_temp_k": row.get("pl_eqt"),
            "eccentricity": row.get("pl_orbeccen"),
            "distance_pc": row.get("sy_dist"),
        }
        planets.append(planet)
    return planets


def _parse_dace_rv(data, star_name):
    """Parse DACE API response into standard RV data format."""
    times = []
    rvs = []
    rv_errs = []
    instruments = []

    # DACE response format varies; handle common structures
    if isinstance(data, dict):
        # Nested by instrument
        for inst_name, measurements in data.items():
            if isinstance(measurements, list):
                for m in measurements:
                    t = m.get("rjd") or m.get("bjd") or m.get("jdb")
                    rv = m.get("rv") or m.get("radvel")
                    err = m.get("rv_err") or m.get("svrad") or m.get("radvel_err")
                    if t is not None and rv is not None:
                        times.append(float(t))
                        rvs.append(float(rv))
                        rv_errs.append(float(err) if err is not None else 1.0)
                        instruments.append(str(inst_name))
    elif isinstance(data, list):
        for m in data:
            t = m.get("rjd") or m.get("bjd") or m.get("jdb")
            rv = m.get("rv") or m.get("radvel")
            err = m.get("rv_err") or m.get("svrad") or m.get("radvel_err")
            inst = m.get("instrument") or m.get("ins_name") or "unknown"
            if t is not None and rv is not None:
                times.append(float(t))
                rvs.append(float(rv))
                rv_errs.append(float(err) if err is not None else 1.0)
                instruments.append(str(inst))

    if not times:
        logger.info("DACE returned data but no parseable RV measurements for '%s'", star_name)
        return None

    times = np.array(times)
    rvs = np.array(rvs)
    rv_errs = np.array(rv_errs)

    # Sort by time
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    rvs = rvs[sort_idx]
    rv_errs = rv_errs[sort_idx]
    instruments = [instruments[i] for i in sort_idx]

    unique_instruments = sorted(set(instruments))
    baseline = float(times[-1] - times[0])

    logger.info("DACE RV for '%s': %d measurements, %.0f day baseline, instruments: %s",
                star_name, len(times), baseline, ", ".join(unique_instruments))

    return {
        "time": times,
        "rv": rvs,
        "rv_err": rv_errs,
        "instrument": instruments,
        "n_measurements": len(times),
        "time_baseline_days": baseline,
        "instruments": unique_instruments,
    }


def _parse_nasa_rv(data, star_name):
    """Parse NASA Exoplanet Archive RV time series response."""
    times = []
    rvs = []
    rv_errs = []
    instruments = []

    for row in data:
        t = row.get("rv_time")
        rv = row.get("rv_value")
        err = row.get("rv_err")
        inst = row.get("rv_instrument", "unknown")
        if t is not None and rv is not None:
            times.append(float(t))
            rvs.append(float(rv))
            rv_errs.append(float(err) if err is not None else 1.0)
            instruments.append(str(inst))

    if not times:
        return None

    times = np.array(times)
    rvs = np.array(rvs)
    rv_errs = np.array(rv_errs)

    sort_idx = np.argsort(times)
    times = times[sort_idx]
    rvs = rvs[sort_idx]
    rv_errs = rv_errs[sort_idx]
    instruments = [instruments[i] for i in sort_idx]

    unique_instruments = sorted(set(instruments))
    baseline = float(times[-1] - times[0])

    logger.info("NASA RV for '%s': %d measurements, instruments: %s",
                star_name, len(times), ", ".join(unique_instruments))

    return {
        "time": times,
        "rv": rvs,
        "rv_err": rv_errs,
        "instrument": instruments,
        "n_measurements": len(times),
        "time_baseline_days": baseline,
        "instruments": unique_instruments,
    }


def _find_peaks(periods, power, n_peaks=5, min_separation_frac=0.1):
    """Find top N peaks in a periodogram, separated by at least min_separation_frac."""
    sorted_idx = np.argsort(power)[::-1]
    peaks = []

    for idx in sorted_idx:
        p = float(periods[idx])
        pw = float(power[idx])

        # Check separation from already-found peaks
        too_close = False
        for existing_p, _ in peaks:
            if abs(p - existing_p) / existing_p < min_separation_frac:
                too_close = True
                break

        if not too_close:
            peaks.append((p, pw))
            if len(peaks) >= n_peaks:
                break

    return [{"period_days": p, "power": pw} for p, pw in peaks]


def _empty_periodogram():
    """Return empty periodogram result."""
    return {
        "periods": np.array([]),
        "power": np.array([]),
        "best_period": None,
        "best_power": None,
        "fap": None,
        "peaks": [],
        "n_measurements": 0,
        "time_baseline_days": 0,
    }
