"""Phase 7b.2: RadVel-based joint Keplerian fitting for multi-planet RV systems.

Replaces the sinusoidal subtraction approach from Phase 7b with a proper
joint Keplerian + per-instrument offset + jitter model via RadVel.

Key improvements over rv_data.rv_residual_analysis():
- Keplerian orbits (not sinusoids): models eccentricity correctly
- Joint fit: all planets + offsets + jitter fitted simultaneously
- Per-instrument gamma (zero-point) and jitter parameters
- MAP optimization via scipy.optimize (RadVel wrapper)

References:
    Fulton et al. (2018), PASP 130, 044504 -- RadVel paper
    Nari et al. (2025), A&A 693, A297 -- HD 20794 3-planet solution
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def fit_keplerian(time, rv, rv_err, instruments, planet_params,
                  exclude_instruments=None, run_mcmc=False,
                  mcmc_nwalkers=50, mcmc_nsteps=1000):
    """Fit a joint Keplerian model to multi-instrument RV data.

    Constructs a RadVel model with N Keplerian orbits plus per-instrument
    gamma (offset) and jitter parameters, then optimizes via MAP.

    Parameters
    ----------
    time : array-like
        BJD or RJD timestamps.
    rv : array-like
        Radial velocities (m/s).
    rv_err : array-like
        RV uncertainties (m/s).
    instruments : list[str]
        Instrument name per measurement.
    planet_params : list[dict]
        Initial guesses for each planet. Each dict has keys:
        period, tc, e, w, k (all floats).
        - period: orbital period (days)
        - tc: time of conjunction (same time system as `time`)
        - e: eccentricity (0 to <1)
        - w: argument of periastron (radians)
        - k: RV semi-amplitude (m/s)
    exclude_instruments : list[str], optional
        Instruments to exclude (e.g., ['CORAVEL-S']).
    run_mcmc : bool
        If True, run MCMC after MAP for uncertainties.
    mcmc_nwalkers : int
        Number of MCMC walkers.
    mcmc_nsteps : int
        Number of MCMC steps per walker.

    Returns
    -------
    dict
        Fit results with keys:
        - planets: list of dicts with period, tc, e, w, k, k_err
        - instruments: dict of {name: {gamma, jitter}}
        - residuals: np.array of residual RVs (m/s)
        - model: np.array of model RVs (m/s)
        - rms_before_ms: RMS of input data (m/s)
        - rms_after_ms: RMS of residuals (m/s)
        - n_measurements: int
        - excluded_instruments: list
        - method: str
        - status: str
    """
    import radvel

    time = np.asarray(time, dtype=float)
    rv = np.asarray(rv, dtype=float)
    rv_err = np.asarray(rv_err, dtype=float)
    instruments = list(instruments)

    if exclude_instruments is None:
        exclude_instruments = []

    # Filter out excluded instruments
    if exclude_instruments:
        inst_arr = np.array(instruments)
        keep = np.ones(len(time), dtype=bool)
        for exc in exclude_instruments:
            keep &= (inst_arr != exc)
        n_removed = int(np.sum(~keep))
        time = time[keep]
        rv = rv[keep]
        rv_err = rv_err[keep]
        instruments = [instruments[i] for i in range(len(keep)) if keep[i]]
        logger.info("Excluded %d measurements from %s",
                     n_removed, exclude_instruments)

    if len(time) == 0:
        logger.error("No measurements remain after exclusion")
        return {"status": "error", "error": "No measurements after exclusion"}

    n_planets = len(planet_params)
    unique_instruments = sorted(set(instruments))
    inst_arr = np.array(instruments)

    rms_before = float(np.std(rv))
    logger.info("Keplerian fit: %d planets, %d measurements, %d instruments, "
                "input RMS=%.2f m/s",
                n_planets, len(time), len(unique_instruments), rms_before)

    # Build RadVel parameter set
    params = radvel.Parameters(n_planets, basis='per tc e w k')

    for i, pp in enumerate(planet_params):
        idx = i + 1  # RadVel uses 1-based indexing
        params[f'per{idx}'] = radvel.Parameter(value=float(pp['period']))
        params[f'tc{idx}'] = radvel.Parameter(value=float(pp['tc']))
        params[f'e{idx}'] = radvel.Parameter(value=float(pp['e']))
        params[f'w{idx}'] = radvel.Parameter(value=float(pp['w']))
        params[f'k{idx}'] = radvel.Parameter(value=float(pp['k']))

    # Trend parameters (no trend)
    params['dvdt'] = radvel.Parameter(value=0.0)
    params['curv'] = radvel.Parameter(value=0.0)

    # Create RV model
    time_base = float(np.median(time))
    mod = radvel.RVModel(params, time_base=time_base)

    # Create per-instrument likelihoods
    likelihoods = []
    for inst in unique_instruments:
        mask = inst_arr == inst
        inst_time = time[mask]
        inst_rv = rv[mask]
        inst_rv_err = rv_err[mask]

        like = radvel.likelihood.RVLikelihood(
            mod, inst_time, inst_rv, inst_rv_err,
            suffix=f'_{inst}',
        )
        # Set initial gamma to instrument median
        like.params[f'gamma_{inst}'] = radvel.Parameter(
            value=float(np.median(inst_rv))
        )
        # Set initial jitter to instrument median error
        like.params[f'jit_{inst}'] = radvel.Parameter(
            value=float(np.median(inst_rv_err))
        )
        likelihoods.append(like)
        logger.info("  Instrument %s: %d pts, median RV=%.1f m/s, "
                     "median err=%.3f m/s",
                     inst, int(np.sum(mask)),
                     float(np.median(inst_rv)),
                     float(np.median(inst_rv_err)))

    # Composite likelihood
    like_composite = radvel.likelihood.CompositeLikelihood(likelihoods)

    # Build posterior with priors
    post = radvel.posterior.Posterior(like_composite)

    # Add eccentricity priors (keep e in [0, 1))
    for i in range(n_planets):
        idx = i + 1
        post.priors.append(
            radvel.prior.EccentricityPrior(n_planets)
        )
        break  # EccentricityPrior covers all planets at once

    # Fix trend to zero (no long-term acceleration)
    post.params['dvdt'].vary = False
    post.params['curv'].vary = False

    # Fix periods (well-constrained from literature) to help convergence
    for i in range(n_planets):
        idx = i + 1
        post.params[f'per{idx}'].vary = False

    # MAP optimization
    logger.info("Running MAP optimization...")
    try:
        post = radvel.fitting.maxlike_fitting(post, verbose=False)
    except Exception as e:
        logger.error("MAP fitting failed: %s", e)
        return {"status": "error", "error": f"MAP fitting failed: {e}"}

    logger.info("MAP optimization complete")

    # Now allow periods to vary and re-optimize
    for i in range(n_planets):
        idx = i + 1
        post.params[f'per{idx}'].vary = True

    try:
        post = radvel.fitting.maxlike_fitting(post, verbose=False)
    except Exception as e:
        logger.warning("Second MAP pass (free periods) failed: %s; "
                       "using fixed-period solution", e)

    # Optional MCMC
    mcmc_chains = None
    if run_mcmc:
        logger.info("Running MCMC (%d walkers x %d steps)...",
                     mcmc_nwalkers, mcmc_nsteps)
        try:
            chains = radvel.mcmc(
                post, nwalkers=mcmc_nwalkers, nrun=mcmc_nsteps,
                savename=None,
            )
            mcmc_chains = chains
            logger.info("MCMC complete")
        except Exception as e:
            logger.warning("MCMC failed: %s", e)

    # Extract results
    fitted_params = post.params

    # Planet parameters
    planets_result = []
    for i in range(n_planets):
        idx = i + 1
        k_val = float(fitted_params[f'k{idx}'].value)
        k_err = None
        if mcmc_chains is not None:
            k_samples = mcmc_chains[f'k{idx}']
            k_err = float(np.std(k_samples))

        planets_result.append({
            "period": float(fitted_params[f'per{idx}'].value),
            "tc": float(fitted_params[f'tc{idx}'].value),
            "e": float(fitted_params[f'e{idx}'].value),
            "w": float(fitted_params[f'w{idx}'].value),
            "k": k_val,
            "k_err": k_err,
        })
        logger.info("  Planet %d: P=%.3f d, K=%.4f m/s, e=%.4f",
                     idx, planets_result[-1]["period"],
                     k_val, planets_result[-1]["e"])

    # Instrument parameters
    instruments_result = {}
    for inst in unique_instruments:
        gamma = float(fitted_params[f'gamma_{inst}'].value)
        jit = float(fitted_params[f'jit_{inst}'].value)
        instruments_result[inst] = {"gamma": gamma, "jitter": jit}
        logger.info("  %s: gamma=%.3f m/s, jitter=%.3f m/s",
                     inst, gamma, jit)

    # Compute model and residuals
    model_rv = np.zeros_like(time)
    residuals = np.zeros_like(time)

    for inst in unique_instruments:
        mask = inst_arr == inst
        inst_time = time[mask]
        # Model = Keplerian + gamma
        model_at_inst = mod(inst_time) + instruments_result[inst]["gamma"]
        model_rv[mask] = model_at_inst
        residuals[mask] = rv[mask] - model_at_inst

    rms_after = float(np.std(residuals))
    logger.info("Keplerian fit: RMS %.2f -> %.2f m/s", rms_before, rms_after)

    # Log gamma differences (instrument offsets)
    if len(unique_instruments) > 1:
        gammas = [instruments_result[inst]["gamma"]
                  for inst in unique_instruments]
        for i, inst_i in enumerate(unique_instruments):
            for j, inst_j in enumerate(unique_instruments):
                if j > i:
                    diff = gammas[j] - gammas[i]
                    logger.info("  Offset %s - %s = %.3f m/s",
                                 inst_j, inst_i, diff)

    return {
        "planets": planets_result,
        "instruments": instruments_result,
        "residuals": residuals,
        "model": model_rv,
        "rms_before_ms": round(rms_before, 4),
        "rms_after_ms": round(rms_after, 4),
        "n_measurements": len(time),
        "excluded_instruments": list(exclude_instruments),
        "method": "radvel_keplerian",
        "status": "ok",
    }


def keplerian_residual_analysis(time, rv, rv_err, instruments, planet_params,
                                 exclude_instruments=None,
                                 min_period=1.0, max_period=None,
                                 run_mcmc=False):
    """Full Keplerian residual analysis: fit + periodogram on residuals.

    Drop-in replacement for rv_data.rv_residual_analysis() that uses
    a joint Keplerian model instead of sinusoidal subtraction.

    Parameters
    ----------
    time : array-like
        BJD/RJD timestamps.
    rv : array-like
        Radial velocities (m/s).
    rv_err : array-like
        RV uncertainties (m/s).
    instruments : list[str]
        Instrument name per measurement.
    planet_params : list[dict]
        Initial guesses per planet (period, tc, e, w, k).
    exclude_instruments : list[str], optional
        Instruments to exclude.
    min_period : float
        Minimum period for residual periodogram (days).
    max_period : float, optional
        Maximum period for residual periodogram (days).
    run_mcmc : bool
        Run MCMC for uncertainties.

    Returns
    -------
    dict
        Backward-compatible with rv_residual_analysis() output, plus
        additional 'keplerian_fit' key with full RadVel results.
        Keys: known_periods_used, offset_correction, sinusoid_subtraction,
        original_periodogram, residual_periodogram, keplerian_fit.
    """
    from src.rv_data import rv_periodogram

    time = np.asarray(time, dtype=float)
    rv = np.asarray(rv, dtype=float)
    rv_err = np.asarray(rv_err, dtype=float)

    known_periods = [pp["period"] for pp in planet_params]
    result = {"known_periods_used": known_periods}

    # Original periodogram for comparison
    logger.info("Keplerian residual analysis: computing original periodogram")
    orig_pg = rv_periodogram(time, rv, rv_err, min_period=min_period,
                              max_period=max_period)
    result["original_periodogram"] = orig_pg

    # Run Keplerian fit
    kep_result = fit_keplerian(
        time, rv, rv_err, instruments, planet_params,
        exclude_instruments=exclude_instruments,
        run_mcmc=run_mcmc,
    )

    result["keplerian_fit"] = kep_result

    if kep_result["status"] != "ok":
        logger.error("Keplerian fit failed: %s", kep_result.get("error"))
        result["offset_correction"] = None
        result["sinusoid_subtraction"] = None
        result["residual_periodogram"] = rv_periodogram(
            time, rv, rv_err, min_period=min_period, max_period=max_period
        )
        return result

    # Map Keplerian results to backward-compatible format
    result["offset_correction"] = {
        "offsets": {inst: data["gamma"]
                    for inst, data in kep_result["instruments"].items()},
        "n_instruments": len(kep_result["instruments"]),
        "method": "radvel_joint_fit",
    }

    result["sinusoid_subtraction"] = {
        "fitted_components": [
            {
                "period_days": p["period"],
                "amplitude_ms": abs(p["k"]),
                "eccentricity": p["e"],
                "method": "keplerian",
            }
            for p in kep_result["planets"]
        ],
        "rms_before_ms": kep_result["rms_before_ms"],
        "rms_after_ms": kep_result["rms_after_ms"],
    }

    # Periodogram on residuals
    # Use only the non-excluded data for residual periodogram
    residuals = kep_result["residuals"]
    n_total = len(time)
    n_fit = kep_result["n_measurements"]

    if n_fit < n_total and exclude_instruments:
        # Recompute mask for time/rv_err alignment
        inst_arr = np.array(list(instruments))
        keep = np.ones(n_total, dtype=bool)
        for exc in exclude_instruments:
            keep &= (inst_arr != exc)
        time_fit = time[keep]
        rv_err_fit = rv_err[keep]
    else:
        time_fit = time
        rv_err_fit = rv_err

    logger.info("Keplerian residual analysis: computing residual periodogram")
    resid_pg = rv_periodogram(time_fit, residuals, rv_err_fit,
                               min_period=min_period, max_period=max_period)
    result["residual_periodogram"] = resid_pg

    logger.info("Keplerian residual analysis complete: original best P=%.2f d, "
                "residual best P=%.2f d, RMS %.2f -> %.2f m/s",
                orig_pg.get("best_period") or 0,
                resid_pg.get("best_period") or 0,
                kep_result["rms_before_ms"],
                kep_result["rms_after_ms"])

    return result
