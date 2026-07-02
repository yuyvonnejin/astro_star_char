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


def _chain_err(chains, key):
    """Half the 16th-84th percentile spread of an MCMC chain column.

    Equivalent to 1-sigma for a Gaussian posterior. Returns None when
    chains are unavailable or the parameter was not varied.
    """
    if chains is None or key not in chains:
        return None
    lo, hi = np.percentile(np.asarray(chains[key], dtype=float),
                           [15.865, 84.135])
    return float((hi - lo) / 2.0)


def _chain_med(chains, key, default):
    """Posterior median of a chain column, or `default` without chains.

    radvel.mcmc leaves the posterior object at an arbitrary final
    walker position, so post-MCMC parameter values must be taken from
    the chains, not from post.params.
    """
    if chains is None or key not in chains:
        return default
    return float(np.median(np.asarray(chains[key], dtype=float)))


def _sync_vary_flags(post):
    """Push params-dict vary flags into the posterior's vector cache.

    RadVel >= 1.4 stores parameter values and vary flags in a numpy
    vector built at Posterior construction, and caches the vary-index
    list on first use. Setting params[...].vary afterward changes only
    the dict; the optimizer keeps using the stale vector. Without this
    sync every parameter stays free regardless of vary flags.
    """
    post.vector.dict_to_vector()
    post.list_vary_params()


# Default quasi-periodic GP hyperparameters, tuned for a solar-type
# magnetic activity cycle (HD 20794: ~3000 d, Nari et al. 2025).
# Values: initial guess; bounds: HardBounds prior ranges.
DEFAULT_GP_HYPERPARAMS = {
    "gp_amp": {"value": 2.0, "bounds": (0.01, 20.0)},
    "gp_explength": {"value": 3000.0, "bounds": (200.0, 30000.0)},
    "gp_per": {"value": 3000.0, "bounds": (500.0, 10000.0)},
    "gp_perlength": {"value": 0.5, "bounds": (0.05, 1.5)},
}

_GP_HNAMES = ["gp_per", "gp_perlength", "gp_explength", "gp_amp"]


def fit_keplerian(time, rv, rv_err, instruments, planet_params,
                  exclude_instruments=None, fix_eccentricities=False,
                  run_mcmc=False, mcmc_nwalkers=50, mcmc_nsteps=1000,
                  mcmc_ensembles=3,
                  use_gp=False, gp_hyperparams=None):
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
    fix_eccentricities : bool
        If True, keep eccentricities fixed at initial values throughout
        the fit. Recommended for sub-m/s signals where the data cannot
        independently constrain e (e.g., HD 20794). Default: False.
    run_mcmc : bool
        If True, run MCMC after MAP for uncertainties.
    mcmc_nwalkers : int
        Number of MCMC walkers.
    mcmc_nsteps : int
        Number of MCMC steps per walker.
    mcmc_ensembles : int
        Number of independent ensembles. radvel's default is 8, and
        with serial=True they run SEQUENTIALLY -- 8x the work. 3 is
        the minimum for cross-ensemble Gelman-Rubin diagnostics.
    use_gp : bool
        If True, model correlated stellar activity noise with a
        quasi-periodic Gaussian Process (RadVel GPLikelihood, pure
        numpy QuasiPer kernel, O(N^3) -- use with nightly-binned
        data). Hyperparameters are shared across instruments.
        Default hyperparameters target a solar-type magnetic cycle.
    gp_hyperparams : dict, optional
        Override DEFAULT_GP_HYPERPARAMS. Keys: gp_amp, gp_explength,
        gp_per, gp_perlength; values: {"value": float,
        "bounds": (lo, hi)}.

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

    # Pre-center: subtract per-instrument median RV so RadVel works
    # near zero instead of ~88 km/s. This avoids numerical precision
    # issues when fitting sub-m/s signals on km/s baselines.
    inst_medians = {}
    rv_centered = rv.copy()
    for inst in unique_instruments:
        mask = inst_arr == inst
        med = float(np.median(rv[mask]))
        inst_medians[inst] = med
        rv_centered[mask] -= med

    logger.info("Keplerian fit: %d planets, %d measurements, %d instruments, "
                "input RMS=%.2f m/s (centered RMS=%.2f m/s)",
                n_planets, len(time), len(unique_instruments),
                rms_before, float(np.std(rv_centered)))

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

    # GP hyperparameters (shared across instruments) must be in the
    # params dict before the model/likelihoods build their vector
    gp_config = {}
    if use_gp:
        gp_config = {k: dict(v) for k, v in DEFAULT_GP_HYPERPARAMS.items()}
        if gp_hyperparams:
            for k, v in gp_hyperparams.items():
                gp_config.setdefault(k, {}).update(v)
        for name in _GP_HNAMES:
            params[name] = radvel.Parameter(value=float(gp_config[name]["value"]))

    # Create RV model
    time_base = float(np.median(time))
    mod = radvel.RVModel(params, time_base=time_base)

    # Create per-instrument likelihoods on centered data
    likelihoods = []
    like_by_inst = {}
    for inst in unique_instruments:
        mask = inst_arr == inst
        inst_time = time[mask]
        inst_rv = rv_centered[mask]
        inst_rv_err = rv_err[mask]

        if use_gp:
            like = radvel.likelihood.GPLikelihood(
                mod, inst_time, inst_rv, inst_rv_err,
                hnames=_GP_HNAMES,
                suffix=f'_{inst}',
                kernel_name="QuasiPer",
            )
        else:
            like = radvel.likelihood.RVLikelihood(
                mod, inst_time, inst_rv, inst_rv_err,
                suffix=f'_{inst}',
            )
        like_by_inst[inst] = like
        # Gamma starts at 0 (data is pre-centered)
        like.params[f'gamma_{inst}'] = radvel.Parameter(value=0.0)
        # Set initial jitter to instrument median error
        like.params[f'jit_{inst}'] = radvel.Parameter(
            value=float(np.median(inst_rv_err))
        )
        likelihoods.append(like)
        logger.info("  Instrument %s: %d pts, median RV=%.1f m/s "
                     "(raw median=%.1f), median err=%.3f m/s",
                     inst, int(np.sum(mask)),
                     float(np.median(inst_rv)),
                     inst_medians[inst],
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

    # GP hyperparameter priors (hard bounds keep the GP from
    # absorbing planetary signals or collapsing to white noise)
    if use_gp:
        for name in _GP_HNAMES:
            lo, hi = gp_config[name]["bounds"]
            post.priors.append(radvel.prior.HardBounds(name, lo, hi))

    # Fix trend to zero (no long-term acceleration)
    post.params['dvdt'].vary = False
    post.params['curv'].vary = False

    # Pass 1: Fix periods AND eccentricities to help convergence.
    # With sub-m/s signals buried in km/s offsets, the optimizer must
    # first find the correct gamma (offset) values before it can
    # meaningfully fit orbital parameters.
    for i in range(n_planets):
        idx = i + 1
        post.params[f'per{idx}'].vary = False
        post.params[f'e{idx}'].vary = False
        post.params[f'w{idx}'].vary = False
    if use_gp:
        # Pass 1: only the GP amplitude floats; timescales wait
        # until offsets and K amplitudes are roughly right
        post.params['gp_amp'].vary = True
        post.params['gp_per'].vary = False
        post.params['gp_explength'].vary = False
        post.params['gp_perlength'].vary = False
    _sync_vary_flags(post)

    logger.info("Running MAP optimization (pass 1: fixed P, e, w)...")
    try:
        post = radvel.fitting.maxlike_fitting(post, verbose=False)
    except Exception as e:
        logger.error("MAP fitting pass 1 failed: %s", e)
        return {"status": "error", "error": f"MAP fitting failed: {e}"}

    # Pass 2: Free eccentricities and w (unless fix_eccentricities),
    # keep periods fixed
    if not fix_eccentricities:
        for i in range(n_planets):
            idx = i + 1
            post.params[f'e{idx}'].vary = True
            post.params[f'w{idx}'].vary = True
        logger.info("Running MAP optimization (pass 2: free e, w)...")
    else:
        # Even with fixed e, free w for phase adjustment
        for i in range(n_planets):
            idx = i + 1
            post.params[f'w{idx}'].vary = True
        logger.info("Running MAP optimization (pass 2: free w, e fixed)...")
    if use_gp:
        # Free GP timescales once orbital phases are roughly set
        post.params['gp_per'].vary = True
        post.params['gp_explength'].vary = True
        post.params['gp_perlength'].vary = True
    _sync_vary_flags(post)

    try:
        post = radvel.fitting.maxlike_fitting(post, verbose=False)
    except Exception as e:
        logger.warning("MAP pass 2 failed: %s; using pass 1 solution", e)

    # Pass 3: Free periods for final refinement
    for i in range(n_planets):
        idx = i + 1
        post.params[f'per{idx}'].vary = True
    _sync_vary_flags(post)

    logger.info("Running MAP optimization (pass 3: free P)...")
    try:
        post = radvel.fitting.maxlike_fitting(post, verbose=False)
    except Exception as e:
        logger.warning("MAP pass 3 (free periods) failed: %s; "
                       "using pass 2 solution", e)

    logger.info("MAP optimization complete")

    # Optional MCMC
    mcmc_chains = None
    if run_mcmc:
        logger.info("Running MCMC (%d ensembles x %d walkers x %d "
                     "steps, serial)...",
                     mcmc_ensembles, mcmc_nwalkers, mcmc_nsteps)
        try:
            # serial=True: radvel's default parallel ensembles use
            # multiprocessing, which fails when called from library
            # code on Windows (spawn re-imports without __main__).
            # ensembles must be passed explicitly: serial mode runs
            # them sequentially, so the default of 8 means 8x work.
            chains = radvel.mcmc(
                post, nwalkers=mcmc_nwalkers, nrun=mcmc_nsteps,
                ensembles=mcmc_ensembles,
                serial=True, headless=True, savename=None,
            )
            mcmc_chains = chains
            logger.info("MCMC complete: %d samples", len(chains))
        except Exception as e:
            logger.warning("MCMC failed: %s", e)

    if mcmc_chains is not None:
        # radvel.mcmc leaves post at an arbitrary final walker
        # position. Restore every sampled parameter to its posterior
        # median so reported values, residuals, and GP predictions
        # are all computed from the same consistent state.
        n_restored = 0
        for key in post.name_vary_params():
            if key in mcmc_chains:
                post.params[key].value = _chain_med(mcmc_chains, key,
                                                    post.params[key].value)
                n_restored += 1
        post.vector.dict_to_vector()
        logger.info("Posterior restored to chain medians (%d params)",
                    n_restored)

    # Extract results
    fitted_params = post.params

    # Planet parameters. With MCMC, post was restored to posterior
    # medians above; without, these are the MAP values.
    planets_result = []
    for i in range(n_planets):
        idx = i + 1
        k_val = float(fitted_params[f'k{idx}'].value)

        planets_result.append({
            "period": float(fitted_params[f'per{idx}'].value),
            "tc": float(fitted_params[f'tc{idx}'].value),
            "e": float(fitted_params[f'e{idx}'].value),
            "w": float(fitted_params[f'w{idx}'].value),
            "k": k_val,
            "k_err": _chain_err(mcmc_chains, f'k{idx}'),
            "period_err": _chain_err(mcmc_chains, f'per{idx}'),
            "tc_err": _chain_err(mcmc_chains, f'tc{idx}'),
            "e_err": _chain_err(mcmc_chains, f'e{idx}'),
            "w_err": _chain_err(mcmc_chains, f'w{idx}'),
        })
        logger.info("  Planet %d: P=%.3f d, K=%.4f m/s, e=%.4f",
                     idx, planets_result[-1]["period"],
                     k_val, planets_result[-1]["e"])

    # Instrument parameters -- add back the pre-subtracted medians
    # so reported gammas are in the original absolute RV frame
    instruments_result = {}
    for inst in unique_instruments:
        gamma_centered = float(fitted_params[f'gamma_{inst}'].value)
        gamma_absolute = gamma_centered + inst_medians[inst]
        jit = float(fitted_params[f'jit_{inst}'].value)
        # Clamp tiny negative jitter at optimizer boundary to zero
        if jit < 0 and jit > -1e-6:
            jit = 0.0
        instruments_result[inst] = {
            "gamma": gamma_absolute,
            "gamma_centered": gamma_centered,
            "jitter": jit,
            "gamma_err": _chain_err(mcmc_chains, f'gamma_{inst}'),
            "jitter_err": _chain_err(mcmc_chains, f'jit_{inst}'),
        }
        logger.info("  %s: gamma=%.3f m/s (centered: %.3f), jitter=%.3f m/s",
                     inst, gamma_absolute, gamma_centered, jit)

    # Compute model and residuals on centered data
    model_rv = np.zeros_like(time)
    residuals = np.zeros_like(time)

    for inst in unique_instruments:
        mask = inst_arr == inst
        inst_time = time[mask]
        # Model = Keplerian + gamma_centered (in centered frame)
        gamma_c = instruments_result[inst]["gamma_centered"]
        model_at_inst = mod(inst_time) + gamma_c
        model_rv[mask] = model_at_inst
        residuals[mask] = rv_centered[mask] - model_at_inst
        if use_gp:
            # Subtract the GP conditional mean (activity model) so
            # residuals reflect what neither planets nor activity
            # explain. Points within an instrument keep input order.
            try:
                gp_mean, _ = like_by_inst[inst].predict(inst_time)
                model_rv[mask] += gp_mean
                residuals[mask] -= gp_mean
            except Exception as e:
                logger.warning("GP prediction failed for %s: %s", inst, e)

    rms_after = float(np.std(residuals))
    logger.info("Keplerian fit: RMS %.2f -> %.2f m/s", rms_before, rms_after)

    gp_result = None
    if use_gp:
        gp_result = {
            name: round(float(fitted_params[name].value), 4)
            for name in _GP_HNAMES
        }
        for name in _GP_HNAMES:
            err = _chain_err(mcmc_chains, name)
            if err is not None:
                gp_result[f"{name}_err"] = round(err, 4)
        logger.info("GP hyperparameters: %s", gp_result)

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
        "method": "radvel_keplerian_gp" if use_gp else "radvel_keplerian",
        "gp": gp_result,
        "status": "ok",
    }


def keplerian_residual_analysis(time, rv, rv_err, instruments, planet_params,
                                 exclude_instruments=None,
                                 fix_eccentricities=False,
                                 min_period=1.0, max_period=None,
                                 run_mcmc=False, use_gp=False,
                                 gp_hyperparams=None):
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
        fix_eccentricities=fix_eccentricities,
        run_mcmc=run_mcmc,
        use_gp=use_gp,
        gp_hyperparams=gp_hyperparams,
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
