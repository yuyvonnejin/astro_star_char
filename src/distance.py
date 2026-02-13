"""Module 1: Distance computation from parallax or Cepheid period-luminosity relation."""

import logging
import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

logger = logging.getLogger(__name__)

PARALLAX_ZERO_POINT_MAS = 0.017
PRIOR_SCALE_LENGTH_PC = 1350.0
GRID_MIN_PC = 1.0
GRID_MAX_PC = 100_000.0
GRID_POINTS = 10_000


def compute_distance(star_dict):
    """Compute distance for a star using either Cepheid Leavitt law or Bayesian parallax.

    Parameters
    ----------
    star_dict : dict
        Star record with at minimum: parallax_mas, parallax_error_mas,
        phot_g_mean_mag, ag_gspphot, is_cepheid, cepheid_period_days.

    Returns
    -------
    dict
        Distance results including distance_pc, distance_method, and optionally
        distance_lower_pc / distance_upper_pc for Bayesian method.
    """
    is_cepheid = star_dict.get("is_cepheid", False)
    period = star_dict.get("cepheid_period_days")

    if is_cepheid and period is not None:
        return _cepheid_distance(star_dict, period)
    else:
        return _bayesian_distance(star_dict)


def _cepheid_distance(star_dict, period):
    """Method B: Cepheid period-luminosity (Leavitt Law) distance."""
    m_v = -2.43 * np.log10(period) - 4.05
    g_0 = star_dict["phot_g_mean_mag"] - star_dict.get("ag_gspphot", 0.0)
    distance_pc = 10.0 ** ((g_0 - m_v + 5.0) / 5.0)

    logger.info(
        "Cepheid distance: period=%.2f d, M_V=%.3f, G_0=%.3f, distance=%.1f pc",
        period, m_v, g_0, distance_pc,
    )
    return {
        "distance_pc": distance_pc,
        "distance_method": "cepheid_leavitt",
    }


def _bayesian_distance(star_dict):
    """Method A: Bayesian parallax inversion with exponentially decreasing space density prior."""
    parallax_obs = star_dict["parallax_mas"] + PARALLAX_ZERO_POINT_MAS
    parallax_error = star_dict["parallax_error_mas"]

    r_grid = np.logspace(
        np.log10(GRID_MIN_PC), np.log10(GRID_MAX_PC), GRID_POINTS
    )

    # Prior: r^2 * exp(-r / L)
    prior = r_grid ** 2 * np.exp(-r_grid / PRIOR_SCALE_LENGTH_PC)

    # Likelihood: Normal(1000/r, parallax_error) evaluated at parallax_obs
    predicted_parallax = 1000.0 / r_grid
    likelihood = norm.pdf(parallax_obs, loc=predicted_parallax, scale=parallax_error)

    # Posterior (unnormalized)
    posterior = prior * likelihood

    # Normalize using trapezoidal integration
    total = trapezoid(posterior, r_grid)
    if total == 0:
        logger.warning("Posterior integrates to zero; returning naive distance estimate")
        naive_dist = 1000.0 / parallax_obs if parallax_obs > 0 else None
        return {
            "distance_pc": naive_dist,
            "distance_method": "parallax_bayesian",
            "distance_lower_pc": None,
            "distance_upper_pc": None,
        }
    posterior /= total

    # Mode (MAP estimate)
    mode_idx = np.argmax(posterior)
    distance_pc = r_grid[mode_idx]

    # Credible interval: 16th and 84th percentiles via CDF
    cdf = np.cumsum(posterior * np.gradient(r_grid))
    cdf /= cdf[-1]
    lower = np.interp(0.16, cdf, r_grid)
    upper = np.interp(0.84, cdf, r_grid)

    logger.info(
        "Bayesian distance: parallax_corr=%.3f mas, d=%.3f pc [%.3f, %.3f]",
        parallax_obs, distance_pc, lower, upper,
    )
    return {
        "distance_pc": distance_pc,
        "distance_lower_pc": lower,
        "distance_upper_pc": upper,
        "distance_method": "parallax_bayesian",
    }
