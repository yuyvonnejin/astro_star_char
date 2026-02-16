// Module 1: Distance computation from parallax or Cepheid period-luminosity relation.
// Ported from src/distance.py

const PARALLAX_ZERO_POINT_MAS = 0.017;
const PRIOR_SCALE_LENGTH_PC = 1350.0;
const GRID_MIN_PC = 1.0;
const GRID_MAX_PC = 100000.0;
const GRID_POINTS = 10000;

function computeDistance(star) {
  if (star.is_cepheid && star.cepheid_period_days != null) {
    return cepheidDistance(star, star.cepheid_period_days);
  }
  return bayesianDistance(star);
}

function cepheidDistance(star, period) {
  const m_v = -2.43 * Math.log10(period) - 4.05;
  const ag = star.ag_gspphot || 0.0;
  const g_0 = star.phot_g_mean_mag - ag;
  const distance_pc = Math.pow(10.0, (g_0 - m_v + 5.0) / 5.0);

  return {
    distance_pc: distance_pc,
    distance_method: "cepheid_leavitt",
    distance_lower_pc: null,
    distance_upper_pc: null,
    // Intermediate values for display
    _intermediates: {
      period_days: period,
      M_V: m_v,
      G_0: g_0,
    },
  };
}

function bayesianDistance(star) {
  const parallax_obs = star.parallax_mas + PARALLAX_ZERO_POINT_MAS;
  const parallax_error = star.parallax_error_mas;

  // Build log-spaced grid
  const r_grid = [];
  const log_min = Math.log10(GRID_MIN_PC);
  const log_max = Math.log10(GRID_MAX_PC);
  for (let i = 0; i < GRID_POINTS; i++) {
    const t = log_min + (log_max - log_min) * i / (GRID_POINTS - 1);
    r_grid.push(Math.pow(10, t));
  }

  // Prior: r^2 * exp(-r / L)
  const prior = r_grid.map(r => r * r * Math.exp(-r / PRIOR_SCALE_LENGTH_PC));

  // Likelihood: Normal(1000/r, parallax_error) evaluated at parallax_obs
  const likelihood = r_grid.map(r => {
    const predicted = 1000.0 / r;
    const z = (parallax_obs - predicted) / parallax_error;
    return Math.exp(-0.5 * z * z) / (parallax_error * Math.sqrt(2 * Math.PI));
  });

  // Posterior (unnormalized)
  const posterior_raw = r_grid.map((_, i) => prior[i] * likelihood[i]);

  // Normalize via trapezoidal integration
  let total = 0;
  for (let i = 1; i < GRID_POINTS; i++) {
    total += 0.5 * (posterior_raw[i] + posterior_raw[i - 1]) * (r_grid[i] - r_grid[i - 1]);
  }

  if (total === 0) {
    const naive = parallax_obs > 0 ? 1000.0 / parallax_obs : null;
    return {
      distance_pc: naive,
      distance_method: "parallax_bayesian",
      distance_lower_pc: null,
      distance_upper_pc: null,
      _intermediates: {parallax_corrected: parallax_obs, naive_distance: naive, posterior_zero: true},
    };
  }

  const posterior = posterior_raw.map(v => v / total);

  // MAP estimate
  let mode_idx = 0;
  let mode_val = posterior[0];
  for (let i = 1; i < GRID_POINTS; i++) {
    if (posterior[i] > mode_val) {
      mode_val = posterior[i];
      mode_idx = i;
    }
  }
  const distance_pc = r_grid[mode_idx];

  // CDF for credible interval
  const cdf = [];
  let cumsum = 0;
  for (let i = 0; i < GRID_POINTS; i++) {
    const dr = i === 0 ? (r_grid[1] - r_grid[0]) : (r_grid[i] - r_grid[i - 1]);
    cumsum += posterior[i] * dr;
    cdf.push(cumsum);
  }
  const cdf_max = cdf[cdf.length - 1];
  for (let i = 0; i < cdf.length; i++) cdf[i] /= cdf_max;

  const lower = linearInterp(0.16, cdf, r_grid);
  const upper = linearInterp(0.84, cdf, r_grid);

  const naive_distance = 1000.0 / parallax_obs;

  return {
    distance_pc: distance_pc,
    distance_lower_pc: lower,
    distance_upper_pc: upper,
    distance_method: "parallax_bayesian",
    _intermediates: {
      parallax_corrected: parallax_obs,
      naive_distance: naive_distance,
      r_grid: r_grid,
      prior: prior,
      likelihood: likelihood,
      posterior: posterior,
    },
  };
}

function linearInterp(target, xp, fp) {
  // Interpolate: given sorted xp array and corresponding fp, find fp value at target x
  for (let i = 1; i < xp.length; i++) {
    if (xp[i] >= target) {
      const t = (target - xp[i - 1]) / (xp[i] - xp[i - 1]);
      return fp[i - 1] + t * (fp[i] - fp[i - 1]);
    }
  }
  return fp[fp.length - 1];
}
