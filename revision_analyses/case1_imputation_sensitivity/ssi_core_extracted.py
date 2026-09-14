# ── Cell 1: Imports ──────────────────────────────────────────────────────────
# Standard scientific stack + scipy.stats for distribution fitting
# and Shapiro-Wilks test. tqdm provides per-station progress bars.

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import shapiro
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore')   # suppress scipy convergence warnings

print('Libraries loaded successfully.')
# ── Cell 3: Define the six candidate distributions ────────────────────────────
# These are the same distribution families commonly used in SPI/SSI literature
# and consistent with the FlexDroughtIndex R package for bottom-truncated data.
#
#   Gamma       — classic SPI distribution; shape + scale
#   Log-Normal  — heavy right tail; common for streamflow
#   Weibull     — flexible; good for low-flow regimes
#   Pearson III — 3-parameter; used in SPEI (unbounded) but valid here too
#   Log-Logistic (Fisk) — often outperforms Gamma for skewed hydro data
#   GEV         — extreme-value family; captures heavy tails

DISTRIBUTIONS = {
    'Gamma'      : stats.gamma,
    'LogNormal'  : stats.lognorm,
    'Weibull'    : stats.weibull_min,
    'PearsonIII' : stats.pearson3,
    'LogLogistic': stats.fisk,
    'GEV'        : stats.genextreme,
}

print('Candidate distributions:')
for name, dist in DISTRIBUTIONS.items():
    print(f'  {name:<12} → scipy.stats.{dist.name}')
# ── Cell 4: Helper functions ──────────────────────────────────────────────────

def doy_circular_distance(doy_array, target_doy):
    """
    Circular distance between an array of DOY values and a target DOY.
    DOY wraps at 365 (e.g., distance between DOY 1 and DOY 360 = 6).
    This handles pooling windows that straddle the year boundary.
    """
    diff = np.abs(doy_array.astype(int) - int(target_doy))
    return np.minimum(diff, 365 - diff)


def fit_and_score(dist, values):
    """
    Fit a scipy distribution to `values` and return (params, W_statistic).

    Steps:
      1. MLE fit via dist.fit(values)
      2. CDF → probability integral transform → standard normal
      3. Shapiro–Wilks W on the transformed series

    Returns (-inf) W if the fit fails or produces degenerate results.
    """
    try:
        params = dist.fit(values)
        p = dist.cdf(values, *params)
        # Clip to avoid ±inf after norm.ppf
        p = np.clip(p, 1e-6, 1 - 1e-6)
        z = stats.norm.ppf(p)
        # Shapiro–Wilks requires 3 ≤ n ≤ 5000; pool size is ~1860 → fine
        w, _ = shapiro(z)
        return params, float(w)
    except Exception:
        return None, -np.inf


print('Helper functions defined.')
# ── Cell 5: Core SSI computation function (single station) ────────────────────
#
# For each unique DOY in the series:
#   a) Pool all Q values within circular DOY ± WINDOW days across all years.
#   b) Handle zero-flow: if proportion p0 > 0, fit only on positive values
#      and shift the CDF: F_adj(q) = p0 + (1-p0)*F(q)  for q > 0
#                         F_adj(0) = p0 * uniform ≈ p0/2  (midpoint)
#   c) Fit 6 distributions on the pool; select best by Shapiro-Wilks W.
#   d) Apply the best-fit CDF to the actual DOY observations → norm.ppf → SSI.
#
# Returns:
#   ssi        : pd.Series aligned to input index
#   stats_log  : dict {doy: {'best': str, 'W': {dist_name: float}}}

WINDOW = 15   # ±15 days pooling window

def compute_ssi_station(series, window=WINDOW):
    """
    Compute daily SSI for a single station.

    Parameters
    ----------
    series : pd.Series
        Daily Q values with a DatetimeIndex. May contain NaNs.
    window : int
        Half-width of the DOY pooling window in days (default 15).

    Returns
    -------
    ssi : pd.Series  (same index as input, float)
    stats_log : dict
    """
    doy_arr  = series.index.dayofyear          # integer 1-366
    ssi      = pd.Series(np.nan, index=series.index, name='SSI')
    stats_log = {}

    for d in sorted(series.index.dayofyear.unique()):

        # ── (a) Build pool ────────────────────────────────────────────────
        mask_pool  = doy_circular_distance(doy_arr, d) <= window
        pool_vals  = series[mask_pool].dropna().values

        if len(pool_vals) < 30:          # safeguard: skip if too few obs
            continue

        # ── (b) Zero-flow probability mass ───────────────────────────────
        p0        = float(np.mean(pool_vals == 0))
        fit_vals  = pool_vals[pool_vals > 0] if p0 > 0 else pool_vals

        if len(fit_vals) < 10:
            continue

        # ── (c) Fit all distributions; select best by Shapiro-Wilks W ───
        best_name, best_params, best_w = None, None, -np.inf
        sw_scores = {}

        for name, dist in DISTRIBUTIONS.items():
            params, w = fit_and_score(dist, fit_vals)
            sw_scores[name] = round(w, 5)
            if w > best_w:
                best_w, best_name, best_params = w, name, params

        stats_log[d] = {'best': best_name, 'W': sw_scores}

        if best_params is None:
            continue

        # ── (d) Apply best CDF to actual DOY observations ────────────────
        mask_doy    = (doy_arr == d)
        actual      = series[mask_doy].dropna()

        if len(actual) == 0:
            continue

        best_dist = DISTRIBUTIONS[best_name]
        p = best_dist.cdf(actual.values, *best_params)

        # Zero-flow adjustment on the output probabilities
        if p0 > 0:
            p = np.where(
                actual.values == 0,
                p0 / 2,                          # midpoint of the zero mass
                p0 + (1.0 - p0) * p              # shift CDF for positive Q
            )

        p = np.clip(p, 1e-6, 1 - 1e-6)
        ssi[actual.index] = stats.norm.ppf(p)

    return ssi, stats_log


print(f'SSI function defined. Pooling window: ±{WINDOW} days.')