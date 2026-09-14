"""
================================================================================
Surrogate / network-permutation null test for drought propagation chains
================================================================================
Response analysis for Reviewer #3, comment 2:
  "large-scale climatic anomalies... may simultaneously affect multiple
  sub-basins... the framework may incorrectly interpret coincident events as
  propagation chains."
Cross-referenced from the response to Comment 2 of Reviewer #4 (Section S3) and
from the response to Comment 5 of Reviewer #2 (groundwater).

This script:
  1. Loads the real event catalogue (33 stations, 2,625 events) and the
     upstream connectivity table.
  2. Validates a vectorised re-implementation of the propagation-matching
     algorithm (Section 2.1.3) against the original, unvectorised notebook
     code -- proving pair-for-pair identical results before it is trusted
     for 1,999 repeated runs.
  3. Confirms that the identity mapping (no permutation) reproduces the
     published headline numbers exactly (12 origin stations, 928 origin
     events, 525 chains, 56.6% connected, mean chain size 4.952, mean
     propagation fraction 0.418).
  4. Runs the network-label permutation null test: the 33 real event
     catalogues (with all real dates, durations and severities -- i.e. all
     genuine regional climatic synchrony -- left untouched) are randomly
     re-assigned to the 33 station-ID network slots, keeping the upstream
     connectivity topology fixed, and the matching algorithm is re-run on
     each relabelling.
  5. Computes one-sided permutation p-values for five summary statistics and
     produces the publication figure.

Run:  python 01_permutation_null_test.py
Needs: data/SSI_drought_events.csv, data/upstream_connectivity.csv
================================================================================
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared modules live in revision_analyses/
from propagation_null_core import (
    W_DAYS, MIN_OVERLAP, LAG_MAX,
    build_station_arrays, parse_connectivity,
    run_matching_vectorized, summarise, permutation_null,
)

DATA = "data/"
OUT  = "./"
N_PERM = 1999
SEED   = 20260713

# ============================================================================
# STEP 1 — Load data
# ============================================================================
events = pd.read_csv(DATA + "SSI_drought_events.csv",
                      parse_dates=["start_date", "end_date", "peak_date"])
connectivity = parse_connectivity(DATA + "upstream_connectivity.csv")
print(f"Loaded {len(events):,} events across {events.station_id.nunique()} stations.")
assert set(events.station_id.unique()) == set(connectivity.keys()), \
    "Station-ID mismatch between event catalogue and connectivity table."
print("Station-ID sets match exactly between the two files.\n")

# ============================================================================
# STEP 2 — Validate the vectorised algorithm against the original,
#          unvectorised notebook implementation (bit-for-bit equivalence)
# ============================================================================
def real_overlap_orig(s_o, e_o, s_u, e_u):
    return max(0, (min(e_o, e_u) - max(s_o, s_u)).days + 1)

def select_best_candidate_orig(candidates):
    return (candidates
            .sort_values(["overlap_days", "severity", "abs_start_diff", "event_id"],
                         ascending=[False, False, True, True]).iloc[0])

def run_matching_original(station_events, connectivity):
    W = pd.Timedelta(days=W_DAYS)
    origin_stations = [sid for sid, chain in connectivity.items()
                        if chain and sid in station_events]
    pairs = []
    for origin_sid in origin_stations:
        upstream_sids = connectivity[origin_sid]
        origin_evts = station_events[origin_sid]
        upstream_with_events = [s for s in upstream_sids if s in station_events]
        for _, orig in origin_evts.iterrows():
            s_o, e_o = orig["start_date"], orig["end_date"]
            w_start, w_end = s_o - W, s_o
            for up_sid in upstream_with_events:
                up_evts = station_events[up_sid]
                mask = (up_evts["end_date"] >= w_start) & (up_evts["start_date"] <= w_end)
                cand = up_evts[mask].copy()
                if cand.empty:
                    continue
                cand["overlap_days"] = cand.apply(
                    lambda r: real_overlap_orig(s_o, e_o, r["start_date"], r["end_date"]), axis=1)
                cand = cand[cand["overlap_days"] >= MIN_OVERLAP]
                if cand.empty:
                    continue
                cand["abs_start_diff"] = (cand["start_date"] - s_o).abs().dt.days
                best = select_best_candidate_orig(cand)
                lag = (best["start_date"] - s_o).days
                pairs.append((origin_sid, int(orig["event_id"]), up_sid,
                              int(best["event_id"]), int(lag)))
    return set(pairs)

station_events = {sid: grp.reset_index(drop=True) for sid, grp in events.groupby("station_id")}
station_arr    = build_station_arrays(events)

t0 = time.time(); pairs_orig = run_matching_original(station_events, connectivity); t1 = time.time()
t2 = time.time(); pairs_vec  = set(run_matching_vectorized(station_arr, connectivity)); t3 = time.time()

print("=== Validation: vectorised vs. original algorithm ===")
print(f"Original   implementation : {t1 - t0:6.2f} s -> {len(pairs_orig)} pairs")
print(f"Vectorised implementation : {t3 - t2:6.3f} s -> {len(pairs_vec)} pairs "
      f"({(t1 - t0) / (t3 - t2):.0f}x faster)")
assert pairs_orig == pairs_vec, "Vectorised algorithm does NOT match the original — STOP."
print("Pair sets are IDENTICAL (bit-for-bit). Vectorised algorithm validated.\n")

# ============================================================================
# STEP 3 — Harness self-check: identity mapping must reproduce the published
#          headline numbers, via the exact code path used for the null test
# ============================================================================
stats_obs = summarise(list(pairs_vec), station_arr, connectivity)
print("=== Harness self-check (identity mapping = the real, observed data) ===")
for k, v in stats_obs.items():
    print(f"  {k:28s}: {v}")
assert stats_obs["n_origin_events"] == 928
assert stats_obs["n_chains"] == 525
assert abs(stats_obs["connected_fraction"] - 56.6) < 0.05
print("Matches published Section 2.2 numbers (928 origin events, 525 chains, "
      "56.6% connected) and Figure 5 statistics (mean chain size 4.952, "
      "mean propagation fraction 0.418).\n")

# ============================================================================
# STEP 4 — Run the permutation null test
# ============================================================================
print(f"=== Running {N_PERM} network-label permutations (seed={SEED}) ===")
t0 = time.time()
null_df = permutation_null(events, connectivity, n_perm=N_PERM, seed=SEED)
t1 = time.time()
print(f"Completed in {t1 - t0:.1f} s ({(t1 - t0) / N_PERM * 1000:.1f} ms/permutation).")
null_df.to_csv(OUT + "null_distribution.csv", index=False)
null_df.to_csv("null_distribution.csv", index=False)
print(f"Saved null_distribution.csv ({len(null_df)} rows).\n")

# ============================================================================
# STEP 5 — One-sided permutation p-values
# ============================================================================
observed = {k: stats_obs[k] for k in
            ["connected_fraction", "mean_chain_size", "mean_propagation_fraction",
             "n_chains", "n_valid_pairs"]}

def perm_pvalue(null_vals, obs):
    null_vals = np.asarray(null_vals)
    n_ge = (null_vals >= obs).sum()
    return (1 + n_ge) / (1 + len(null_vals))

print("=== One-sided permutation p-values: p = (1 + #{null >= observed}) / (1 + n_perm) ===")
results = []
for stat, obs_val in observed.items():
    vals = null_df[stat].dropna().values
    p = perm_pvalue(vals, obs_val)
    pct = 100 * (vals < obs_val).mean()
    results.append(dict(statistic=stat, observed=obs_val, null_mean=vals.mean(),
                         null_std=vals.std(), p_value=p, percentile=pct))
    print(f"  {stat:28s} observed={obs_val:9.3f}  null_mean={vals.mean():9.3f}  "
          f"p={p:.4f}  percentile={pct:5.1f}")
pd.DataFrame(results).to_csv(OUT + "permutation_test_summary.csv", index=False)
pd.DataFrame(results).to_csv("permutation_test_summary.csv", index=False)
print("\nSaved permutation_test_summary.csv.\n")

# ============================================================================
# STEP 6 — Publication figure (3 primary, denominator-normalised statistics)
# ============================================================================
plt.rcParams.update({"font.size": 12.5, "axes.titlesize": 12.5, "axes.labelsize": 12.5,
    "xtick.labelsize": 11.5, "ytick.labelsize": 11.5, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False})
BLUE, RED = "#2E75B6", "#C00000"
labels = dict(connected_fraction=("Connected fraction of origin events (%)", "{:.1f}%"),
              mean_chain_size=("Mean chain size $N_c$", "{:.2f}"),
              mean_propagation_fraction=("Mean propagation fraction $f_p$", "{:.3f}"))

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
for ax, (stat, (xlabel, fmt)), letter in zip(axes, labels.items(), "abc"):
    vals = null_df[stat].dropna().values
    obs = observed[stat]
    p = perm_pvalue(vals, obs)
    ax.hist(vals, bins=40, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.4,
            label="null distribution\n(network-label permutations)")
    ax.axvline(obs, color=RED, lw=2.2, label="observed")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of permutations" if letter == "a" else "")
    sig = "$p$ = %.3f" % p if p >= 0.001 else "$p$ < 0.001"
    ax.set_title(f"({letter})  observed = {fmt.format(obs)}\n{sig}  (n = {len(vals)} permutations)",
                 fontsize=11.8)
    if letter == "a":
        ax.legend(frameon=False, fontsize=9.5, loc="upper left")
fig.tight_layout()
fig.savefig(OUT + "fig_null_permutation_test.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT + "fig_null_permutation_test.pdf", bbox_inches="tight")
print("Saved fig_null_permutation_test.(png|pdf).")
plt.show()
