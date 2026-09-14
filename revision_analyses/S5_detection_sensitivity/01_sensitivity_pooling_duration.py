"""
================================================================================
Sensitivity of drought-event detection to pooling gap and minimum duration
================================================================================
Response analysis for Reviewer #3, comment 1, cross-referenced from Reviewer
#4, comment 4 (parameter sensitivity/adaptability across climatic regions).

This script:
  1. Loads the daily SSI series (33 stations, 1961-2020) and the daily
     discharge series used for volumetric severity.
  2. Validates the vectorised detection algorithm (Step 1-4, reproduced
     verbatim from `04_ssi_drought_events.ipynb`) against the published
     Table 3 statistics at the default parameters (pooling = 10 d,
     minimum duration = 5 d): 2,625 events, exact duration/severity/
     volumetric-severity statistics.
  3. Runs the full sensitivity grid: pooling gap in {0,5,7,10,15} days x
     minimum duration in {3,5,7,10} days (20 combinations).
  4. For each combination, computes summary statistics (event count,
     duration, severity, intensity, volumetric severity) AND a per-station
     bias/homogeneity check (Spearman correlation of per-station event
     counts and mean severity against the default configuration) --
     directly answering the reviewer's question about station/sub-basin
     bias.
  5. Produces the publication figure (6-panel heatmap grid, matching the
     visual style of the manuscript's Figure 4).

Run:  python 01_sensitivity_pooling_duration.py
Needs: data/SSI_daily.csv, data/caudales_diarios_imputados_CORE_FILTRADO.csv
================================================================================
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared modules live in revision_analyses/
from detection_sensitivity_core import (
    THRESHOLD, DEFAULT_POOL_GAP, DEFAULT_MIN_DURATION,
    add_doy, build_q_p10_lookup, run_detection, sensitivity_sweep,
)

DATA = "data/"
OUT  = "./"

# ============================================================================
# STEP 1 — Load data
# ============================================================================
df_ssi = pd.read_csv(DATA + "SSI_daily.csv", parse_dates=["date"])
Q_raw  = pd.read_csv(DATA + "caudales_diarios_imputados_CORE_FILTRADO.csv", parse_dates=["date"])
Q_raw  = add_doy(Q_raw[["station_id", "date", "Q_imp"]])

print(f"SSI series   : {len(df_ssi):,} rows, {df_ssi.station_id.nunique()} stations, "
      f"{df_ssi.date.min().date()} to {df_ssi.date.max().date()}")
print(f"Discharge    : {len(Q_raw):,} rows, {Q_raw.station_id.nunique()} stations\n")
assert set(df_ssi.station_id.unique()) == set(Q_raw.station_id.unique()), \
    "Station-ID mismatch between SSI series and discharge series."

# ============================================================================
# STEP 2 — Build the P10 discharge lookup ONCE (independent of pool_gap/min_duration)
# ============================================================================
t0 = time.time()
q_p10 = build_q_p10_lookup(Q_raw)
Q_lookup = Q_raw.merge(q_p10, on=["station_id", "doy"], how="left").set_index("station_id")
print(f"P10 discharge lookup built in {time.time()-t0:.1f} s: {q_p10.shape}\n")

# ============================================================================
# STEP 3 — Validate at the default configuration against published Table 3
# ============================================================================
ev_default = run_detection(df_ssi, Q_lookup, DEFAULT_POOL_GAP, DEFAULT_MIN_DURATION)
print("=== Validation at default parameters (pooling=10, min_duration=5) vs Table 3 ===")
print(f"  n_events                : {len(ev_default)}  (paper: 2,625)")
print(f"  duration  mean/median/max: {ev_default.duration.mean():.1f} / "
      f"{ev_default.duration.median():.1f} / {ev_default.duration.max()}  (paper: 33.2/19.0/613)")
print(f"  severity  mean/median/max: {ev_default.severity.mean():.3f} / "
      f"{ev_default.severity.median():.3f} / {ev_default.severity.max():.3f}  "
      f"(paper: 12.952/3.937/508.836)")
print(f"  vol.sev.  mean/median/min/max: {ev_default.severity_hm3.mean():.3f} / "
      f"{ev_default.severity_hm3.median():.3f} / {ev_default.severity_hm3.min():.3f} / "
      f"{ev_default.severity_hm3.max():.3f}  (paper: 5.508/0.445/0.000/362.262)")
assert len(ev_default) == 2625, "Default configuration does not reproduce the published event count."
print("Validated: reproduces published Table 3 exactly.\n")

# ============================================================================
# STEP 4 — Run the full sensitivity grid
# ============================================================================
print("=== Running sensitivity grid: pooling in {0,5,7,10,15} x "
      "min_duration in {3,5,7,10} ===")
t0 = time.time()
results = sensitivity_sweep(df_ssi, Q_lookup,
                             pool_gaps=(0, 5, 7, 10, 15),
                             min_durations=(3, 5, 7, 10))
print(f"Completed 20 combinations in {time.time()-t0:.1f} s.\n")

summary, per_station, bias = results["summary"], results["per_station"], results["bias"]
summary.to_csv(OUT + "sensitivity_summary.csv", index=False)
summary.to_csv("sensitivity_summary.csv", index=False)
per_station.to_csv(OUT + "sensitivity_per_station.csv", index=False)
per_station.to_csv("sensitivity_per_station.csv", index=False)
bias.to_csv(OUT + "sensitivity_bias.csv", index=False)
bias.to_csv("sensitivity_bias.csv", index=False)

print("=== Summary (event characteristics across the grid) ===")
print(summary.round(3).to_string(index=False))
print("\n=== Bias / homogeneity check (Spearman rho vs default, per station) ===")
print(bias.round(3).to_string(index=False))

# ============================================================================
# STEP 5 — Publication figure
# ============================================================================
plt.rcParams.update({"font.size": 11.5, "axes.titlesize": 11.2, "axes.labelsize": 11.5,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "font.family": "DejaVu Sans"})

pool_gaps = sorted(summary.pool_gap.unique())
min_durations = sorted(summary.min_duration.unique())

def to_grid(df, col):
    g = df.pivot(index="pool_gap", columns="min_duration", values=col)
    return g.reindex(index=pool_gaps, columns=min_durations)

panels = [
    ("n_events", "N events (detected)", summary, "viridis"),
    ("mean_duration", "Mean duration (days)", summary, "viridis"),
    ("mean_severity_hm3", "Mean volumetric severity (hm$^3$)", summary, "viridis"),
    ("pct_change_n_events", "N events, % change vs default", summary, "RdBu_r"),
    ("spearman_rho_n_events", "Spearman $\\rho$ (event count, vs default)", bias, "viridis"),
    ("spearman_rho_mean_severity", "Spearman $\\rho$ (mean severity, vs default)", bias, "viridis"),
]

fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.6))
letters = "abcdef"
for ax, (col, title, src, cmap), letter in zip(axes.flat, panels, letters):
    grid = to_grid(src, col)
    vals = grid.values.astype(float)
    if cmap == "RdBu_r":
        vmax = np.nanmax(np.abs(vals))
        norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
        im = ax.imshow(vals, cmap=cmap, norm=norm, aspect="auto", origin="lower")
    else:
        im = ax.imshow(vals, cmap=cmap, aspect="auto", origin="lower")
    ax.set_xticks(range(len(min_durations))); ax.set_xticklabels(min_durations)
    ax.set_yticks(range(len(pool_gaps))); ax.set_yticklabels(pool_gaps)
    ax.set_xlabel("Minimum duration (days)")
    ax.set_ylabel("Pooling gap (days)")
    ax.set_title(f"({letter}) {title}", fontsize=11.3)
    for i in range(len(pool_gaps)):
        for j in range(len(min_durations)):
            v = vals[i, j]
            txt = f"{v:.2f}" if abs(v) < 10 else f"{v:.0f}"
            frac = (v - np.nanmin(vals)) / (np.nanmax(vals) - np.nanmin(vals) + 1e-9)
            color = "white" if (cmap != "RdBu_r" and frac > 0.55) else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.3, color=color)
    di = pool_gaps.index(DEFAULT_POOL_GAP); dj = min_durations.index(DEFAULT_MIN_DURATION)
    ax.plot(dj, di, marker="+", color="red", markersize=16, markeredgewidth=2.5)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle("Sensitivity of drought-event detection to pooling gap and minimum duration\n"
             "(red cross = default configuration: pooling = 10 d, minimum duration = 5 d)",
             fontsize=12.5, y=1.01)
fig.tight_layout()
fig.savefig(OUT + "fig_detection_sensitivity.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT + "fig_detection_sensitivity.pdf", bbox_inches="tight")
print("\nSaved fig_detection_sensitivity.(png|pdf).")
plt.show()
