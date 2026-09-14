"""
================================================================================
Events that do not propagate: isolated downstream droughts and
non-propagating upstream events
--------------------------------------------------------------------------------
Response analysis for Reviewer #2, comment 3:
  "To better understand upstream-downstream drought propagation, it would
  also be important to assess: a. upstream drought events that did not
  propagate downstream; and b. downstream drought events that occurred
  without any upstream trigger, if such cases existed."

This reuses the validated, vectorised matching algorithm (propagation_null_
core.py -- identical to the one behind the published 525 chains, validated
in the R3.2 null-test analysis to reproduce the manuscript's headline
numbers exactly) WITHOUT modification to the accepted-match logic. It only
adds two diagnostic passes:

(b) Downstream (origin) events without an upstream trigger: origin events
    (at the 12 stations with upstream connectivity) that receive ZERO
    accepted upstream matches. This is the direct complement of the
    published 56.6% connected fraction (928 origin events, 525 connected
    -> 403 unconnected, 43.4%).

(a) Upstream events that do not propagate: among all events that were ever
    a WINDOW-FILTERED CANDIDATE for at least one origin event (i.e.
    genuinely "in play" -- not simply any event at an upstream-role
    station, since many such events are temporally far from any origin
    event and were never real candidates), how many were NEVER selected as
    the accepted match for any origin -- either because their overlap fell
    below the minimum threshold M, or because they lost the deterministic
    tie-break to a better candidate at the same station?

Inputs:  data/SSI_drought_events.csv, data/upstream_connectivity.csv
Outputs: origin_events_connectivity_R2_3.csv, upstream_candidate_events_R2_3.csv,
         fig_non_propagating_events.png/pdf, printed statistics.
================================================================================
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared modules live in revision_analyses/
from propagation_null_core import (
    W_DAYS, MIN_OVERLAP, build_station_arrays, parse_connectivity,
    run_matching_vectorized, summarise,
)

DATA = "data/"
OUT = "./"

# ============================================================================
# STEP 1 — Load data and reproduce the published, accepted matches
# ============================================================================
events = pd.read_csv(DATA + "SSI_drought_events.csv", parse_dates=["start_date", "end_date", "peak_date"])
connectivity = parse_connectivity(DATA + "upstream_connectivity.csv")
station_arr = build_station_arrays(events)
events_idx = events.set_index(["station_id", "event_id"])
origin_stations = [sid for sid, chain in connectivity.items() if chain and sid in station_arr]

accepted_pairs = run_matching_vectorized(station_arr, connectivity)
stats = summarise(accepted_pairs, station_arr, connectivity)
print(f"Validation: {stats['n_chains']} connected origin events out of "
      f"{stats['n_origin_events']} ({stats['connected_fraction']:.1f}%) -- "
      f"must match the published 525/928 (56.6%).")
assert stats["n_chains"] == 525 and stats["n_origin_events"] == 928

connected_origin_keys = {(o_sid, o_eid) for (o_sid, o_eid, u_sid, u_eid, lag) in accepted_pairs}
propagated_upstream_keys = {(u_sid, u_eid) for (o_sid, o_eid, u_sid, u_eid, lag) in accepted_pairs}

# ============================================================================
# STEP 2 — PART (b): downstream origin events WITHOUT an upstream trigger
# ============================================================================
origin_rows = []
for origin_sid in origin_stations:
    o = station_arr[origin_sid]
    for i in range(len(o["start"])):
        o_eid = int(o["eid"][i])
        row = events_idx.loc[(origin_sid, o_eid)]
        origin_rows.append(dict(
            origin_station=origin_sid, origin_event_id=o_eid,
            duration=row["duration"], severity=row["severity"], intensity=row["intensity"],
            connected=(origin_sid, o_eid) in connected_origin_keys,
        ))
origin_df = pd.DataFrame(origin_rows)
origin_df.to_csv(OUT + "origin_events_connectivity_R2_3.csv", index=False)
origin_df.to_csv("origin_events_connectivity_R2_3.csv", index=False)

unconnected = origin_df[~origin_df.connected]
connected = origin_df[origin_df.connected]
print(f"\n=== PART (b): downstream origin events without an upstream trigger ===")
print(f"Unconnected: {len(unconnected)} / {len(origin_df)} ({100*len(unconnected)/len(origin_df):.1f}%)")
print(f"Duration  - connected: mean={connected.duration.mean():.1f} median={connected.duration.median():.1f} "
      f"| unconnected: mean={unconnected.duration.mean():.1f} median={unconnected.duration.median():.1f}")
print(f"Severity  - connected: mean={connected.severity.mean():.2f} median={connected.severity.median():.2f} "
      f"| unconnected: mean={unconnected.severity.mean():.2f} median={unconnected.severity.median():.2f}")

by_station = origin_df.groupby("origin_station").agg(
    n_total=("connected", "size"), n_unconnected=("connected", lambda x: (~x).sum()))
by_station["pct_unconnected"] = 100 * by_station.n_unconnected / by_station.n_total
by_station = by_station.sort_values("pct_unconnected", ascending=False)
by_station.to_csv(OUT + "table_S3_unconnected_by_station.csv")
by_station.to_csv("table_S3_unconnected_by_station.csv")
print("\nUnconnected fraction by origin station:")
print(by_station.to_string())

# ============================================================================
# STEP 3 — PART (a): upstream events that never propagate downstream
# ============================================================================
candidate_keys = set()
for origin_sid in origin_stations:
    upstream_sids = [s for s in connectivity[origin_sid] if s in station_arr]
    o = station_arr[origin_sid]
    for i in range(len(o["start"])):
        s_o = o["start"][i]
        w_start = s_o - W_DAYS
        for up_sid in upstream_sids:
            u = station_arr[up_sid]
            mask = (u["end"] >= w_start) & (u["start"] <= s_o)
            for eid in u["eid"][mask]:
                candidate_keys.add((up_sid, int(eid)))

print(f"\n=== PART (a): upstream events that never propagate downstream ===")
print(f"Events ever a window-candidate for >=1 origin: {len(candidate_keys)}")
non_propagating = candidate_keys - propagated_upstream_keys
print(f"Never selected as accepted match ('non-propagating'): {len(non_propagating)} "
      f"({100*len(non_propagating)/len(candidate_keys):.1f}%)")

cand_rows = []
for (sid, eid) in candidate_keys:
    row = events_idx.loc[(sid, eid)]
    cand_rows.append(dict(station_id=sid, event_id=eid, duration=row["duration"],
                           severity=row["severity"], intensity=row["intensity"],
                           propagated=(sid, eid) in propagated_upstream_keys))
cand_df = pd.DataFrame(cand_rows)
cand_df.to_csv(OUT + "upstream_candidate_events_R2_3.csv", index=False)
cand_df.to_csv("upstream_candidate_events_R2_3.csv", index=False)

prop = cand_df[cand_df.propagated]
nonprop = cand_df[~cand_df.propagated]
print(f"Duration  - propagated: mean={prop.duration.mean():.1f} median={prop.duration.median():.1f} "
      f"| non-propagating: mean={nonprop.duration.mean():.1f} median={nonprop.duration.median():.1f}")
print(f"Severity  - propagated: mean={prop.severity.mean():.2f} median={prop.severity.median():.2f} "
      f"| non-propagating: mean={nonprop.severity.mean():.2f} median={nonprop.severity.median():.2f}")

# ============================================================================
# STEP 4 — Publication figure
# ============================================================================
plt.rcParams.update({"font.size": 11.5, "axes.titlesize": 11.5, "axes.labelsize": 11.3,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE = "#2E75B6", "#E1812C"

fig = plt.figure(figsize=(13, 8.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1.1, 1])

# (a) bar chart: % unconnected by station (wide panel spanning row 1)
ax_a = fig.add_subplot(gs[0, :])
stations = by_station.index.astype(str)
ax_a.bar(stations, by_station.pct_unconnected, color=ORANGE, edgecolor="k", linewidth=0.4)
ax_a.set_ylabel("Origin events without\nupstream connection (%)")
ax_a.set_xlabel("Origin station")
ax_a.set_title("(a) Isolated downstream events, by origin station", loc="left", fontsize=11.5)
for i, (n, pct) in enumerate(zip(by_station.n_total, by_station.pct_unconnected)):
    ax_a.text(i, pct + 1.5, f"n={n}", ha="center", fontsize=8, color="grey")
ax_a.set_ylim(0, 90)

def box_panel(ax, data_a, data_b, label_a, label_b, ylabel, title, logscale=False):
    bp = ax.boxplot([data_a, data_b], tick_labels=[label_a, label_b], patch_artist=True,
                     widths=0.55, showfliers=True, flierprops=dict(marker="o", markersize=3, alpha=0.35))
    for patch, color in zip(bp["boxes"], [BLUE, ORANGE]):
        patch.set_facecolor(color); patch.set_alpha(0.55)
    for median in bp["medians"]:
        median.set_color("black"); median.set_linewidth(1.5)
    if logscale:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10.8)

# (b) duration: connected vs unconnected origin events
ax_b = fig.add_subplot(gs[1, 0])
box_panel(ax_b, connected.duration, unconnected.duration, "Connected", "Unconnected",
          "Duration (days)", "(b) Origin events\nduration")

# (c) severity: connected vs unconnected origin events
ax_c = fig.add_subplot(gs[1, 1])
box_panel(ax_c, connected.severity, unconnected.severity, "Connected", "Unconnected",
          "Severity", "(c) Origin events\nseverity", logscale=True)

# (d) duration: propagating vs non-propagating upstream events
ax_d = fig.add_subplot(gs[1, 2])
box_panel(ax_d, prop.duration, nonprop.duration, "Propagates", "Does not\npropagate",
          "Duration (days)", "(d) Upstream events\nduration")

# (e) severity: propagating vs non-propagating upstream events
ax_e = fig.add_subplot(gs[1, 3])
box_panel(ax_e, prop.severity, nonprop.severity, "Propagates", "Does not\npropagate",
          "Severity", "(e) Upstream events\nseverity", logscale=True)

fig.suptitle("Events that do not propagate: isolated downstream droughts (a-c) and\n"
             "non-propagating upstream events (d-e)", fontsize=12.5, y=1.03)
fig.tight_layout()
fig.savefig(OUT + "fig_non_propagating_events.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT + "fig_non_propagating_events.pdf", bbox_inches="tight")
print("\nSaved fig_non_propagating_events.(png|pdf)")
plt.show()
