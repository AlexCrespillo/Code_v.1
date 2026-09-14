"""
================================================================================
Diagnostic check: how many candidate connections are excluded by requiring
temporal overlap (M >= 5 days), as opposed to allowing a temporal GAP between
the end of an upstream event and the start of the origin event?
--------------------------------------------------------------------------------
Response analysis for Reviewer #2, comment 2 ("Why is it not plausible... for
an upstream-downstream drought propagation event to occur without temporal
overlap?").

This reuses the validated, vectorised algorithm core from the R3.2 null-test
module (propagation_null_core.py) without modification to its accepted-match
logic; it only ADDS a diagnostic pass that looks, for every (origin_event,
upstream_station) pair that currently has NO accepted match (overlap >= M),
for the best (smallest) temporal GAP candidate within the existing search
window, and records its properties. It does not alter the published matching
algorithm or results in any way.

Definitions:
  - overlap = 0 candidate: the upstream event's interval does not intersect
    the origin event's interval at all.
  - gap = s_o - end_u (days): the number of days between the end of the
    upstream event and the start of the origin event, for candidates with
    overlap = 0 and end_u < s_o (upstream event already finished).
  - A gap of 0 means the upstream event ended the day before the origin
    event started (adjacent, no overlap). Only candidates within the
    existing search window (end_u >= s_o - W) are considered, exactly as
    in the published algorithm -- no new window is introduced.
================================================================================
"""
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # shared modules live in revision_analyses/
from propagation_null_core import (
    W_DAYS, MIN_OVERLAP, build_station_arrays, parse_connectivity,
    run_matching_vectorized,
)

DATA = "data/"
OUT = "./"
GAP_THRESHOLDS = [0, 5, 10, 15, 20, 30, 45]  # days; 45 = full window W

events = pd.read_csv(DATA + "SSI_drought_events.csv", parse_dates=["start_date", "end_date", "peak_date"])
connectivity = parse_connectivity(DATA + "upstream_connectivity.csv")
station_arr = build_station_arrays(events)

# ------------------------------------------------------------------------
# 1. Reproduce the published, ACCEPTED matches (overlap >= M) -- unchanged.
# ------------------------------------------------------------------------
accepted_pairs = run_matching_vectorized(station_arr, connectivity)
accepted_keys = {(o_sid, o_eid, u_sid) for (o_sid, o_eid, u_sid, u_eid, lag) in accepted_pairs}
print(f"Published (accepted) origin-event x upstream-station connections: {len(accepted_keys)}")

# ------------------------------------------------------------------------
# 2. For every (origin_event, upstream_station) WITHOUT an accepted match,
#    scan the SAME search window for the best (smallest) gap candidate
#    with overlap = 0 (i.e. the upstream event already ended).
# ------------------------------------------------------------------------
origin_stations = [sid for sid, chain in connectivity.items() if chain and sid in station_arr]

gap_candidates = []  # one row per (origin_event, upstream_station) that has no accepted match
for origin_sid in origin_stations:
    upstream_sids = [s for s in connectivity[origin_sid] if s in station_arr]
    o = station_arr[origin_sid]
    n_o = len(o["start"])
    for i in range(n_o):
        s_o, e_o = o["start"][i], o["end"][i]
        w_start = s_o - W_DAYS
        o_eid = int(o["eid"][i])
        for up_sid in upstream_sids:
            if (origin_sid, o_eid, up_sid) in accepted_keys:
                continue  # already has a valid overlap>=M match; not relevant here
            u = station_arr[up_sid]
            # same window filter as the published algorithm
            mask = (u["end"] >= w_start) & (u["start"] <= s_o)
            if not mask.any():
                continue
            s2, e2, sv2, id2 = u["start"][mask], u["end"][mask], u["sev"][mask], u["eid"][mask]
            overlap = np.maximum(0, np.minimum(e_o, e2) - np.maximum(s_o, s2) + 1)
            zero_ov = overlap == 0
            if not zero_ov.any():
                continue
            # true gap = days between end of upstream event and start of origin event
            gaps = s_o - e2[zero_ov]
            gaps = gaps[gaps >= 0]  # keep only genuinely-ended-before candidates
            if len(gaps) == 0:
                continue
            best_idx_local = np.argmin(gaps)
            best_gap = int(gaps[best_idx_local])
            # retrieve the corresponding event's own properties for characterisation
            ev_idx = np.where(zero_ov)[0][np.where(e2[zero_ov] == (s_o - best_gap))[0][0]]
            gap_candidates.append(dict(
                origin_station=origin_sid, origin_event_id=o_eid, upstream_station=up_sid,
                best_gap_days=best_gap, upstream_severity=float(sv2[ev_idx]),
                upstream_duration=int(e2[ev_idx] - s2[ev_idx] + 1),
            ))

gap_df = pd.DataFrame(gap_candidates)
gap_df.to_csv(OUT + "gap_candidates_R2_2.csv", index=False)
print(f"Currently-unmatched (origin_event, upstream_station) pairs with a "
      f"zero-overlap gap candidate in the search window: {len(gap_df)}")

# ------------------------------------------------------------------------
# 3. For each gap tolerance G, count how many NEW pairs would be gained,
#    and characterise them against the currently accepted matches.
# ------------------------------------------------------------------------
accepted_sev = np.array([sv for (o, e, u, ue, l) in [] ])  # placeholder, filled below properly
# Recompute accepted severities properly from the matching output
acc_df = pd.DataFrame(accepted_pairs, columns=["origin_station", "origin_event_id",
                                                "upstream_station", "upstream_event_id", "lag_days"])
# attach upstream severity/duration for the accepted matches, for comparison
sev_lookup = events.set_index(["station_id", "event_id"])[["severity", "duration"]]
acc_df = acc_df.join(sev_lookup, on=["upstream_station", "upstream_event_id"])

print("\n=== Sensitivity of connectivity to allowing a temporal GAP (no overlap) ===")
print(f"{'G (days)':>9} | {'new pairs':>9} | {'% of accepted (n=2075)':>22} | "
      f"{'median gap sev.':>15} | {'median accepted sev.':>20}")
rows = []
for G in GAP_THRESHOLDS:
    sub = gap_df[gap_df.best_gap_days <= G]
    pct = 100 * len(sub) / len(accepted_keys)
    med_gap_sev = sub.upstream_severity.median() if len(sub) else np.nan
    print(f"{G:9d} | {len(sub):9d} | {pct:21.1f}% | {med_gap_sev:15.2f} | "
          f"{acc_df.severity.median():20.2f}")
    rows.append(dict(G_days=G, n_new_pairs=len(sub), pct_of_accepted=pct,
                      median_severity_new=med_gap_sev, median_duration_new=sub.upstream_duration.median() if len(sub) else np.nan))

summary = pd.DataFrame(rows)
summary.to_csv(OUT + "gap_sensitivity_summary_R2_2.csv", index=False)
print("\nAccepted-match (published) upstream severity: mean=%.2f median=%.2f" %
      (acc_df.severity.mean(), acc_df.severity.median()))
print("Saved gap_candidates_R2_2.csv and gap_sensitivity_summary_R2_2.csv")
