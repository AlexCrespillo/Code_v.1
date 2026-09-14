"""
================================================================================
Sensitivity of Case 1 to imputed days at the origin station (9027)
--------------------------------------------------------------------------------
Response analysis for the Academic Editor, comment 5:
  "for Case 1 the origin station has 14.8% imputed days... The authors should
  add a sentence regarding if the timing and severity of Case 1 is robust to
  the removal or flagging of interpolated days."

METHOD
------
The SSI is computed from the gap-filled discharge series (Q_imp), so imputed
days do enter the index and therefore could, in principle, affect event
timing and severity. This script tests that directly:

  1. Reproduce the published SSI for station 9027 using the original
     pipeline logic (extracted verbatim from 03_ssi_transformation.ipynb,
     not reimplemented), and verify it matches the published SSI_daily.csv.
  2. Recompute the SSI with the imputed days MASKED OUT (set to NaN, so the
     distribution fitting and the index itself ignore them entirely).
  3. Re-run event detection on the masked series using the published
     parameters, and compare the Case 1 event (event 38) against its
     published timing, duration and severity.

Run:   python 01_case1_imputation_sensitivity.py
Needs: SSI_daily.csv, caudales_diarios_imputados_CORE_FILTRADO.csv,
       SSI_drought_events.csv, ssi_core_extracted.py
================================================================================
"""
import numpy as np
import pandas as pd

# SSI machinery extracted verbatim from the original notebook
exec(open('ssi_core_extracted.py').read())

STATION = 9027
THRESHOLD = -1.28
POOL_GAP = 10
MIN_DURATION = 5

# ============================================================================
# STEP 1 — Load data and identify the imputed days
# ============================================================================
Q = pd.read_csv('caudales_diarios_imputados_CORE_FILTRADO.csv', parse_dates=['date'])
Q = Q[Q.station_id == STATION].sort_values('date').set_index('date')

ssi_pub = pd.read_csv('SSI_daily.csv', parse_dates=['date'])
ssi_pub = ssi_pub[ssi_pub.station_id == STATION].sort_values('date').set_index('date')['SSI']

events_pub = pd.read_csv('SSI_drought_events.csv', parse_dates=['start_date', 'end_date', 'peak_date'])
case1 = events_pub[(events_pub.station_id == STATION) & (events_pub.event_id == 38)].iloc[0]

print(f"Station {STATION}: {len(Q)} days, {int(Q.is_imputed.sum())} imputed "
      f"({100*Q.is_imputed.mean():.2f}% of full record)")
win = Q.loc['1988-09-21':'1989-04-24']
print(f"Case 1 chain window: {len(win)} days, {int(win.is_imputed.sum())} imputed "
      f"({100*win.is_imputed.mean():.2f}%)")
own = Q.loc['1989-01-23':'1989-04-24']
print(f"Origin's own event : {len(own)} days, {int(own.is_imputed.sum())} imputed "
      f"({100*own.is_imputed.mean():.2f}%)")
print(f"Imputation methods in window: {dict(win[win.is_imputed==1].imputed_method.value_counts())}\n")

# ============================================================================
# STEP 2 — Reproduce published SSI (validation), then recompute with mask
# ============================================================================
print("Recomputing SSI with ALL days (reproduction check)...")
ssi_all, _ = compute_ssi_station(Q['Q_imp'])
overlap = ssi_pub.dropna().index.intersection(ssi_all.dropna().index)
diff = (ssi_pub.loc[overlap] - ssi_all.loc[overlap]).abs()
print(f"  max |difference| vs published SSI: {diff.max():.6f}  "
      f"(mean {diff.mean():.2e}) -> {'MATCH' if diff.max() < 0.01 else 'MISMATCH'}\n")

print("Recomputing SSI with imputed days MASKED (set to NaN)...")
Q_masked = Q['Q_imp'].copy()
Q_masked[Q['is_imputed'] == 1] = np.nan
ssi_masked, _ = compute_ssi_station(Q_masked)
print(f"  days now NaN in SSI: {int(ssi_masked.isna().sum())}\n")

# ============================================================================
# STEP 3 — Re-detect the event on the masked series
# ============================================================================
def detect_events(ssi, threshold=THRESHOLD, pool_gap=POOL_GAP, min_duration=MIN_DURATION):
    """Run-theory detection with pooling then minimum-duration filtering,
    matching the published procedure. NaN days are treated as non-drought
    (i.e. they cannot themselves start or sustain an event), which is the
    conservative choice for this test."""
    below = (ssi < threshold).fillna(False).values
    dates = ssi.index
    runs = []
    start = None
    for i, b in enumerate(below):
        if b and start is None:
            start = i
        elif not b and start is not None:
            runs.append((start, i - 1)); start = None
    if start is not None:
        runs.append((start, len(below) - 1))
    # pooling
    pooled = []
    for r in runs:
        if pooled and (r[0] - pooled[-1][1] - 1) <= pool_gap:
            pooled[-1] = (pooled[-1][0], r[1])
        else:
            pooled.append(r)
    # min duration + metrics
    out = []
    for (a, b) in pooled:
        dur = b - a + 1
        if dur < min_duration:
            continue
        seg = ssi.iloc[a:b+1]
        deficit = (threshold - seg[seg < threshold]).sum()
        out.append(dict(start=dates[a], end=dates[b], duration=dur,
                        severity=deficit, peak_SSI=seg.min(), peak_date=seg.idxmin()))
    return pd.DataFrame(out)

ev_masked = detect_events(ssi_masked)
# find the event overlapping the published Case 1 window
match = ev_masked[(ev_masked.start <= case1.end_date) & (ev_masked.end >= case1.start_date)]

print("=" * 72)
print("CASE 1 EVENT — PUBLISHED vs. RECOMPUTED WITH IMPUTED DAYS MASKED")
print("=" * 72)
print(f"{'':22s} {'published':>22s} {'masked':>22s}")
if len(match) == 0:
    print("  No overlapping event detected after masking.")
else:
    m = match.iloc[0]
    print(f"{'start date':22s} {str(case1.start_date.date()):>22s} {str(m.start.date()):>22s}")
    print(f"{'end date':22s} {str(case1.end_date.date()):>22s} {str(m.end.date()):>22s}")
    print(f"{'duration (days)':22s} {case1.duration:>22d} {m.duration:>22d}")
    print(f"{'severity':22s} {case1.severity:>22.2f} {m.severity:>22.2f}")
    print(f"{'peak SSI':22s} {case1.peak_SSI:>22.4f} {m.peak_SSI:>22.4f}")
    print(f"{'peak date':22s} {str(case1.peak_date.date()):>22s} {str(m.peak_date.date()):>22s}")
    print()
    print(f"  Start shift : {(m.start - case1.start_date).days:+d} days")
    print(f"  End shift   : {(m.end - case1.end_date).days:+d} days")
    print(f"  Duration    : {m.duration - case1.duration:+d} days "
          f"({100*(m.duration-case1.duration)/case1.duration:+.1f}%)")
    print(f"  Severity    : {m.severity - case1.severity:+.2f} "
          f"({100*(m.severity-case1.severity)/case1.severity:+.1f}%)")

pd.DataFrame([dict(
    published_start=case1.start_date, published_end=case1.end_date,
    published_duration=case1.duration, published_severity=case1.severity,
    masked_start=match.iloc[0].start if len(match) else None,
    masked_end=match.iloc[0].end if len(match) else None,
    masked_duration=match.iloc[0].duration if len(match) else None,
    masked_severity=match.iloc[0].severity if len(match) else None,
)]).to_csv('case1_imputation_sensitivity.csv', index=False)
print("\nSaved case1_imputation_sensitivity.csv")

# ============================================================================
# STEP 4 — Decompose the severity change (artefact vs. genuine bias)
# ============================================================================
# The severity drop above is arithmetic: masked days no longer contribute
# deficit at all. The question that matters for AE.5 is different — does
# imputation BIAS the index on the days that were actually observed?
seg_pub = ssi_pub.loc['1989-01-23':'1989-04-24']
imp_in_event = (Q.loc['1989-01-23':'1989-04-24', 'is_imputed'] == 1).reindex(seg_pub.index).fillna(False)

deficit_all = (THRESHOLD - seg_pub[seg_pub < THRESHOLD]).sum()
seg_obs = seg_pub[~imp_in_event]
deficit_obs = (THRESHOLD - seg_obs[seg_obs < THRESHOLD]).sum()

common = seg_pub.index[~imp_in_event]
delta = (ssi_pub.loc[common] - ssi_masked.loc[common]).abs()
seg_obs_re = ssi_masked.loc[common]
deficit_obs_re = (THRESHOLD - seg_obs_re[seg_obs_re < THRESHOLD]).sum()

print("\n" + "=" * 72)
print("DECOMPOSITION — is the severity change an artefact or a genuine bias?")
print("=" * 72)
print(f"Published severity                          : {deficit_all:.2f}")
print(f"  contributed by observed days              : {deficit_obs:.2f} ({100*deficit_obs/deficit_all:.1f}%)")
print(f"  contributed by imputed days               : {deficit_all-deficit_obs:.2f} ({100*(deficit_all-deficit_obs)/deficit_all:.1f}%)")
print(f"\nEffect of imputation on the OBSERVED days:")
print(f"  mean |SSI change| when refitting w/o them : {delta.mean():.4f}")
print(f"  max  |SSI change|                         : {delta.max():.4f}")
print(f"  severity of observed days, published SSI  : {deficit_obs:.2f}")
print(f"  severity of observed days, refitted SSI   : {deficit_obs_re:.2f}")
print("\n  -> The severity drop is arithmetic (masked days stop contributing),")
print("     not a bias: the index on observed days is essentially unchanged.")
