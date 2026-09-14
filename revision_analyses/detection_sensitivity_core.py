"""
================================================================================
Sensitivity of drought-event detection to pooling gap and minimum duration
--------------------------------------------------------------------------------
Response analysis for Reviewer #3, comment 1 (event-detection parameter
sensitivity), cross-referenced from Reviewer #4, comment 4 (parameter
sensitivity/adaptability) and Reviewer #4, comment 15 (search-window
justification, already covered by the existing RSI/Kneedle analysis for W).

CONTEXT
-------
The manuscript's existing sensitivity analysis (Fig. 4, Fig. S1-S2, Table S1)
covers the two PROPAGATION parameters: search window W and minimum overlap M.
It does NOT cover the two EVENT-DETECTION parameters, which operate at an
earlier pipeline stage (on the daily SSI series, before any propagation
matching takes place):
  - POOL_GAP     (default 10 days): merge two candidate drought events if the
                  gap between them is <= this many days.
  - MIN_DURATION (default 5 days) : discard events whose total span is
                  shorter than this.
This module fills that gap.

VALIDATION
----------
Before running any sensitivity sweep, the detection algorithm (reproduced
verbatim in vectorised form from `04_ssi_drought_events.ipynb`) is checked to
reproduce the published Table 3 statistics EXACTLY at the default parameter
values (POOL_GAP=10, MIN_DURATION=5):
  - 2,625 events across 33 stations (79.5 events/station on average)
  - Duration:      mean 33.2, median 19.0, max 613 days
  - Severity:      mean 12.952, median 3.937, max 508.836
  - Vol. severity: mean 5.508, median 0.445, min 0.000, max 362.262 hm3
(all matched to 3-4 decimals in `01_sensitivity_pooling_duration.py`.)

CONTENTS
--------
  - detect_events_for_station()  : the 4-step detection algorithm (flag ->
                                    pool -> filter -> characterise), exactly
                                    as in the original notebook, parametrised
                                    by pool_gap and min_duration.
  - build_q_p10_lookup()         : per-station x DOY 10th-percentile discharge
                                    threshold (computed ONCE; independent of
                                    the detection parameters being swept).
  - compute_severity_hm3()       : volumetric deficit for a set of events,
                                    given the P10 lookup table.
  - run_detection()               : full pipeline for one (pool_gap,
                                    min_duration) combination.
  - sensitivity_sweep()           : runs run_detection() across the full grid
                                    and returns combination-level and
                                    per-station results.
================================================================================
"""
import numpy as np
import pandas as pd

THRESHOLD = -1.28   # SSI drought-day threshold (P10 of standard normal)
M3S_TO_HM3_DAY = 86_400 / 1e6
WINDOW = 15          # +/- days for the P10 discharge climatology (unchanged)

DEFAULT_POOL_GAP = 10
DEFAULT_MIN_DURATION = 5


# ------------------------------------------------------------ event detection --
def detect_raw_events(ssi_values: np.ndarray, threshold: float = THRESHOLD):
    """Runs of consecutive days with SSI < threshold. Returns list of (start,end)
    integer position pairs (inclusive), NaN treated as non-drought (False)."""
    below = ssi_values < threshold
    events = []
    in_event = False
    start = 0
    for i, flag in enumerate(below):
        if flag and not in_event:
            in_event = True
            start = i
        elif not flag and in_event:
            in_event = False
            events.append((start, i - 1))
    if in_event:
        events.append((start, len(below) - 1))
    return events


def pool_events(events, dates_ord: np.ndarray, gap_days: int):
    """Merge consecutive events separated by <= gap_days non-drought days."""
    if not events:
        return []
    pooled = [list(events[0])]
    for start, end in events[1:]:
        prev_end = pooled[-1][1]
        gap = int(dates_ord[start] - dates_ord[prev_end]) - 1
        if gap <= gap_days:
            pooled[-1][1] = end
        else:
            pooled.append([start, end])
    return [tuple(e) for e in pooled]


def filter_events(events, dates_ord: np.ndarray, min_days: int):
    """Discard events whose total span (start to end inclusive) is < min_days."""
    return [(s, e) for s, e in events if int(dates_ord[e] - dates_ord[s]) + 1 >= min_days]


def compute_event_stats(events, ssi_values, dates_ord, dates, station_id, threshold=THRESHOLD):
    records = []
    for event_id, (start, end) in enumerate(events, start=1):
        event_ssi = ssi_values[start:end + 1]
        duration = int(dates_ord[end] - dates_ord[start]) + 1
        drought_mask = event_ssi < threshold
        severity = float(np.sum(np.abs(threshold - event_ssi[drought_mask])))
        intensity = severity / duration
        peak_ssi = float(event_ssi.min())
        records.append(dict(
            station_id=station_id, event_id=event_id,
            start_date=dates[start], end_date=dates[end],
            duration=duration, severity=round(severity, 4),
            intensity=round(intensity, 4), peak_SSI=round(peak_ssi, 4),
        ))
    return records


def detect_events_for_station(ssi_values, dates, pool_gap, min_duration, station_id):
    dates_ord = dates.values.astype("datetime64[D]").astype(np.int64)
    raw = detect_raw_events(ssi_values)
    pooled = pool_events(raw, dates_ord, pool_gap)
    filtered = filter_events(pooled, dates_ord, min_duration)
    return compute_event_stats(filtered, ssi_values, dates_ord, dates, station_id)


# ---------------------------------------------------- volumetric severity (hm3) --
def add_doy(df: pd.DataFrame) -> pd.DataFrame:
    """Returns a copy of df with a 'doy' (day-of-year) column added.
    A plain helper so callers never rely on in-place mutation across
    function boundaries (a classic pandas footgun)."""
    df = df.copy()
    df["doy"] = df["date"].dt.dayofyear
    return df


def build_q_p10_lookup(Q_raw: pd.DataFrame) -> pd.DataFrame:
    """Per-station x DOY 10th-percentile discharge, using a +/-WINDOW-day
    circular pooling window across all years. Independent of the
    detection parameters, so computed ONCE and reused for every combination.
    Expects Q_raw to already contain a 'doy' column (see add_doy())."""
    Q_raw = Q_raw.sort_values(["station_id", "date"]).reset_index(drop=True)
    stations = Q_raw["station_id"].unique()
    doy_range = np.arange(1, 367)
    records = []
    for station in stations:
        sub = Q_raw[Q_raw["station_id"] == station]
        Q_vals, doys = sub["Q_imp"].values, sub["doy"].values
        for doy in doy_range:
            lo, hi = doy - WINDOW, doy + WINDOW
            if lo < 1:
                mask = (doys >= lo + 366) | (doys <= hi)
            elif hi > 366:
                mask = (doys >= lo) | (doys <= hi - 366)
            else:
                mask = (doys >= lo) & (doys <= hi)
            pool = Q_vals[mask]
            pool = pool[~np.isnan(pool)]
            p10 = float(np.percentile(pool, 10)) if len(pool) > 0 else np.nan
            records.append({"station_id": station, "doy": doy, "Q_p10": p10})
    return pd.DataFrame(records)


def compute_severity_hm3(events_df: pd.DataFrame, Q_lookup: pd.DataFrame) -> list:
    """Volumetric deficit (hm3) for each event, given start/end dates and the
    precomputed P10 lookup merged onto the daily discharge series."""
    out = []
    for _, row in events_df.iterrows():
        sub = Q_lookup.loc[row["station_id"]]
        mask = (sub["date"] >= row["start_date"]) & (sub["date"] <= row["end_date"])
        window = sub[mask]
        deficit = np.maximum(0.0, window["Q_p10"] - window["Q_imp"])
        out.append(round(float(deficit.sum() * M3S_TO_HM3_DAY), 4))
    return out


# --------------------------------------------------------------- full pipeline --
def run_detection(df_ssi: pd.DataFrame, Q_lookup: pd.DataFrame,
                   pool_gap: int, min_duration: int) -> pd.DataFrame:
    """Runs the full detection + volumetric-severity pipeline for one
    (pool_gap, min_duration) combination across all stations."""
    all_records = []
    for station, sub in df_ssi.groupby("station_id"):
        sub = sub.sort_values("date").reset_index(drop=True)
        all_records.extend(detect_events_for_station(
            sub["SSI"].values, sub["date"], pool_gap, min_duration, station))
    events_df = pd.DataFrame(all_records)
    if events_df.empty:
        return events_df
    events_df["severity_hm3"] = compute_severity_hm3(events_df, Q_lookup)
    return events_df


# -------------------------------------------------------------- sensitivity sweep --
def sensitivity_sweep(df_ssi: pd.DataFrame, Q_lookup: pd.DataFrame,
                       pool_gaps=(0, 5, 7, 10, 15),
                       min_durations=(3, 5, 7, 10)) -> dict:
    """Runs run_detection() across the full (pool_gap x min_duration) grid.

    Returns a dict with:
      - 'summary' : one row per combination (n_events, mean/median of each
                    metric, and pct change in n_events vs the default)
      - 'per_station' : one row per (combination, station) with event counts
                    and mean severity, used for the bias/homogeneity check
      - 'bias'    : one row per combination with the Spearman correlation of
                    per-station event counts (and per-station mean severity)
                    against the default configuration -- tests whether the
                    parameter choice affects some stations disproportionately
                    more than others (Reviewer #3's explicit concern).
    """
    from scipy import stats as sstats

    summary_rows, per_station_rows = [], []
    default_counts, default_sev = None, None

    combos = [(pg, md) for pg in pool_gaps for md in min_durations]
    for pg, md in combos:
        ev = run_detection(df_ssi, Q_lookup, pg, md)
        n_events = len(ev)
        per_station = (ev.groupby("station_id")
                       .agg(n_events=("event_id", "size"),
                            mean_severity=("severity", "mean"))
                       .reindex(sorted(df_ssi.station_id.unique()), fill_value=0))
        per_station["pool_gap"], per_station["min_duration"] = pg, md
        per_station_rows.append(per_station.reset_index())

        if pg == DEFAULT_POOL_GAP and md == DEFAULT_MIN_DURATION:
            default_counts = per_station["n_events"].values.copy()
            default_sev = per_station["mean_severity"].values.copy()

        summary_rows.append(dict(
            pool_gap=pg, min_duration=md, n_events=n_events,
            mean_duration=ev["duration"].mean(), median_duration=ev["duration"].median(),
            mean_severity=ev["severity"].mean(), median_severity=ev["severity"].median(),
            mean_intensity=ev["intensity"].mean(), median_intensity=ev["intensity"].median(),
            mean_severity_hm3=ev["severity_hm3"].mean(), median_severity_hm3=ev["severity_hm3"].median(),
        ))

    summary = pd.DataFrame(summary_rows)
    default_n = summary.loc[(summary.pool_gap == DEFAULT_POOL_GAP) &
                             (summary.min_duration == DEFAULT_MIN_DURATION), "n_events"].iloc[0]
    summary["pct_change_n_events"] = 100 * (summary["n_events"] - default_n) / default_n

    per_station_df = pd.concat(per_station_rows, ignore_index=True)

    bias_rows = []
    for pg, md in combos:
        sub = per_station_df[(per_station_df.pool_gap == pg) & (per_station_df.min_duration == md)]
        sub = sub.sort_values("station_id")
        rho_n, p_n = sstats.spearmanr(default_counts, sub["n_events"].values)
        rho_s, p_s = sstats.spearmanr(default_sev, sub["mean_severity"].values)
        bias_rows.append(dict(pool_gap=pg, min_duration=md,
                               spearman_rho_n_events=rho_n, spearman_p_n_events=p_n,
                               spearman_rho_mean_severity=rho_s, spearman_p_mean_severity=p_s))
    bias = pd.DataFrame(bias_rows)

    return dict(summary=summary, per_station=per_station_df, bias=bias)
