# Hydrological Drought Propagation Analysis

This document describes the algorithm implemented in `05_ssi_propagation.ipynb` (and, with identical results, in its vectorised version `05b_ssi_propagation_vectorised.ipynb`), which detects and characterises **drought propagation chains** across river networks — sequences of drought events that travel from upstream tributaries down to a downstream outlet station.

---

## Overview

A **drought propagation chain** links one drought event at a downstream (origin/outlet) station to all drought events at upstream stations that are temporally compatible with it. The goal is to determine: *when a drought occurs at a river outlet, how many upstream stations were already experiencing drought, and how far in advance did those upstream droughts begin?*

The analysis produces two outputs:
- A **pair table**: every matched upstream–downstream event pair.
- A **chain summary table**: one row per origin drought event, with aggregate metrics.

> **Design note — no chain merging.** Each origin event produces exactly one propagation chain. No post-hoc merging of consecutive or overlapping chains is applied. See the [rationale](#no-merge-rationale) below.

---

## Inputs

| File | Description |
|---|---|
| `SSI_drought_events.csv` | Drought event catalogue (one row per event per station). Produced by `04_ssi_drought_events.ipynb`. |
| `upstream_connectivity.csv` | For each station, the list of all hydraulically upstream gauging stations. |

The event catalogue contains, for each event: `station_id`, `event_id`, `start_date`, `end_date`, `duration`, `severity` (SSI-based), and `severity_hm3` (volume-based).

The connectivity file is semicolon-separated; the `upstream_chain` column is a comma-separated list of station IDs.

---

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `W_DAYS` | 45 | Days subtracted from `origin_start` to define the left boundary of the search window: `[origin_start − W_DAYS, origin_start]` |
| `MIN_OVERLAP` | 5 | Minimum real overlap (days) required between an upstream event and the origin event |
| `LAG_MAX` | 0 | Maximum lag allowed — only upstream events that **started before or on the same day** as the origin event are kept |

---

## Algorithm: Step by Step

### Step 1 — Load and index data

The event catalogue is loaded and indexed into a dictionary `station_events` keyed by station ID. This allows O(1) lookups of all events at any station during the main loop.

The connectivity file is parsed into a dictionary `connectivity` mapping each station to its list of upstream stations.

**Example:**
```
station 9002 → upstream stations: [9010, 9015, 9020, 9025]
```

---

### Step 2 — Build the asymmetric search window

For each drought event at an origin (outlet) station, the search window is defined as:

```
w_start = origin_start − W_DAYS   (left side: expanded to capture precursors)
w_end   = origin_start             (right side: capped at origin start)
```

**Why this exact shape?** The lag filter (Step 6) requires `upstream_start ≤ origin_start`. This means:
- Any upstream event starting *after* `origin_start` is physically inadmissible — it cannot represent a drought that propagated to the outlet.
- Setting `w_end = origin_start` eliminates this class of candidates at the window stage, before any overlap or tie-breaking computation.
- The window is therefore exactly aligned with the lag constraint: every candidate that enters the pipeline can, in principle, survive all subsequent filters.

The left-side expansion by `W_DAYS` is essential: it captures precursor upstream droughts that began weeks before the origin event onset but overlapped with it.

**Example:**
```
Origin event at station 9002:
  start_date = 1990-06-01
  end_date   = 1990-08-15

Search window:
  w_start = 1990-04-17  (June 1 − 45 days)  ← expanded to catch precursors
  w_end   = 1990-06-01  (= origin start)     ← capped: no post-onset starts admitted
```

---

### Step 3 — Find candidate upstream events

For each upstream station in the network, find all events whose interval intersects the expanded window.

```
Search window: 1990-04-17 → 1990-06-01

Upstream station 9010 has events:
  Event A: 1990-05-10 → 1990-07-20  ✓ started before origin_start, ends after w_start
  Event B: 1990-07-01 → 1990-09-30  ✗ starts after origin_start (w_end)
  Event C: 1989-01-01 → 1990-04-16  ✗ ends before w_start
```

---

### Step 4 — Apply the real overlap filter

After the expanded window pre-filter, only candidates with at least `MIN_OVERLAP` = 5 days of **real** (non-expanded) overlap with the origin event are kept.

The real overlap is computed as:
```
overlap = max(0, min(origin_end, upstream_end) − max(origin_start, upstream_start) + 1 day)
```

The +1 accounts for inclusive date counting (Jan 1 to Jan 3 = 3 days, not 2).

**Example:**
```
Origin event  :  1990-06-01 → 1990-08-15
Upstream event:  1990-05-10 → 1990-07-20

Real overlap = min(Aug 15, Jul 20) − max(Jun 1, May 10) + 1
             = Jul 20 − Jun 1 + 1
             = 50 days  ✓ (≥ 5)
```

```
Upstream event:  1990-08-12 → 1990-09-30

Real overlap = min(Aug 15, Sep 30) − max(Jun 1, Aug 12) + 1
             = Aug 15 − Aug 12 + 1
             = 4 days  ✗ (< 5, rejected)
```

---

### Step 5 — Select the best candidate per upstream station

If multiple upstream events at the same station pass the overlap filter, exactly one is chosen using a four-level deterministic tie-breaking rule:

| Priority | Criterion | Direction |
|---|---|---|
| 1 | `overlap_days` | Highest first |
| 2 | `severity` | Highest first |
| 3 | `abs_start_diff` | Closest start date to origin |
| 4 | `event_id` | Lowest ID (reproducible tie-break) |

**Example:**
```
Station 9010 has two qualifying events:
  Event A: overlap=50 days, severity=−1.8  ← winner (more overlap)
  Event B: overlap=30 days, severity=−2.5
```

---

### Step 6 — Compute the temporal lag

For the winning upstream event, the lag is:

```
lag_days = upstream_start_date − origin_start_date
```

A **negative lag** means the upstream drought started before the downstream one — the expected pattern for drought propagation. A lag of 0 means simultaneous onset.

**Example:**
```
Origin start   : 1990-06-01
Upstream start : 1990-05-10

lag_days = May 10 − Jun 1 = −22 days
```
The upstream drought began 22 days before the downstream event.

---

### Step 7 — Apply the lag filter

Only pairs with `lag_days ≤ 0` are retained. Because the search window already caps `w_end = origin_start`, all candidates entering this step have `upstream_start ≤ origin_start` by construction, so the lag filter is always satisfied. It is kept as an explicit safeguard to catch any edge cases (e.g. same-day starts due to timestamp ties).

```
Before lag filter: N pairs   (all should already have lag ≤ 0 by window design)
After lag filter:  M pairs   (M ≈ N in practice)
```

---

### Step 8 — Compute chain metrics per origin event

For each origin event, two metrics summarise the spatial reach of the drought:

| Metric | Formula | Meaning |
|---|---|---|
| `chain_size` | Number of upstream stations with a valid (lag ≤ 0) match | Raw count of upstream stations in the chain |
| `n_upstream_with_events` | Upstream stations present in the event catalogue | Denominator for normalisation |
| `propagation_fraction` | `chain_size / n_upstream_with_events` | Fraction of available upstream network affected [0–1] |

**Why normalise?** A chain size of 4 means different things for a station with 5 upstream neighbours vs. one with 20. The propagation fraction makes stations comparable.

**Example:**
```
Origin station 9002 has 4 upstream stations in the catalogue.
After lag filter, 3 of them matched the origin event.

chain_size = 3
n_upstream_with_events = 4
propagation_fraction = 3/4 = 0.75
```

---

### Step 9 — Build the chain summary table

Each origin event becomes one row in the `chains` summary table. The chain temporal window is constructed as follows:

| Field | Definition |
|---|---|
| `chain_start` | Earliest upstream event start date (always ≤ `origin_start`, since all lags ≤ 0) |
| `chain_end` | Latest date across all upstream event ends and the origin event end |
| `chain_duration` | `chain_end − chain_start + 1` (days) |
| `chain_size` | Upstream stations matched + 1 (origin station itself) |
| `lag_mean` | Mean lag across all upstream–origin pairs in the chain |
| `chain_severity` | Sum of SSI-based severity across all events (upstream + origin) |
| `chain_severity_hm3` | Sum of volume-based severity (hm³) across all events |

**Example:**
```
Chain 9002_7 (origin station 9002, event 7):

  origin_start : 1990-06-01    origin_end : 1990-08-15
  upstream matches:
    9010: 1990-05-10 → 1990-07-20  lag = −22 days
    9015: 1990-05-25 → 1990-08-30  lag = −7 days
    9020: 1990-06-01 → 1990-09-10  lag = 0 days

chain_start    = 1990-05-10  (earliest upstream start)
chain_end      = 1990-09-10  (latest end across all events)
chain_duration = 124 days
chain_size     = 3 upstream + 1 origin = 4 stations
lag_mean       = (−22 + −7 + 0) / 3 = −9.7 days
```

---

## No-merge rationale {#no-merge-rationale}

An earlier version of this notebook included a Step 10 that merged chains whose `[chain_start, chain_end]` intervals overlapped, with the intention of consolidating the same drought episode captured by multiple consecutive origin events. This step was removed for two reasons:

**1. Merging operated at the wrong level.** Chain intervals overlap because of long-duration *upstream* events, not because of continuity at the *origin* station. Two origin drought events that are months apart (e.g. April and November of the same year) can produce overlapping chain windows simply because a single upstream gauge remained in continuous drought throughout. Merging them creates a synthetic, multi-month episode that has no physical counterpart in the origin station record.

**2. Loss of origin event integrity.** After merging, the `origin_end` field was set to the end of the *last* sub-chain's origin event. This meant the merged "origin event" encompassed long recovery periods above the SSI threshold, directly contradicting the event definition used in `04_ssi_drought_events.ipynb`. For example, a merged chain at station 9027 spanning April 2011 to January 2013 (696 days) actually comprised four separate origin events separated by months of normal flow — none of which individually exceeded a few weeks to months.

The drought event catalogue already guarantees non-overlapping, continuous events at each station. There is therefore no true double-counting: each origin event is a distinct, well-defined, physically consistent observation. Computing all metrics independently per event is the most transparent and reproducible approach for peer-reviewed publication.

---

## Outputs

### `propagation_W45_minov5.csv`
All compatible upstream–downstream pairs before the lag filter. Useful for sensitivity analysis.

| Column | Description |
|---|---|
| `origin_station` | Downstream outlet station |
| `origin_event_id` | Event ID at the origin station |
| `origin_start/end` | Origin event dates |
| `upstream_station` | Upstream station |
| `upstream_event_id` | Best matched upstream event |
| `upstream_start/end` | Upstream event dates |
| `overlap_days` | Real overlap with origin event |
| `lag_days` | `upstream_start − origin_start` |
| `upstream_severity` | SSI-based severity of the upstream event |

### `propagation_W45_minov5_lagneg.csv`
Lag-filtered pairs (lag ≤ 0), with chain metrics appended.
Same columns as above, plus `chain_size`, `n_upstream_with_events`, `propagation_fraction`.

### `chains_summary_W45_minov5.csv` / `Chains_final.csv`
**Primary analytical output.** One row per origin event. `Chains_final.csv` is an alias used by downstream analysis and visualisation scripts.

| Column | Description |
|---|---|
| `chain_id` | `{origin_station}_{origin_event_id}` |
| `origin_station` | Outlet station ID |
| `origin_start/end` | Origin event dates (single continuous drought event) |
| `chain_start` | Earliest upstream event start |
| `chain_end` | Latest event end in the chain |
| `chain_duration` | Chain span in days |
| `chain_size` | Total stations involved (upstream + origin) |
| `propagation_fraction` | Fraction of upstream network affected |
| `lag_mean` | Mean upstream lag (days) |
| `chain_severity` | Total SSI-based severity |
| `chain_severity_hm3` | Total volume-based severity (hm³) |

---

## Sanity Checks

After the lag filter, the notebook runs four automated checks:

1. All retained `lag_days` are ≤ 0.
2. All `overlap_days` are ≥ `MIN_OVERLAP`.
3. All `propagation_fraction` values are in [0, 1].
4. No duplicate `(origin_station, origin_event_id, upstream_station)` triples.

---

## File Structure

```
05_ssi_propagation.ipynb                 ← Main notebook (this algorithm)
05b_ssi_propagation_vectorised.ipynb     ← Vectorised version, identical results
04_ssi_drought_events.ipynb              ← Upstream dependency: produces the event catalogue
data/SSI_drought_events.csv              ← Input: event catalogue
data/upstream_connectivity.csv           ← Input: network topology
data/propagation_W45_minov5.csv          ← Output: all pairs (before lag filter)
data/propagation_W45_minov5_lagneg.csv   ← Output: lag-filtered pairs + chain metrics
data/chains_summary_W45_minov5.csv       ← PRIMARY OUTPUT: one chain per origin event
data/Chains_final.csv                    ← Alias of chains_summary for downstream scripts
```

---

## Quick Reference: Key Decisions

| Decision | Rationale |
|---|---|
| Asymmetric window: `[origin_start − 45d, origin_start]` | Left expansion captures precursor upstream droughts; right boundary capped at `origin_start` so the window is exactly consistent with the lag constraint — no upstream event starting after `origin_start` is even evaluated |
| Require ≥5 days real overlap | Eliminates spurious near-misses while retaining genuine co-occurring events |
| Keep only lag ≤ 0 | Enforces physical causality: upstream must precede downstream |
| Normalise by `n_upstream_with_events` | Makes propagation fraction comparable across basins with different network sizes |
| No chain merging | Each origin event is a distinct, continuous drought episode. Merging at the chain level conflates temporally distant events and violates origin event integrity. See [rationale](#no-merge-rationale) |
