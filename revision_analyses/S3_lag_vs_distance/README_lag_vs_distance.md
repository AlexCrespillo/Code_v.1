# Propagation lag vs. along-channel distance — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #4, comment 2** ("insufficient physical support; incorporate flow travel time and hydraulic distance") and basis for reframing the *no-flow-routing* limitation (**comment 16**). It is a methodological inflection point: the result decides whether a travel-time/routing constraint should be imposed on the event-matching, or whether the deliberately non-routing design is justified.

**Headline result.** Across matched upstream–origin pairs, the propagation lag is **at most weakly, and not significantly, related to along-channel distance** (station-pair Spearman ρ = 0.14, p = 0.16; ρ = 0.16, p = 0.12 excluding the regulated station). Distance explains ≈1–2 % of the variance in lag. Drought-propagation timing in this basin is therefore controlled by catchment **storage, memory and regulation**, not by hydraulic travel time — consistent with the manuscript's definition of propagation as a spatiotemporal (not flow-routing) process.

---

## 1. Rationale

For flood waves, the downstream lag scales with hydraulic distance through wave celerity, so a travel-time constraint is physically meaningful. Hydrological drought is different: the downstream signal reflects the **depletion and recovery of stores** (soil, aquifers, reservoirs), so the onset lag need not track distance at all. Testing this directly (a) answers the reviewer's request for physical support with data rather than assertion, and (b) tells us whether imposing a routing/travel-time filter would be appropriate. If lag does **not** scale with distance, such a filter would force a relationship the data do not support and would be inappropriate for drought (as opposed to flood) propagation.

## 2. Inputs

| File | Content | Notes |
|---|---|---|
| `Locations.csv` | `station_id, xutm, yutm, NEAR_FID` | `;`-separated, decimal comma; UTM metres; `NEAR_FID` = nearest network segment from the ArcGIS *Near* tool |
| `Grafo.csv` | `arcid, from_node, to_node, ORIG_FID, Shape_Length` | node–arc list exported from `river_connectivity_singlepart`; lengths in metres, decimal comma |
| `Propagation_raw.csv` | candidate upstream→origin **pair-events** | includes `origin_event_id, upstream_event_id, overlap_days, upstream_severity, lag_days`. **Pre-tie-break** (see §3.3) |

`lag_days` is defined as `upstream_start − origin_start ≤ 0`; `|lag|` is the number of days the upstream event began *before* the origin event.

## 3. Method and justification, step by step

### 3.1 River network as an undirected weighted graph
Nodes are network junctions; edges are river segments weighted by their length (m). Parallel edges (rare) keep the shortest length. **Correctness guarantee:** the graph has 978 nodes, 974 edges and 4 connected components, i.e. `edges = nodes − components` → the graph is a **forest** (`networkx.is_forest = True`). On a forest there is **exactly one path** between any two connected nodes, so the shortest path *is* the unique along-channel path and flow direction is irrelevant for distance. All 101 station pairs are reachable (no `NaN` distances; enforced by an assertion in the code).

### 3.2 Station → node mapping, and its validation
Each station is mapped to a network node via `NEAR_FID → ORIG_FID` (0-based), i.e. the station's nearest segment; the station is assigned to that segment's downstream endpoint (`to_node`).

Two independent checks confirm the mapping is correct:
- **Cross-check against straight-line distance:** along-channel distance correlates with Euclidean distance at **r = 0.97**, with a median network/Euclidean ratio of **1.65** (river paths ~65 % longer than the straight line — physically sensible sinuosity). A wrong key would not produce r = 0.97.
- **Node-endpoint sensitivity:** assigning stations to `to_node` vs `from_node` changes the station-pair correlation only from ρ = 0.14 to ρ = 0.12 (both non-significant), with corr(network, Euclidean) = 0.97 in both cases. The endpoint choice is immaterial to the conclusion.

The correct key was determined empirically by testing candidate keys (`arcid`, `ORIG_FID`, row index): only `ORIG_FID` (= row index) gave r = 0.97; `arcid` gave r = 0.40 and was rejected.

### 3.3 Tie-break to the published matched set (critical correction)
`Propagation_raw.csv` is the **pre-tie-break candidate set**: a single origin event can have several candidate upstream *events* at the *same* upstream station (here 466 such cases, up to 4 candidates each; 2,075 candidate rows in total). The published methodology (Section 2.1.3) retains **one event per upstream station** through a deterministic ranking:

1. maximum temporal overlap with the origin event;
2. highest drought severity;
3. minimum onset difference (= minimum `|lag|`);
4. event identifier as final deterministic tie-break.

This ranking is reproduced in the code, yielding **1,512 matched pair-events across 101 station pairs**, consistent with the 525 reconstructed chains. **Why it matters:** using the raw candidate set instead of the tie-broken set is inconsistent with the chains and changes the numbers — the station-pair correlation moves from ρ = 0.06 (raw) to ρ = 0.14 (correct). The qualitative conclusion (no significant lag–distance relationship) is unchanged, but all reported statistics use the tie-broken set.

### 3.4 Distance metrics
Two metrics were computed per **station pair** (distance is a property of the pair, not the event): the straight-line (Euclidean) distance in UTM, used only as a first pass and mapping validation, and the **along-channel distance** (shortest path on the forest), which is the metric the reviewer asked for and the one used for all reported results.

### 3.5 Two analysis levels (and why both)
- **Event level (n = 1,512):** every matched pair-event is a point. High power but **pseudoreplicated** — the many events of a given station pair share one distance and are not independent — so its p-value is optimistic.
- **Station-pair level (n = 101):** one point per pair, using the **median `|lag|`** of its events. This removes the pseudoreplication and is the level on which the conclusion is based.

Because distance is constant within a station pair, the event-level panel shows vertical bands (same x, many y); this is expected, not an artefact.

### 3.6 Within-origin partial correlation
To separate a genuine distance effect from between-sub-basin heterogeneity, `|lag|` and distance were de-meaned by origin station and re-correlated (a within-origin/fixed-effects estimator), and per-origin Spearman coefficients were computed. This isolates whether, *for a given downstream station*, more distant upstream contributors show longer lags.

### 3.7 Regulation
Station 9101 (immediately downstream of the Yesa reservoir, documented in the manuscript) is flagged as regulated. This flag is **illustrative, not a complete regulation inventory**; it is used only to show that regulated pairs are identifiable outliers and to check that the main result is not driven by them.

## 4. Results (tie-broken set, along-channel distance)

| Test | n | Spearman ρ | p | Reading |
|---|---|---|---|---|
| Event level | 1,512 | +0.09 | 4×10⁻⁴ | significant only through pseudoreplication; effect ≈ 1 % |
| **Station-pair level** | **101** | **+0.14** | **0.16** | **not significant** |
| Pair level, excl. 9101 | 97 | +0.16 | 0.12 | not significant; result not driven by regulation |
| Within-origin partial | 1,512 | +0.10 | 7×10⁻⁵ | weak but significant |
| Per-origin (6 origins) | — | median +0.09 | — | 2/6 significant (9011 ρ=0.22; 9002 ρ=0.12); others ≈0 or negative |
| Regulation signature | — | — | — | mean `|lag|` at 9101 = 175 d vs 32 d elsewhere (**5.6×**) |

## 5. Interpretation and conclusion

The propagation lag is, at most, **weakly and non-significantly** related to along-channel distance (~1–2 % of variance; not significant at the station-pair level, with or without the regulated station). A faint within-origin tendency exists, concentrated in two origins (notably 9011), but it is far too small to justify a travel-time constraint.

This is physically informative rather than a null: unlike flood propagation, hydrological drought propagation is not governed by flow celerity but by store depletion and recovery, so the lag reflects **storage, catchment memory and regulation** rather than distance. Reservoir regulation is the clearest single control identified (station 9101, 5.6× longer lags). Imposing a hydraulic-distance/travel-time constraint would therefore force a relationship the data do not support and would be inappropriate for drought propagation.

**Consequences for the manuscript**
- **R4.2:** the absence of a lag–distance relationship *justifies* the non-routing design; a supplementary figure and a main-text sentence report it.
- **R4.16:** the routing limitation is reframed as a deliberate, physically motivated simplification; the natural future extension is storage-/memory-based timing, not hydraulic routing.
- **Important caveat for the response:** the physical validity of the chains does **not** rest on a lag–distance relationship, but on directional connectivity, temporal ordering (`lag ≤ 0`), and the surrogate/null test (Reviewer #3, comment 2). That test remains the primary evidence that chains exceed chance co-occurrence.

## 6. Reproducibility

**Dependencies:** Python ≥3.9 with `pandas`, `numpy`, `networkx`, `scipy`, `matplotlib`.
**Run:** place the three input CSVs beside the script/notebook, set `UP`/`OUT`, then

```bash
python lag_distance_analysis.py
```

or run `lag_vs_distance_analysis.ipynb` top-to-bottom.

**Outputs:** `fig_lag_vs_distance.png` / `.pdf` (the supplementary figure), `pairwise_lag_distance.csv` (the 101 station pairs with distances and median lag), and the printed statistics table.

## 7. Caveats of the analysis itself

- Stations are assigned to a segment endpoint rather than to the exact snapped point on the segment; the residual error is bounded by segment length and is immaterial here (endpoint-choice sensitivity: ρ = 0.12–0.14).
- The regulated-station set is illustrative (9101/Yesa); a full regulation inventory is out of scope, and the main result does not depend on it.
- The within-origin de-meaning combines a linear operation with a rank correlation; it is reported as a robustness indicator alongside the per-origin Spearman coefficients, not as a primary statistic.

## 8. Audit log

- ✅ Graph verified as a **forest** → unique along-channel paths; all pairs reachable.
- ✅ `NEAR_FID → ORIG_FID` mapping validated (r = 0.97 vs Euclidean) and endpoint-choice tested.
- ✅ **Tie-break bug found and fixed:** analysis now uses the post-tie-break set (1,512 pairs) consistent with the 525 chains; statistics updated (station-pair ρ 0.06 → 0.14; conclusion unchanged).
- ✅ Distances checked for `NaN` (none; assertion enforced).
- ✅ Result robust to removing the regulated station and to the node-endpoint choice.
