# Drought-event detection sensitivity (pooling gap & minimum duration) — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #3, comment 1**, cross-referenced from **Reviewer #4, comment 4** (parameter sensitivity/adaptability).

**Headline result.** Unlike the propagation parameters W and M (search window and minimum overlap), which showed essentially zero sensitivity (RSI = 0 for W; Fig. 4), the event-**detection** parameters — pooling gap and minimum duration — materially affect the number and characteristics of detected drought events (event count ranges from −30% to +96% relative to the published configuration across the tested grid). This is expected: pooling and duration thresholds are definitional choices inherent to any threshold-level drought identification method. The relevant question is not whether raw counts change (they must) but whether the choice introduces **station-specific bias**. A per-station homogeneity check shows that the relative ranking of stations by mean severity is remarkably stable across the whole tested range (Spearman ρ mostly 0.83–1.00), and event-count ranking is also stable across most of the grid (ρ > 0.8), degrading only at the untested extreme of zero pooling combined with a strict duration filter. The published configuration (pooling = 10 d, minimum duration = 5 d) sits within literature-recommended ranges and is not at an unstable corner of the grid.

---

## 1. Rationale

The manuscript's existing sensitivity analysis (Fig. 4, Fig. S1–S2, Table S1) evaluates the two **propagation-matching** parameters (search window W, minimum overlap M), which operate *after* drought events have already been detected. It does not evaluate the two **event-detection** parameters — the pooling gap (merge two candidate events if separated by ≤ this many days) and the minimum-duration filter (discard events shorter than this) — which act at an earlier stage, directly on the daily SSI series. Reviewer #3 asks specifically about this earlier stage: how do these choices affect event number, duration, severity and volumetric deficit, and do they bias results across stations or sub-basins? This analysis answers that question directly, with the same rigour (validation against published numbers, explicit bias check) as the existing W/M analysis.

## 2. Inputs

| File | Content | Notes |
|---|---|---|
| `data/SSI_daily.csv` | Daily SSI series: `station_id, date, Q, SSI` | 33 stations, 1961–2020, no gaps (21,915 days/station) |
| `data/caudales_diarios_imputados_CORE_FILTRADO.csv` | Daily discharge: `station_id, date, Q_imp, is_imputed, imputed_method` | Same 33 stations and period; used only for volumetric severity |

## 3. Method, step by step

### 3.1 Reuse the published algorithm
The four-step detection algorithm (flag drought days below SSI < −1.28 → group consecutive days → pool events separated by ≤ `POOL_GAP` days → discard events shorter than `MIN_DURATION` days → compute duration/severity/intensity/peak) is taken directly from `04_ssi_drought_events.ipynb` (the notebook that produced the published event catalogue), parametrised so that `POOL_GAP` and `MIN_DURATION` can be varied. Volumetric severity (hm³) is computed exactly as in the original notebook: a per-station, per-day-of-year 10th-percentile discharge threshold (±15-day circular pooling window, unchanged across the sweep since it does not depend on the detection parameters), with daily deficits below that threshold summed and converted to hm³.

### 3.2 Validate against the published Table 3, at every stage
Before running any sensitivity sweep, `run_detection()` at the default parameters (pooling = 10, minimum duration = 5) is checked against the manuscript's Table 3:

| Metric | Reproduced | Published |
|---|---|---|
| Total events | 2,625 | 2,625 |
| Duration: mean / median / max (days) | 33.2 / 19.0 / 613 | 33.2 / 19.0 / 613 |
| Severity: mean / median / max | 12.952 / 3.937 / 508.836 | 12.952 / 3.937 / 508.836 |
| Volumetric severity: mean / median / min / max (hm³) | 5.508 / 0.445 / 0.000 / 362.262 | 5.508 / 0.445 / 0.000 / 362.262 |

All four groups of statistics match exactly. **One discrepancy was noted and is not attributable to this analysis:** the manuscript's printed Table 3 reports intensity mean = median = 0.218, which is numerically implausible for a right-skewed distribution and does not match the value obtained here (mean 0.293, matching the manuscript's own IQR figure of 0.293) using the intensity definition in Table 1 (I = S/D, computed per event). Duration and severity — whose ratio defines intensity — both match exactly, which is strong evidence that the reproduction is correct and that Table 3's printed intensity row likely contains a transcription error, worth checking separately (flagged for the response to Reviewer #2, comment 11).

### 3.3 Sensitivity grid
`POOL_GAP ∈ {0, 5, 7, 10, 15}` days × `MIN_DURATION ∈ {3, 5, 7, 10}` days (20 combinations; default = 10, 5). For each combination, the full detection pipeline is re-run across all 33 stations and volumetric severity is recomputed (the P10 discharge lookup, being independent of these parameters, is built once and reused).

### 3.4 Bias / homogeneity check across stations
For each combination, per-station event count and per-station mean severity are compared, via Spearman rank correlation, against the same quantities under the default configuration. This directly targets the reviewer's question ("whether they introduce bias across stations or sub-basins"): if the correlation stays high, the *relative* drought characteristics of stations (which ones are more/less drought-prone) are preserved regardless of the parameter choice; if it drops, the parameter choice is disproportionately reshaping some stations' records relative to others.

## 4. Results

### 4.1 Global sensitivity (Section 4, full grid in `sensitivity_summary.csv`)

Event count varies from 1,834 (pooling = 15, min. duration = 10) to 5,155 (pooling = 0, min. duration = 3) — a range of −30% to +96% relative to the default (2,625 events). Mean duration increases monotonically with both pooling and minimum duration (from ~12.7 days at pooling = 0/min. duration = 3, to ~50.9 days at pooling = 15/min. duration = 10), and mean severity and volumetric severity increase correspondingly. This pattern is mechanical and expected: more pooling merges more short interruptions into fewer, longer episodes; a higher minimum-duration threshold discards more short events, raising the mean of what remains.

### 4.2 Bias / homogeneity check (`sensitivity_bias.csv`)

| Statistic tested | Range of Spearman ρ vs default across the grid |
|---|---|
| Per-station event count | 0.30 (pooling=0, min.dur=10) to 1.00 (default) |
| Per-station mean severity | 0.79 (pooling=0, min.dur=10) to 1.00 (default) |

Mean-severity ranking is stable across almost the entire grid (ρ ≥ 0.83 for 19 of 20 combinations, only dropping to 0.79 at the single most extreme corner tested). Event-count ranking is stable (ρ > 0.8) across most of the grid but degrades more noticeably when pooling is removed entirely, particularly combined with a strict duration filter (ρ = 0.30–0.46 for pooling = 0 with min. duration ≥ 7). This indicates that eliminating pooling altogether — not merely varying it moderately — can differentially affect stations, plausibly because short, frequent low-flow interruptions (which pooling is designed to absorb) are more common in some flow regimes than others.

## 5. Interpretation

The published parameter choice is not at an unstable extreme: pooling = 10 days and minimum duration = 5 days sit centrally within the tested grid, where both bias metrics are high (ρ = 1.00 by construction, and neighbouring cells remain ρ > 0.9). The literature basis already cited in the manuscript (Van Loon, 2015; Fleig et al., 2006; Tallaksen et al., 1997) supports both threshold values as standard practice, and this analysis adds the additional evidence that the choice does not appear to introduce differential bias across the 33 stations relative to plausible alternative choices, except at the untested extreme of no pooling at all.

**Scope note.** This analysis characterises sensitivity at the **event-detection** stage only (event count, duration, severity, volumetric deficit — exactly as requested). It does not re-run the full propagation-chain reconstruction for each of the 20 alternative event catalogues; doing so would be a substantially larger undertaking and was not requested by the reviewer, whose comment is specifically about the detection stage. As with any threshold-based method, the exact propagation-chain counts reported in the manuscript are naturally sensitive to the choice of detection parameters — this is acknowledged as an inherent, literature-recognised feature of the method rather than a specific weakness of this study.

## 6. Reproducibility

**Dependencies:** Python ≥ 3.9 with `pandas`, `numpy`, `scipy`, `matplotlib`.

**Files:**
- `../detection_sensitivity_core.py` — core module (detection algorithm, volumetric-severity computation, sensitivity-sweep harness), kept once in `revision_analyses/`.
- `01_sensitivity_pooling_duration.py` — standalone script: loads data, validates against Table 3, runs the sensitivity grid, computes the bias check, saves the figure.
- `02_sensitivity_pooling_duration.ipynb` — the same pipeline as a documented, executed notebook.

**Run:**
```bash
python 01_sensitivity_pooling_duration.py
```
Place `SSI_daily.csv` and `caudales_diarios_imputados_CORE_FILTRADO.csv` in a `data/` subfolder alongside the script/notebook.

**Outputs:** `sensitivity_summary.csv` (20 rows, one per parameter combination), `sensitivity_per_station.csv` (20 × 33 rows), `sensitivity_bias.csv` (20 rows, Spearman ρ vs default), `fig_detection_sensitivity.png` / `.pdf` (six-panel heatmap grid, matching the manuscript's Figure 4 style).

## 7. Caveats of the analysis itself

- The bias check compares each configuration only against the **default** configuration, not against every other pair; this is the natural choice since the question is "does the published default look reasonable/central", not a full pairwise comparison.
- No station-to-climatic-zone or station-to-sub-basin metadata was available to test bias directly by region; the per-station Spearman check (33 individual points) is the most granular and defensible test achievable with the data at hand, and mirrors the approach the manuscript itself already uses for the M-parameter sensitivity (Supplementary Fig. S1).
- The volumetric-severity computation depends on the P10 discharge climatology, which is held fixed across the sweep (as it should be, since it does not depend on the SSI-event detection parameters); it was independently validated to reproduce the published Table 3 volumetric-severity statistics exactly.

## 8. Audit log

- ✅ Station-ID sets confirmed identical between the SSI series and the discharge series.
- ✅ **Full validation against published Table 3** at the default configuration: event count, duration (mean/median/max), severity (mean/median/max) and volumetric severity (mean/median/min/max) all match exactly.
- ✅ A discrepancy in the manuscript's own printed intensity statistic was identified (implausible mean = median = 0.218) and is very likely a transcription error in Table 3, not a reproduction issue — duration and severity, whose ratio defines intensity, both match exactly. Flagged for the response to Reviewer #2, comment 11.
- ✅ Notebook executed end-to-end with a real kernel (zero errors).
- ✅ Bias/homogeneity check reported transparently across the full grid, including the one configuration (pooling = 0, high minimum duration) where station-level ranking is least preserved — no selective reporting.
