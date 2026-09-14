# Events that do not propagate — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #2, comment 3**: *"assess (a) upstream drought events that did not propagate downstream; and (b) downstream drought events that occurred without any upstream trigger."*

**Headline result.** Both directions show coherent, physically sensible patterns. **(b)** 403 of 928 origin events (43.4%) have no upstream connection; these isolated events are only moderately less severe than connected ones, and the probability of isolation is driven mainly by **network position** (0% at the best-connected origin station, 78% at the most limited one) rather than by event characteristics. **(a)** Of 1,558 upstream events that were genuine candidates for at least one origin event, 602 (38.6%) never propagate; these non-propagating events are roughly a third the duration and a quarter the severity of events that do propagate — evidence that the matching algorithm correctly discriminates substantial, sustained droughts as propagation drivers rather than selecting arbitrarily.

---

## 1. Rationale

The published catalogue of 525 chains (Figure 5) characterises *connected* events only. Reviewer #2 asks for the complementary picture: what happens to events that do **not** connect, in either direction? This matters for interpretation because it tests whether the matching algorithm behaves sensibly at the margins — whether isolated downstream events are simply weak/negligible, and whether non-propagating upstream events are arbitrary or systematically different from those that do propagate.

## 2. Inputs

Same as the R3.2 null-test analysis: `data/SSI_drought_events.csv` (2,625 events, 33 stations) and `data/upstream_connectivity.csv`.

## 3. Method

### 3.1 Reuse, don't reconstruct
`propagation_null_core.py` — the exact vectorised matching algorithm validated in the R3.2 null-test analysis (confirmed there to reproduce the published headline numbers and Figure 5 statistics to three decimals) — is reused **without modification** to its accepted-match logic. Two diagnostic passes are added on top.

### 3.2 Part (b): isolated downstream events
For every origin event (at the 12 stations with upstream connectivity), record whether it received at least one accepted upstream match. Unconnected origin events are the direct complement of the published 56.6% connected fraction.

### 3.3 Part (a): non-propagating upstream events
The relevant population is not "any event at a station that serves as upstream for some origin" — many such events are temporally distant from any origin event and were never genuinely in contention. Instead, for every origin event and every connected upstream station, the same window filter used by the published algorithm (`end >= origin_start - W`, `start <= origin_start`) is applied to identify every upstream event that was ever a **genuine candidate**. Among this population of 1,558 candidate events, those never selected as the accepted match for any origin (whether because their overlap fell below M, or because they lost the deterministic tie-break to a stronger candidate at the same station) are classified as non-propagating.

## 4. Results

### 4.1 Isolated downstream events (Table S3)

| | Connected (n=525) | Unconnected (n=403) |
|---|---|---|
| Mean duration (days) | 35.1 | 32.7 |
| Mean severity | 13.88 | 11.22 |

Unconnected fraction by origin station ranges from **0%** (station 9011, 22 upstream contributors) to **78.0%** (station 9004). The regulated station (9101, downstream of the Yesa reservoir) is the third most isolated (74.3%), consistent with regulation-driven decoupling discussed elsewhere in the revision (response to Reviewer #4, Comment 3).

### 4.2 Non-propagating upstream events

| | Propagates (n=956) | Does not propagate (n=602) |
|---|---|---|
| Mean duration (days) | 54.3 | 19.3 |
| Mean severity | 24.10 | 5.64 |

Events that propagate are ~2.8× longer and ~4.3× more severe, on average, than those that do not.

## 5. Interpretation

**Isolated downstream events are not simply negligible.** Their severity and duration are only moderately lower than connected events, meaning many represent genuine local droughts that lack a detectable upstream driver — most likely because they occur at stations with limited upstream network complexity, rather than because they are trivial fluctuations. This is directly supported by the strong station-to-station variation in isolation rate, which tracks network position (number and configuration of upstream contributors) far more than it tracks any property of the events themselves.

**Non-propagating upstream events behave exactly as expected physically.** The ~3–4-fold gap in duration and severity between propagating and non-propagating events indicates the algorithm is not selecting propagation partners arbitrarily: it preferentially captures substantial, sustained droughts as chain drivers, consistent with the physical expectation that only sufficiently prolonged or severe deficits generate a signal detectable downstream. This is a reassuring validation of the matching algorithm's behaviour at the margins, complementing the network-permutation null test (Section S4) and the lag-vs-distance analysis (Section S3).

## 6. Reproducibility

**Dependencies:** Python ≥ 3.9 with `pandas`, `numpy`, `matplotlib`.

**Files:**
- `../propagation_null_core.py` — the validated matching algorithm (unchanged from the R3.2 analysis), shared module kept once in `revision_analyses/`.
- `01_non_propagating_events.py` — standalone script: validates against published numbers, computes both diagnostic passes, saves the figure.
- `02_non_propagating_events.ipynb` — the same pipeline as a documented, executed notebook.

**Run:**
```bash
python 01_non_propagating_events.py
```
Place `SSI_drought_events.csv` and `upstream_connectivity.csv` in a `data/` subfolder alongside the script/notebook.

**Outputs:** `origin_events_connectivity_R2_3.csv` (928 origin events with connectivity status), `upstream_candidate_events_R2_3.csv` (1,558 candidate upstream events with propagation status), `table_S3_unconnected_by_station.csv`, `fig_non_propagating_events.png` / `.pdf` (five-panel figure).

## 7. Caveats of the analysis itself

- "Non-propagating" (part a) is defined relative to events that were genuine window-filtered candidates, not all events at upstream-role stations; this is a deliberate, more precise choice (see Section 3.3), but means the 38.6% figure should not be read as "38.6% of all upstream droughts never contribute to any chain" — only of those that were ever temporally in contention.
- An upstream event classified as "non-propagating" may still have lost a close tie-break to a stronger candidate at the same station for a specific origin event, rather than having no relationship to the origin event at all; the diagnostic does not distinguish "never a plausible candidate" from "always outcompeted."
- Station-level isolation rates (Table S3) reflect the specific, fixed network topology of the Ebro basin and are not necessarily generalisable to other basins with different upstream connectivity structures.

## 8. Audit log

- ✅ Reused, unmodified matching algorithm re-validated to reproduce the published 525/928 (56.6%) connected origin events before running any diagnostic pass.
- ✅ The "candidate" population for part (a) is explicitly restricted to genuine window-filtered candidates, not all events at upstream-role stations, avoiding an inflated or misleading denominator.
- ✅ Notebook executed end-to-end with a real kernel (zero errors after fixing an escaping typo in the generation script).
- ✅ Figure colour channels verified programmatically (not just "non-blank") to confirm both series are actually rendered.
- ✅ Both directions (isolated downstream, non-propagating upstream) reported with full distributional statistics, not just headline percentages, and the station-level breakdown is reported in full (Table S3) rather than only summarised.
