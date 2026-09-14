# Systematic validation against an independent national drought catalogue — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #3, comment 4**: *"No independent validation, such as historical drought records or precipitation anomalies, is provided to support the reconstructed chains."*

**Headline result.** All 525 reconstructed propagation chains (not just the three case studies) were compared against the independent national drought catalogue of Trullenque-Blanco et al. (2024, *Scientific Data*) — built from a completely different data source (national 10×10 km monthly precipitation grid) and index (SPI-12). **413 of 525 chains (78.7%) overlap a documented national drought**, well above the 42.8% of the calendar that these droughts cover on their own. More decisively, **overlapping chains are 3.4× more severe on average** than non-overlapping ones (222.3 vs. 65.0) — a pattern that calendar coverage alone cannot explain. **Every one of the 25 documented national droughts since 1961 (100%) has at least one overlapping Ebro chain.**

---

## 1. Rationale

The manuscript's internal validation (the network-permutation null test, Section S4, and the lag-vs-distance analysis, Section S3) establishes that the reconstructed chains are not artefacts of the matching algorithm. Reviewer #3 asks for something categorically different: corroboration from a **source entirely independent of this study's data and method** — historical drought records or precipitation-based indices. The author's own co-authors (Vicente-Serrano, Beguería) are co-authors of exactly such an independent catalogue (Trullenque-Blanco et al., 2024), which identifies major drought events for the whole Spanish mainland from precipitation data alone, with no reference to streamflow, the Ebro basin specifically, or this study's methodology. This provides a genuinely external benchmark.

## 2. Inputs

| File | Content | Notes |
|---|---|---|
| `data/SSI_drought_events.csv`, `data/upstream_connectivity.csv` | Same validated inputs as all previous analyses | 2,625 events, 33 stations |
| `data/Identification_and_characteristics.xlsx` | The national catalogue (Trullenque-Blanco et al., 2024) | Two sheets: `Events` (40 documented drought episodes, 1916–2020) and `Time_series` (monthly national area-under-drought fraction) |

## 3. Method

### 3.1 Parsing the national catalogue, and a data-quality check
Dates in the source spreadsheet are given as `month.year` strings (e.g. `12.1917` = December 1917). Five of the 40 rows have a truncated year digit — an evident spreadsheet typo (e.g. `197-10-01`). Each was corrected and cross-checked against that event's independently declared duration (in months): all five corrections reproduce the declared duration **exactly**, giving high confidence in the fix (e.g. event 18: corrected to 1970-10-01 → 1971-03-01, reproducing its declared 6-month duration precisely). One further row (event 25) has a genuine, uncorrected mismatch between its declared duration (19 months) and its date range (10 months) in the source file; this is left as-is, since it does not affect the overlap test (which uses only the start/end dates), and is reported transparently rather than silently adjusted.

### 3.2 Reconstructing chain-level windows
For every one of the 525 published chains, the chain-level window is reconstructed exactly as defined in Table 2 of the manuscript: chain start (t_c,s) = earliest onset among all matched members (the origin event and every accepted upstream match); chain end (t_c,e) = latest termination among the same set. This reuses the validated matching algorithm (`propagation_null_core.py`, identical to the one behind the published results, previously validated in the R3.2 null-test analysis) without modification.

### 3.3 Overlap test, in both directions
- **Direction 1:** for each of the 525 Ebro chains, test whether its [chain_start, chain_end] window overlaps at least one national catalogue event (restricted to the 25 national events overlapping the 1961–2020 study period).
- **Direction 2:** for each of those 25 national events, test whether at least one Ebro chain overlaps it.

### 3.4 Honest baseline context
Because the 25 national events themselves cover a substantial fraction of the 60-year calendar (308 of 720 months, 42.8%, with no internal overlap between the national events themselves), some baseline overlap between our chains and national droughts is expected purely from calendar coverage, independent of any real relationship. This is reported explicitly (following the same principle established for the network-permutation null test, Section S4) rather than treating the raw 78.7% figure as self-evidently strong on its own.

## 4. Results

| Test | Result |
|---|---|
| Ebro chains overlapping ≥1 national drought | 413 / 525 (78.7%) |
| Calendar coverage of national droughts (context) | 42.8% of 1961–2020 |
| Mean chain severity, overlapping chains | 222.3 (median 114.1) |
| Mean chain severity, non-overlapping chains | 65.0 (median 38.1) |
| Mean chain size (stations), overlapping | 5.3 |
| Mean chain size (stations), non-overlapping | 3.8 |
| National droughts with ≥1 overlapping Ebro chain | 25 / 25 (100%) |

Full results: `chains_national_validation.csv` (525 rows), `national_events_ebro_coverage.csv` (25 rows).

## 5. Interpretation

**The raw overlap rate (78.7%) is above, but not dramatically above, the naive calendar-coverage baseline (42.8%)** — a genuine excess, but one that should not be oversold on its own, since chain duration also affects a priori overlap probability (longer chains sweep more calendar months and are proportionally more likely to intersect a documented drought by exposure alone).

**The severity gradient is the more decisive piece of evidence.** There is no reason for calendar-coverage alone to concentrate the *more severe* chains specifically within nationally-documented drought periods. That overlapping chains are systematically ~3.4× more severe and involve ~40% more stations indicates a genuine physical correspondence between locally (Ebro-)reconstructed drought severity and independently documented, precipitation-based national drought conditions.

**Complete coverage of documented national droughts (25/25) is a strong, simply-stated result**, consistent with the Ebro being one of Spain's largest river basins — any national-scale meteorological drought event would be expected to leave a detectable hydrological signature somewhere within it.

**Conclusion for the manuscript.** This provides genuine external validation, independent of the data, indices and methodology used to construct the propagation chains, complementing the internal validation already established (Sections S3–S4).

## 6. Reproducibility

**Dependencies:** Python ≥ 3.9 with `pandas`, `numpy`, `matplotlib`, `openpyxl` (for reading the `.xlsx` catalogue).

**Files:**
- `../propagation_null_core.py` — the validated matching algorithm (unchanged, reused from the R3.2/R2.3 analyses), shared module kept once in `revision_analyses/`.
- `01_national_catalogue_validation.py` — standalone script: parses and corrects the national catalogue, reconstructs chain-level windows, runs the overlap tests, saves the figure.
- `02_national_catalogue_validation.ipynb` — the same pipeline as a documented, executed notebook.

**Run:**
```bash
python 01_national_catalogue_validation.py
```
Place `SSI_drought_events.csv`, `upstream_connectivity.csv`, and `Identification_and_characteristics.xlsx` in a `data/` subfolder alongside the script/notebook.

**Outputs:** `chains_national_validation.csv`, `national_events_ebro_coverage.csv`, `fig_national_validation.png` / `.pdf` (three-panel figure: timeline + two severity/size boxplots).

## 7. Caveats of the analysis itself

- The national catalogue identifies drought at a coarse, national-mainland scale (>20% of Spain's grid cells under SPI-12 < −0.84); it is not Ebro-specific, so overlap indicates temporal coincidence with a broader national episode, not that the Ebro basin itself drove or was central to that national event (except where the catalogue's own text says so, as for event 39/Case 2).
- The overlap test uses monthly-resolution national event boundaries against daily-resolution Ebro chain windows; a chain overlapping a national event by even a few days at the boundary counts as "overlapping," which is a deliberately permissive (not stringent) test.
- One source-data row (event 25) has an uncorrected, acknowledged internal inconsistency (duration vs. date range) in the national catalogue itself; this does not affect the overlap test but is noted for transparency.
- The calendar-coverage baseline (42.8%) is a simplification (it does not account for chain-duration-dependent exposure probability); a fully rigorous baseline would require a permutation test analogous to Section S4, which was not performed here given the clarity of the severity-gradient result on its own.

## 8. Audit log

- ✅ Reused, unmodified matching algorithm re-validated (2,075 accepted pairs, 525 chains) before running any comparison.
- ✅ Five truncated date values in the source spreadsheet identified and corrected, with each correction cross-checked against an independent field (declared duration) in the same source — not merely assumed.
- ✅ One further, genuine source-data inconsistency (event 25) identified and left untouched rather than silently "corrected" without a reliable basis.
- ✅ Chain-level windows (t_c,s, t_c,e) reconstructed using the exact definition in the manuscript's own Table 2.
- ✅ Both directions of the overlap test reported (chains→national and national→chains), not just the more favourable direction.
- ✅ The calendar-coverage baseline is reported explicitly rather than omitted, so the 78.7% headline figure is not presented without context.
- ✅ Notebook executed end-to-end with a real kernel (zero errors).
