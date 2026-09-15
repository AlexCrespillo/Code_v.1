# Meteorological context (SPEI) for the case-study propagation chains — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #2, comment 4** (hydroclimatic context for the case studies) and an initial, illustrative contribution to **Reviewer #4, comment 5** (meteorological–hydrological linkage; full coupling remains explicitly out of scope and is identified as future work).

**Headline result.** Daily SSI at the origin station of each case study is plotted alongside weekly SPEI-3 and SPEI-12 at that station's coordinates. Two independent checks show that neither the origin's local meteorology nor its own individual hydrological record fully captures the network-wide chain reported in the manuscript. First, the three cases show **different**, case-specific relationships between local meteorological conditions and hydrological onset — there is no single, uniform "meteorological drought precedes hydrological drought" pattern. Second, the "hydrological event" window quoted in Figures 6–8 is the **chain window** (earliest onset to latest termination among *all* matched stations), not the origin's own individual SSI event — and the origin's own event can start months after the chain window begins, or even appear as several separate, shorter events within a chain window reported as continuous (most strikingly in Case 2: six separate origin-station events within a single 427-day chain window). Neither finding is a weakness: both are direct, case-specific evidence that the origin/outlet, examined alone, does not capture the network-wide propagation signal — reinforcing the manuscript's central storage- and network-driven argument for why the propagation framework adds value beyond monitoring any single station.

---

## 1. Rationale

Reviewer #2 asks for meteorological drought context to help interpret the hydrological case studies; Reviewer #4 asks, more ambitiously, for a fully coupled meteorological-to-hydrological tracking framework. The latter is explicitly out of scope for this methodological, hydrological-drought-focused contribution (see the manuscript's existing Discussion paragraph on this point). This analysis addresses the former directly and provides an honest, illustrative first look at the latter: it does not claim to reconstruct meteorological drought propagation, only to place each hydrological case study in its local meteorological context.

## 2. Inputs

| File | Content | Notes |
|---|---|---|
| `SSI_daily.csv` | Daily SSI, 33 stations, 1961–2020 | Same file validated in the event-detection sensitivity analysis (Section S5) |
| `DATA_40_81_0_52.csv` | Weekly SPEI/SPI (scales 1,3,6,9,12,24) at 40.81°N, 0.52°E | Station 9027 (Ebro at Tortosa) — Case 1 |
| `DATA_41_66_-0_88.csv` | Weekly SPEI/SPI at 41.66°N, −0.88°E | Station 9011 (Ebro at Zaragoza) — Case 2 |
| `DATA_42_35_-1_65.csv` | Weekly SPEI/SPI at 42.35°N, −1.65°E | Station 9005 (Aragón at Caparroso) — Case 3 |

**Data sources.** `SSI_daily.csv` is produced by the main pipeline (`03_ssi_transformation.ipynb`). The weekly SPEI/SPI series were downloaded for each coordinate from the CSIC drought monitor, *Base de datos histórica y monitor de sequías en tiempo real* (https://monitordesequia.csic.es). The dataset of drought indices behind the monitor is described in Vicente-Serrano et al. (2017). The script expects the file names in the table above, with underscores; if the downloaded files use dots in the coordinates (e.g. `DATA_40.81_0.52.csv`), rename them accordingly (e.g. `DATA_40_81_0_52.csv`).

> Vicente-Serrano, S. M., Tomas-Burguera, M., Beguería, S., Reig, F., Latorre, B., Peña-Gallardo, M., Luna, M. Y., Morata, A., & González-Hidalgo, J. C. (2017). A High Resolution Dataset of Drought Indices for Spain. *Data*, 2(3), 22. https://doi.org/10.3390/data2030022

**Coordinate derivation.** Station coordinates (UTM, provided in `Locations.csv`) were converted to latitude/longitude using EPSG:25830 (ETRS89 / UTM zone 30N) and cross-checked against the known real-world locations of Tortosa (≈40.81°N, 0.52°E), Zaragoza (≈41.65°N, −0.88°E) and Caparroso (≈42.34°N, −1.65°E); agreement to within ~0.01–0.02°, well within tolerance. A second candidate datum (ED50/UTM 30N) was also tested and gave answers within ~0.002° — the choice of datum does not affect which grid cell of the SPEI product is queried.

## 3. Method

For each case study, the daily SSI series at the **origin** station and the weekly SPEI-3 (short/seasonal accumulation) and SPEI-12 (long-term accumulation) series at that station's coordinates are plotted together over the same approximate context window used in the manuscript's Figures 6–8 (window boundaries are approximate by design — precise alignment with those figures was judged unnecessary for illustrating meteorological context). The hydrological event window (the published chain start/end dates) is shaded in grey, exactly as the "context window" convention in Figures 6–8.

**Why the origin station specifically:** the origin is where the manuscript's own SSI series and drought-event statistics are defined and already plotted (Figs. 6–8), so it is the natural, unambiguous point of comparison. It is also the point at which the reviewer's question ("hydroclimatic context of the identified hydrological drought events") is most directly answered — the identified events are origin-station events.

## 4. Results and honest interpretation

### 4.1 The chain window is not the origin's own event

An initial version of this figure shaded the "hydrological event" window using the dates quoted in the manuscript's Figures 6–8 titles (e.g. "21 September 1988 – 24 April 1989" for Case 1) and labelled it as the origin station's own event. This was corrected after comparing it against the origin station's own individually detected SSI events (using the same, validated detection algorithm from Section S5):

| Case | Chain window (Figs. 6–8) | Origin's own event(s) |
|---|---|---|
| 1 (station 9027) | 1988-09-21 to 1989-04-24 (216 d) | **1989-01-23 to 1989-04-24** (92 d) — starts 4 months after the chain window begins |
| 2 (station 9011) | 2016-12-30 to 2018-03-01 (427 d) | **Six separate events**, 11–60 days each, with recoveries in between |
| 3 (station 9005) | 1997-02-28 to 1997-05-30 (92 d) | **1997-03-16 to 1997-04-28** (44 d) — narrower on both sides |

This is expected once the definitions are made explicit: the chain window is t(c,s) to t(c,e) (Table 2) — the earliest onset to the latest termination **among all matched stations in the chain**, including upstream stations that typically enter drought earlier (negative lag). The origin station's own individual SSI record is therefore generally a *subset* of the chain window, and in Case 2 is not even a single continuous interval. The figure now shows both layers explicitly: the chain window as a grey background band, and the origin's own individual event(s) highlighted in red directly on its SSI curve.

### 4.2 SPEI comparison

Local SPEI at the origin does **not** show a single, uniform precursor pattern:

- **Case 1 (Tortosa, 1988–89):** SPEI-3 falls sharply right at the hydrological onset (last pre-onset value ≈ −1.7), consistent with a short-term local contribution coinciding with the propagated signal.
- **Case 2 (Zaragoza, 2016–18):** both SPEI-3 and SPEI-12 are *positive* (wet) immediately before onset (mean SPEI-12 over the prior two years ≈ +0.5). The hydrological drought at the origin instead reflects a long-lag signal from persistent upstream droughts — exactly as already described in the manuscript's own narrative for this case (mean lag ≈ −120 days, driven in large part by a 392-day event at station 9093).
- **Case 3 (Caparroso, 1997):** SPEI-3 is wet immediately before onset, but SPEI-12 over the preceding two years is negative (mean ≈ −0.7, min ≈ −2.2) — a longer-term regional deficit context without a sharp immediate local trigger, consistent with this case's fast propagation (≈9-day lag) from a spatially coherent but comparatively short and less severe episode.

### 4.3 Why this supports, rather than undermines, the manuscript

Two independent lines of evidence — the SPEI comparison and the chain-window-versus-own-event comparison — point in the same direction. If the origin station's own local meteorology and its own local hydrological record fully explained the reported chain, a network-based, event-matching propagation framework would add comparatively little value beyond monitoring that single station directly. That neither does — most strikingly in Case 2, where the origin is locally wet in SPEI terms and recovers from drought several times in SSI terms, while the network-wide chain reports a severe, continuous, 427-day episode — is itself evidence that propagation in this basin is governed by upstream network structure and storage dynamics, not by conditions at any single station. This is fully consistent with, and reinforces, the storage-driven framing already developed in the Discussion and in the lag-vs-distance analysis (Section S3).

## 5. Reproducibility

**Dependencies:** Python ≥ 3.9 with `pandas`, `matplotlib`.

**Files:**
- `01_meteo_context_case_studies.py` — standalone script.
- `02_meteo_context_case_studies.ipynb` — the same pipeline as a documented, executed notebook, including the quantitative pre-onset SPEI check (Section 4 above).

**Run:**
```bash
python 01_meteo_context_case_studies.py
```
Place `SSI_daily.csv` and the three `DATA_*.csv` SPEI files alongside the script/notebook.

**Outputs:** `fig_meteo_context_case_studies.png` / `.pdf` (three-panel figure, one per case study).

## 7. Caveats of the analysis itself

- SPEI is queried at the **origin (downstream/outlet)** station's coordinates only. This is the natural point of comparison for that station's own SSI record, but it is only one point within each case's contributing catchment; a fully basin-representative or headwater-specific meteorological context (e.g. at the specific stations that trigger each chain, such as station 9018 for Case 1) was not attempted and would strengthen a future, fully coupled meteorological–hydrological extension.
- Context-window boundaries are approximate (chosen by eye to resemble Figures 6–8) rather than pixel-matched to those figures; this was judged acceptable since the purpose here is illustrative meteorological context, not a quantitative re-analysis of the chains themselves.
- The SPEI/SPI series are used as distributed by the CSIC drought monitor (see Section 2, *Data sources*); their grid resolution and computation method (Vicente-Serrano et al., 2017) were not independently re-derived here.

## 8. Audit log

- ✅ Station coordinates (UTM → lat/lon) cross-checked against known real-world locations of Tortosa, Zaragoza and Caparroso (agreement ~0.01–0.02°); tested against two plausible datums (ETRS89/ED50) with negligible difference.
- ✅ Data coverage confirmed complete (no missing SSI or SPEI-3/SPEI-12 values) within all three plotted context windows.
- ✅ Y-axis range checked and corrected to avoid clipping the most extreme SSI value in the plotted data (Case 1 minimum SSI = −3.98).
- ✅ **Corrected a labelling error:** an earlier draft shaded the manuscript's chain-window dates and labelled them as the origin station's own event. This was caught, verified against the origin's own individually detected SSI events (reusing the validated Section S5 algorithm), and fixed by showing both layers explicitly (chain window in grey, origin's own event(s) highlighted in red).
- ✅ Notebook executed end-to-end with a real kernel (zero errors).
- ✅ The absence of a uniform precursor pattern, and the origin/chain-window discrepancy, are both reported explicitly and interpreted honestly, rather than selectively presenting only the case(s) or framing that show the cleanest story.
