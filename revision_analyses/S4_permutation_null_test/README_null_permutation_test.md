# Surrogate / network-permutation null test — analysis README

**Manuscript:** *Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework* (Journal of Hydrology, major revision).
**Role of this analysis:** direct response to **Reviewer #3, comment 2** ("large-scale climatic anomalies... may cause the framework to incorrectly interpret coincident events as propagation chains"). Cross-referenced from the response to **Reviewer #4, comment 2** (Section S3) and **Reviewer #2, comment 5** (groundwater). This is the analysis that underpins the claim, made in those two responses, that "the detected chains substantially exceed chance co-occurrence."

**Headline result.** Of five summary statistics, three are significant at p < 0.05 under a network-label permutation null (connected fraction of origin events, p = 0.020; mean chain size, p = 0.012; total number of valid matched pairs, p = 0.004); the remaining two trend consistently in the same direction without reaching conventional significance (mean propagation fraction, 84th percentile; raw chain count, 88th percentile). **Important nuance:** the null distribution itself is not centred near zero — under preserved regional climatic synchrony but random network placement, a connected fraction of ~50 % is already expected by chance. The true network placement adds a modest but statistically detectable increment (~6 percentage points) on top of that climatic baseline. The chains are not statistically indistinguishable from chance co-occurrence, but a substantial part of the raw matching rate does reflect genuine, basin-wide climatic synchrony — consistent with the manuscript's own discussion of this phenomenon.

---

## 1. Rationale

The matching criterion (temporal overlap ≥ M, lag ≤ 0) verifies that an upstream event precedes or coincides with a downstream event and overlaps it sufficiently. It does **not**, by construction, rule out the possibility that both events are independent responses to the same regional climatic driver (e.g. a basin-wide dry year) rather than a genuine upstream-to-downstream transmission. Reviewer #3 raises exactly this concern. Testing it directly requires a null model in which regional climatic synchrony is preserved (so the test does not simply penalise the existence of synchrony, which is real and expected) while the correspondence between a station's drought history and its true network position is destroyed. If the real network placement still produces more/larger chains than random placement of the same histories, that is evidence of genuine network structure over and above climatic coincidence.

## 2. Inputs

| File | Content | Notes |
|---|---|---|
| `data/SSI_drought_events.csv` | Full event catalogue: `station_id, event_id, start_date, end_date, duration, severity, intensity, peak_SSI, peak_date, severity_hm3` | 2,625 events across 33 stations — matches Table 3 of the manuscript exactly |
| `data/upstream_connectivity.csv` | Per-station list of all hydraulically upstream stations | `;`-separated; `upstream_chain` is a comma-separated string of station IDs |

**Station-ID consistency (checked in code):** the 33 station IDs in the two files are identical sets. This is a necessary condition for the permutation (which reassigns whole catalogues between station-ID "slots") to be well defined.

## 3. Method, step by step

### 3.1 Reuse the published algorithm, not a reconstruction
The matching algorithm (asymmetric search window → real-overlap filter → four-level deterministic tie-break → lag) is taken directly from `05_ssi_propagation.ipynb` (Section 2.1.3 of the manuscript), supplied by the author. This avoids the risk of subtly reimplementing the method incorrectly.

### 3.2 Vectorise for speed, and validate the vectorisation is exact
The original implementation takes ~13 s per run (mostly `pandas.apply` calls with `Timestamp` arithmetic), which would make ~2,000 repetitions impractical (~4 hours). A numpy-vectorised re-implementation (`run_matching_vectorized` in `propagation_null_core.py`) reduces this to ~0.08 s (**~160–170× faster**).

**Before trusting the fast version for the permutation loop, it is validated against the original, unvectorised code on the real data: the two implementations return the exact same set of 2,075 (origin station, origin event, upstream station, upstream event, lag) tuples — not just similar summary statistics, but bit-for-bit identical matched pairs.** This check is included as an executable assertion in both the script and the notebook (it will halt execution if it ever fails, e.g. after a future change to the algorithm).

**A subtlety that was explicitly checked before vectorising:** the asymmetric search-window pre-filter is *not* logically redundant with the final overlap/lag filters. A constructed counter-example (an origin event long enough to overlap with a late-occurring, high-overlap upstream candidate that starts *after* the origin, alongside an earlier, lower-overlap true precursor) shows that dropping the window filter as a "shortcut" would silently change which candidate is tie-broken as "best" at a station, and hence the final result. The window filter is therefore preserved exactly in the vectorised code.

### 3.3 Harness self-check: identity mapping reproduces the published numbers
Before running any permutation, the harness is run once with **no shuffling** (the identity mapping) through the exact same code path (`build_station_arrays` → `run_matching_vectorized` → `summarise`) used for every permutation. This reproduces, without any adjustment: **12** origin stations, **928** origin events, **525** chains, **56.6 %** connected — the published Section 2.2 numbers — and, further, mean chain size **4.952** and mean propagation fraction **0.418**, matching the manuscript's Figure 5 statistics to three decimal places. This is strong independent confirmation that the summary-statistic definitions used here (`chain_size_Nc = N_up + 1`, `propagation_fraction = N_up / N_avail`) are exactly the ones underlying the published figures.

### 3.4 Null model: network-label permutation
For each of `N_PERM = 1,999` permutations (seed = 20260713, `numpy.random.default_rng`):

1. A random bijection of the 33 station IDs is drawn.
2. Each network "slot" (a station ID, used only for connectivity lookups) is assigned the event catalogue that originally belonged to a different, randomly chosen station ID. The **network topology** (`upstream_connectivity`, i.e. who is upstream of whom) is **never altered**. The **event catalogues** (dates, durations, severities) are **never altered** — only their assignment to network slots is randomised.
3. The validated matching algorithm is re-run on this relabelled configuration, and the five summary statistics are recorded.

This preserves regional climatic synchrony (the same 33 real drought histories, with their real dates, are always present somewhere in the network) while destroying the correspondence between a specific history and its true topological position.

### 3.5 Test statistics
Three primary, denominator-normalised statistics (comparable across permutations despite the varying number of origin events landing in the 12 origin-connected slots each time):

- **Connected fraction** (%) — mirrors the "56.6 % of origin events... exhibit at least one upstream connection" statement (Section 2.2 / Supplementary S1).
- **Mean chain size** (N_c = N_up + 1) — mirrors Figure 5b.
- **Mean propagation fraction** (f_p = N_up / N_avail) — mirrors Figure 5c.

Two secondary, non-normalised statistics are reported for completeness: raw chain count and total number of valid matched pairs (edges).

### 3.6 Significance test
One-sided exact permutation p-value, since the alternative hypothesis is directional (the true configuration produces *more* structure than random): 

```
p = (1 + #{null permutations with statistic >= observed}) / (1 + N_PERM)
```

With `N_PERM = 1,999`, the smallest attainable p-value is 1/2000 = 0.0005.

## 4. Results

| Statistic | Observed | Null mean ± SD | p-value | Percentile |
|---|---|---|---|---|
| Connected fraction (%) | 56.6 | 50.3 ± 3.1 | **0.020** | 98.0 |
| Mean chain size (N_c) | 4.95 | 4.28 ± 0.31 | **0.012** | 98.9 |
| Mean propagation fraction (f_p) | 0.418 | 0.393 ± 0.025 | 0.156 | 84.4 |
| Chain count (raw) | 525 | 480 ± 38 | 0.123 | 87.7 |
| Valid pairs (raw) | 2,075 | 1,578 ± 204 | **0.004** | 99.6 |

Full permutation-level results: `null_distribution.csv` (1,999 rows). Summary table: `permutation_test_summary.csv`.

## 5. Interpretation

Three of five statistics reach conventional significance (p < 0.05); the remaining two are not significant at that threshold but sit consistently in the upper 84th–88th percentile of the null distribution — directionally consistent, and most plausibly diluted by the fact that raw counts (unlike the normalised statistics) are also driven by how many total origin events happen to land in the 12 origin-connected network slots in a given permutation, which varies considerably (776–1,125 across permutations).

**The honest and important nuance is in the null distribution itself, not just the p-values.** Under a null that preserves real regional climatic synchrony but randomises network placement, a connected fraction of ~50 % already arises by chance — confirming that basin-wide, synchronous drought episodes (which the manuscript itself discusses) genuinely inflate the raw matching rate. The true network placement adds a further, smaller increment (~6 percentage points for connected fraction; a comparable relative increment for chain size) on top of that climatic baseline, and this increment is what the permutation test shows to be statistically real. This is a moderate, not an overwhelming, result — and it should be reported as such rather than oversold.

**Conclusion for the manuscript.** The reconstructed propagation chains are not statistically indistinguishable from what regional climatic synchrony alone would produce; the true upstream–downstream network structure contributes real, additional, detectable connectivity beyond that baseline. This supports the physical validity of the chains as reflecting genuine network transmission and not pure coincidence, while also validating part of the reviewer's concern — a non-trivial share of matched pairs likely does reflect co-occurring but not causally connected drought — a point already partly addressed by the Section S3 (lag vs. distance) analysis and worth stating explicitly as a shared limitation.

## 6. Reproducibility

**Dependencies:** Python ≥ 3.9 with `pandas`, `numpy`, `matplotlib`, `tqdm` (only needed if re-running the original notebook cell for comparison).

**Files:**
- `../propagation_null_core.py` — core module (algorithm, validated vectorisation, permutation harness), shared with S7 and S8 and kept once in `revision_analyses/`. Import this in any downstream script or notebook.
- `01_permutation_null_test.py` — standalone script: loads data, validates the vectorised algorithm against the original, runs the harness self-check, runs the permutation test, computes p-values, saves the figure.
- `02_permutation_null_test.ipynb` — the same pipeline as a documented, executed notebook (recommended entry point for inspecting intermediate outputs).

**Run:**
```bash
python 01_permutation_null_test.py
# or
jupyter nbconvert --to notebook --execute --inplace 02_permutation_null_test.ipynb
```

Place `SSI_drought_events.csv` and `upstream_connectivity.csv` in a `data/` subfolder alongside the script/notebook. With the fixed seed (20260713), results are exactly reproducible.

**Outputs:** `null_distribution.csv` (permutation-level results), `permutation_test_summary.csv` (the table in Section 4), `fig_null_permutation_test.png` / `.pdf` (the supplementary figure).

## 7. Caveats of the analysis itself

- The permutation reassigns **whole catalogues** between stations; it does not test whether the *specific* pairing of individual events within a catalogue matters (a finer within-catalogue permutation is possible but was not needed to answer the reviewer's question, which is about station-level correspondence).
- Raw (non-normalised) statistics are more variable across permutations because the composition of the 12 origin-connected slots changes every permutation; they are reported as secondary, corroborating evidence rather than primary statistics.
- As with any permutation test, statistical significance indicates the observed configuration is unlikely under the null model tested; it does not by itself quantify *what fraction* of any individual chain is "real" propagation versus coincidence — that question is better addressed qualitatively, per case, alongside the lag-vs-distance analysis (Section S3) and the discussion of co-occurring, disconnected droughts already present in the manuscript.

## 8. Audit log

- ✅ Station-ID sets confirmed identical between the event catalogue and the connectivity table.
- ✅ **Vectorised algorithm validated bit-for-bit** against the original, unvectorised notebook implementation (set-equality over 2,075 matched-pair tuples), after explicitly checking that the search-window filter is not a redundant shortcut (constructed counter-example).
- ✅ **Harness self-check passed:** identity mapping reproduces the published 928 / 525 / 56.6 % and Figure 5's mean chain size (4.952) and mean propagation fraction (0.418) to three decimals, via the same code path used for every permutation.
- ✅ Notebook executed end-to-end with a real kernel (zero errors); results identical to the standalone script under the same seed.
- ✅ Five test statistics reported transparently, including two that do not reach conventional significance — no selective reporting.
- ✅ Null-distribution baseline (~50 % expected by chance under preserved synchrony) reported explicitly rather than omitted, to avoid overstating the result.
