# Revision analyses

Analyses developed during peer review, each corresponding to a section of the
Supplementary Material. The six first-round folders (`S3`–`S8`) contain a
standalone script, a documented notebook, a README describing method and
caveats, and the resulting data files. The two second-round folders
(`case1_imputation_sensitivity/`, `excluded_stations/`) contain the standalone
script and its results; their method and outcome are described in the script
docstring.

| Folder | Supplementary section | Purpose |
|---|---|---|
| `S3_lag_vs_distance/` | S3 | Tests whether propagation lag relates to along-channel hydraulic distance |
| `S4_permutation_null_test/` | S4, S1.4 | Network-permutation null model; temporal-overlap vs gap sensitivity |
| `S5_detection_sensitivity/` | S5 | Sensitivity of the event catalogue to pooling gap and minimum duration |
| `S6_meteorological_context/` | S6 | SPEI context and origin-station records for the three case studies |
| `S7_non_propagating_events/` | S7 | Isolated downstream events and non-propagating upstream events |
| `S8_national_catalogue/` | S8 | External corroboration against the national drought catalogue |
| `case1_imputation_sensitivity/` | S6 | Robustness of Case 1 to imputed days at the origin station |
| `excluded_stations/` | S9 | Metadata and reservoir-proximity test for excluded stations |

## Shared modules

Two modules are used by several analyses and are kept once, in this folder:

- `propagation_null_core.py` reproduces the published propagation-matching
  algorithm in vectorised form (used by `S4`, `S7` and `S8`). It is validated
  pair-for-pair against the original implementation in
  `S4_permutation_null_test/01_permutation_null_test.py`, and is the same
  matching engine as `05b_ssi_propagation_vectorised.ipynb` in the main
  pipeline.
- `detection_sensitivity_core.py` reproduces the event-detection pipeline
  (used by `S5`), validated against the published Table 3.

Scripts and notebooks add this folder to the Python path before importing
them, so no copies are needed in the subfolders. `ssi_core_extracted.py` is
used only by `case1_imputation_sensitivity/` and stays there.

## Running an analysis

Run each script or notebook from its own subfolder, for example:

```bash
cd revision_analyses/S4_permutation_null_test
python 01_permutation_null_test.py
```

Outputs are written to that same subfolder.

## Data

The inputs are the files of the main pipeline (`SSI_daily.csv`,
`SSI_drought_events.csv`, `upstream_connectivity.csv`, etc.) plus a few
external datasets, and they are not included in the repository (see the main
README). `S4`, `S5`, `S7` and `S8` read them from a `data/` subfolder inside
the analysis folder; `S3`, `S6`, `case1_imputation_sensitivity` and
`excluded_stations` read them from the analysis folder itself. Each README
(or script docstring) lists the exact files required.
