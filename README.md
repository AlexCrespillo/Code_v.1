# Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework

This repository contains the complete analytical pipeline for the paper:

> **Tracking Hydrological Drought Propagation in River Networks: An Event-Based Framework**

The framework identifies hydrological drought events at individual gauging stations using the Standardised Streamflow Index (SSI) and reconstructs propagation chains by linking temporally consistent upstream-downstream events based on river network connectivity.

## Study area

The method is applied to daily streamflow data from **33 gauging stations** in the **Ebro River Basin** (Spain) over the period **1961-2020**.

## Repository structure

```
.
├── data/                        # Input data (not tracked by Git — see below)
│   ├── spatial/                 # Shapefiles, DEM raster
│   ├── caudales_raw.csv         # Raw CEDEX discharge records
│   ├── upstream_connectivity.csv
│   └── ...
├── output/                      # Reproducible outputs (not tracked by Git)
│   └── figures/
├── 01_pre_imputation.ipynb      # Step 1: Quality control & station filtering
├── 02_imputation.ipynb          # Step 2: Hierarchical gap-filling
├── 03_ssi_transformation.ipynb  # Step 3: Daily SSI computation
├── 04_ssi_drought_events.ipynb  # Step 4: Drought event detection
├── 05_ssi_propagation.ipynb     # Step 5: Propagation chain construction
├── 06_ssi_propagation_sensitivity.ipynb  # Step 6: Parameter sensitivity analysis
├── 07_combined_chain_figure.ipynb        # Figure: Combined chain visualisation
├── 08_discharge_distribution_figure.ipynb # Figure: Discharge distribution
├── 09_drought_event_figure.ipynb          # Figure: Single drought event
├── 10_study_area_figure.ipynb             # Figure: Study area map + distribution
├── 11_chain_metrics_distribution_figure.py # Figure: Chain metric distributions
├── requirements.txt
└── README.md
```

## Workflow

The notebooks are numbered sequentially and should be executed in order. Notebooks 1-6 implement the analytical pipeline; notebooks 7-11 generate publication figures.

```
Raw discharge data
  │
  ▼
01  Quality control & station filtering
  │   → Filters 302 stations down to 45 (then 33 after SSI screening)
  │   → Thresholds: ≤30% NaN, ≤365-day gaps, ≤120-day zero-flow runs
  ▼
02  Hierarchical gap-filling (imputation)
  │   → Level 1: Log-interpolation (gaps ≤ 15 days)
  │   → Level 2: Donor-based log-regression (10 nearest stations, r > 0.75)
  │   → Level 3: DOY climatology fallback (station-specific median)
  ▼
03  Daily SSI computation
  │   → 6 candidate distributions (Gamma, LogNormal, Weibull, PearsonIII, LogLogistic, GEV)
  │   → DOY ±15-day pooling window (~1,860 obs per fit)
  │   → Selection by Shapiro-Wilks W statistic
  │   → Zero-flow probability mass adjustment
  ▼
04  Drought event detection
  │   → SSI < -1.28 threshold (P10)
  │   → Pooling: merge events separated by ≤ 10 days
  │   → Minimum duration: 5 days
  │   → Volumetric severity (hm³) via P10 discharge threshold
  ▼
05  Propagation chain construction
  │   → Asymmetric search window (W = 45 days before origin onset)
  │   → Minimum temporal overlap (M = 5 days)
  │   → Deterministic tie-breaking (overlap → severity → proximity → ID)
  │   → Lag filter: retain only upstream_start ≤ origin_start
  ▼
06  Parameter sensitivity analysis
      → Grid search: W ∈ [0, 120] days, M ∈ [1, 30] days (750 combinations)
      → Metrics: RSI, marginal gain, Kneedle elbow detection
```

## Data

Input data files are **not included** in this repository due to size constraints. To reproduce the analysis:

1. Place raw discharge data (`caudales_raw.csv`) and station metadata files in the `data/` folder.
2. Place spatial data (DEM, shapefiles) in `data/spatial/`.
3. Run the notebooks in order (01 → 06), which will generate all intermediate CSV files in `data/`.

The raw discharge data was obtained from the [CEDEX](https://ceh.cedex.es/anuarioaforos/default.asp) hydrological monitoring network.

## Installation

```bash
# Clone the repository
git clone https://github.com/<your-username>/drought-propagation-framework.git
cd drought-propagation-framework

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

### Geospatial dependencies

Notebooks 07 and 10 require `cartopy`, `geopandas`, and `rasterio`, which depend on system-level libraries (GEOS, PROJ, GDAL). On some systems you may need to install these via conda:

```bash
conda install -c conda-forge cartopy geopandas rasterio
```

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| SSI threshold | -1.28 | P10 of standard normal (moderate drought onset) |
| Pooling gap | 10 days | Max inter-event gap for merging |
| Min duration | 5 days | Minimum event span after pooling |
| Search window (W) | 45 days | Look-back window for upstream precursors |
| Min overlap (M) | 5 days | Minimum temporal overlap for event matching |

## Citation

If you use this framework, please cite:

> [Paper citation will be added upon publication]

## License

[License to be specified]
