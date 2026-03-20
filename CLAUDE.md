# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R research project benchmarking 7 multivariate outlier detection algorithms for an academic article. The novel methods are DSO, DSO-SDC, and DSO-SDM, compared against baselines (REPPlab, HDO, PP-MCD, ICS) on simulated and real datasets.

## Running the Code

All scripts are run via `source()` in R or RStudio. There is no build system or CLI entrypoint.

```r
# One-time: generate simulation data (produces data/raw/simulations_data.rds)
source("R/Generate_Simulated_Data.R")

# Main benchmarks (independent, run as needed)
source("notebooks/Simulaciones_Art.R")  # Simulated data (~50 replications × 150 scenarios)
source("notebooks/Real_Data.R")         # Real datasets (Wine, Satimage-2, Breastw)
```

Results are written to `outputs/run_YYYYMMDD_HHMMSS/` — CSV and RDS files. The `outputs/` directory is gitignored.

## Environment

Uses `renv` for reproducibility (R 4.4.2). The `.Rprofile` auto-activates renv. To restore the environment:

```r
renv::restore()
```

## Architecture

### Source vs. Notebooks

- **`R/`** — Reusable source files, sourced by notebooks
- **`notebooks/`** — Analysis entry points; source `R/utils_analysis.R` at the top (and `R/OL_skew_self_V2.R` where needed)

### Key Files

| File | Role |
|------|------|
| `R/OL_skew_self_V2.R` | Core algorithm primitives: skewness projection (`max_skew`), adjusted outlyingness (`adj.outly`), data generators (`GenAtip`, `GenAtip_FICM`) |
| `R/utils_analysis.R` | Shared utilities: `cluster_to_labels`, `calculate_metrics`, `compute_confusion_metrics`, `save_consolidated_results`, `create_run_directory` |
| `notebooks/Simulaciones_Art.R` | Benchmark over simulated data; defines `pena_prieto_transformation`, `SDC`, `SDM`, `calculate_outlyingness`, `calculate_r_i`, `quantil_id_robust`, `calculate_ab_robust`, `robust_outlier_detection` inline |
| `notebooks/Real_Data.R` | Same 7 methods applied to real `.mat` datasets; defines the same helpers inline plus `component_selection`, `huber` |
| `data/raw/simulations_data.rds` | Pre-generated simulation cache (gitignored) |

### Inline helpers (NOT in shared files)

Both notebooks define their own local copies of `pena_prieto_transformation`, `SDC`, `SDM`, `calculate_outlyingness`, `calculate_r_i`, `convert_to_matrix`, and `cm_to_row`. These are **not** imported from `R/utils_analysis.R`. The commented-out duplicate versions inside the notebooks are superseded dead code and should not be restored.

### Outlier Detection Pipeline (per method)

1. Transform/project the data (method-specific)
2. Compute an outlyingness score per observation
3. Convert scores to binary labels via `cluster_to_labels()` (1D k-means; cluster with higher mean = outlier cluster 2)
4. Evaluate with `compute_confusion_metrics()` / `calculate_metrics()`

### Novel Methods (DSO family) — pipeline differences between notebooks

All three share the first two steps regardless of notebook:
1. **DOBIN** rotation (from the `dobin` package)
   - `Simulaciones_Art.R`: uses `dobin(X)$rotation`
   - `Real_Data.R`: uses `dob$basis` with fallback to `dob$rotation`
2. **Peña-Prieto transformation** (`pena_prieto_transformation`) — orthogonal complement projection; takes `X` and the DOBIN rotation matrix (not the full `dobin` object)

After step 2, the notebooks diverge:

| Step | `Simulaciones_Art.R` | `Real_Data.R` |
|------|----------------------|---------------|
| Component selection | **Skipped** — all columns of T3 used | `component_selection()` keeps components with kurtosis+skewness above median |
| Before SDC/SDM | Raw `selected_data` passed directly | `huber(selected_data)` applied first (clips MAD-standardized values at 97.5th percentile columnwise) |
| DSO-SDC | `SDC(selected_data, r_i)` | `SDC(selected_data_h, r_i)` where `r_i` computed on huberized data |
| DSO-SDM | `SDM(selected_data)` | `SDM(selected_data_h)` |

**DSO**: simple componentwise outlyingness via MAD (`calculate_outlyingness`)
**DSO-SDC**: adaptive weighting based on scaled direction calibration (`SDC()`, takes `x` and precomputed `r_i`)
**DSO-SDM**: adaptive weighting with scaling modulation (`SDM()`, computes its own internal `r_i_SDM`, distinct from the `r_i` used in SDC)

### HDO implementation (both notebooks)

Both notebooks now use the same custom random-projection method:
- `calculate_ab_robust(d, n, alpha=0.05)` — simulates 100 random projections of a standard normal dataset to compute thresholds `a0` and `b0`
- `robust_outlier_detection(X, a0, b0, num_projections=50)` — flags observation `i` as outlier if any of its 50 random projections exceeds `b0`

`depthTools` and `doRNG` are **not used** and must not be added back to the package lists.

### Simulation Design

`simulations_data.rds` contains a nested list organized by:
- Contamination model: FDCM (casewise) and FICM (cellwise)
- Simulation mode: usual / type A / type B (modes 1–3)
- Dimensions: p ∈ {10, 30, 50}
- Contamination levels: ε ∈ {0.05, 0.10, 0.20, 0.30, 0.40}
- Outlier distance: d ∈ {3, 4}
- 50 replications each

Key lookup pattern: `sim_data[["FDCM_mode1_p10_alpha0.05_delta3_c0.5"]][[rep_index]]`

`Real_Data.R` iterates `simulations_to_run <- 2:4`, corresponding to Wine (2), Satimage-2 (3), Breastw (4).

### Parallelism

`Simulaciones_Art.R` registers `doParallel` with `detectCores() - 1` workers, but the inner replication loop currently uses a regular `for` loop (not `foreach`). The outer scenario loop is also sequential.

## Important Conventions

- **Seeds**: `set.seed(50)` for simulation generation; `set.seed(123)` in notebooks
- **Thresholds**: REPPlab uses 95th percentile; PP-MCD uses χ²(0.95) Mahalanobis cutoff; ICS uses `ics.outlier(level.test = 0.3)`; HDO uses simulation-based `b0` threshold (see HDO section)
- **Error handling**: PP-MCD calls are wrapped in `tryCatch`; failed replications return `NA` metrics
- **Output format**: `confusion_ALL.csv` and `metrics_ALL.csv` are the consolidated result files consumed for the article; intermediate per-scenario CSVs are also written during the run
