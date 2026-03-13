# Longitudinal Regimes of Arts and Cultural Engagement and Frailty Among Older Adults in the United States: A G-Formula Approach

This repository contains the analysis code accompanying the manuscript of the same title (currently under peer review).

## Citation

Manuscript currently under peer review. Citation details will be added upon publication.

## Data

This analysis uses data from the [Health and Retirement Study (HRS)](https://hrs.isr.umich.edu/), which is publicly available upon registration.

**Important:** The input Stata files (`.dta`) required by Script 1 are **not included** in this repository. These files were derived from raw HRS data through a separate pre-processing pipeline that is not part of this repository. To reproduce the analysis, place the following files in the `data/` folder:

| File | Description |
|------|-------------|
| `cams_covars00-18.dta` | Longitudinal covariates from HRS core interviews |
| `frailty_full.dta` | Frailty index across waves |
| `cams05_17.dta` | CAMS arts/cultural engagement items (2005–2017) |
| `trk2022tr_r.dta` | 2022 HRS Tracker file (eligibility and vital status) |

The first three files were prepared from raw HRS data through a separate pre-processing step. For questions about the data preparation, please contact Dr Feifei Bu (f.bu@ucl.ac.uk). For questions about the pipeline provided here, please contact the author of this repository (Martin Danka, martin.danka.21@ucl.ac.uk). The tracker file (`trk2022tr_r.dta`) is available directly from the [HRS website](https://hrs.isr.umich.edu/).

## Requirements

- **R** v4.4.1 (other versions may work but have not been tested)
- **Operating system:** macOS, Linux, or Windows (see [Parallelisation](#parallelisation) for OS-specific notes)
- **RAM:** ~32 GB recommended for the full imputation pipeline (50 imputations x 40 iterations)
- **Computation time:** The imputation and g-formula steps may require several hours of computation

## Setup

1. Clone or download the repository.
2. Place the required `.dta` files in the `data/` folder.
3. Set your working directory to the project root.
4. Restore package dependencies:

```r
install.packages("renv")  # if not already installed
renv::restore()
```

Alternatively, install packages manually. Each script loads its dependencies at the top.

## Repository Structure

| File | Description |
|---|---|
| `code/1-data-cleaning.R` | Imports pre-processed and raw Stata files, applies exclusions, derives variables, reshapes to wide format |
| `code/2-imputations.R` | Multiple imputation for missingness |
| `code/3-gformulaMI.R` | Effect estimation via GformulaMI |
| `code/4-gformulaICE.R` | Sensitivity analysis using gformula via iterative conditional expectation |
| `code/5-diagnostics.R` | Outcome model diagnostics (residual and calibration plots) |
| `code/helpers.R` | Shared utility functions sourced by all scripts |
| `renv.lock` | Snapshot of exact package versions |

## Running the Code

Run scripts **sequentially** (1 through 6). Each script saves intermediate `.rds` files to `data/`, which subsequent scripts load. On re-runs, you can comment out the heavy computation blocks and load from saved results instead.

Shared functions are stored in `helpers.R` and imported by each script via `source("code/helpers.R")`.

## Parallelisation

Scripts 2 and 4 use parallel computing. By default, they are configured for **macOS** using `future::plan(multicore)` with `ncore <- 5`. **On Windows**, or if `multicore` is unavailable, use `plan(multisession)` instead. Each script contains a commented-out "General set up for other machines" block that you can swap in. This uses `parallel::detectCores() - 2` to adapt to your hardware.

Results may differ slightly across parallelisation backends due to how random seeds interact with parallel processing.
