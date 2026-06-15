# Trans-GAM Simulation Code

This repository contains the simulation and model-fitting code for the source-score residual transfer GAM estimator used in the manuscript. The release is intended for reviewer inspection and reproducibility of the simulation workflows. Plotting scripts, manuscript tables, local result files, and protected real-data processing code are not included.

## Contents

- `R/transfer_gam_functions.R`: model-fitting functions, IPCW utilities, transfer estimators, calibration/selection variants, and baseline learners.
- `R/sim_setup_weibull.R`: Weibull simulation data-generating mechanisms and evaluation helpers.
- `R/sim_setup_aftlogistic.R`: AFT-logistic simulation data-generating mechanisms and evaluation helpers.
- `R/server_common.R`: shared command-line parsing, method lists, output handling, and summary helpers.
- `scripts/main_simulation_server.R`: command-line entry point for running simulations.
- `scripts/submit_simulation_array.sh`: example Slurm array job.
- `env/environment.yml`: conda environment specification.

## Installation

Create the R environment with conda:

```bash
conda env create -f env/environment.yml
conda activate transgam-r
```

If the environment already exists:

```bash
conda env update -f env/environment.yml --prune
conda activate transgam-r
```

## Quick Smoke Test

From the repository root:

```bash
Rscript scripts/main_simulation_server.R \
  --model weibull \
  --methods stat,transfer \
  --settings aligned_linear \
  --n-rep 1 \
  --n-s 200 \
  --n-t 50 \
  --n-test 100 \
  --tq "{0.50}" \
  --out-root results/smoke \
  --cores 1 \
  --no-parallel
```

The run writes:

- `all_results_<model>_nonlinear.csv`
- `summary_table_<model>_nonlinear.csv`
- `eta_alignment_<model>.csv`
- run metadata files recording arguments, methods, and prediction horizons.

## Main Simulation Examples

Weibull simulation:

```bash
Rscript scripts/main_simulation_server.R \
  --model weibull \
  --methods stat,transfer,gam_gam_basis,oracle \
  --settings all \
  --n-rep 100 \
  --n-s 2000 \
  --n-t 300 \
  --n-test 1000 \
  --tq "{0.25,0.50,0.75}" \
  --out-root results/weibull_300 \
  --cores 8
```

AFT-logistic simulation:

```bash
Rscript scripts/main_simulation_server.R \
  --model aftlogistic \
  --methods stat,transfer,gam_gam_basis,oracle \
  --settings all \
  --n-rep 100 \
  --n-s 2000 \
  --n-t 300 \
  --n-test 1000 \
  --tq "{0.25,0.50,0.75}" \
  --out-root results/aftlogistic_300 \
  --cores 8
```

Available setting keys:

- `aligned_linear`
- `aligned_nonlinear`
- `misaligned_linear`
- `misaligned_nonlinear`
- `all`

Common method groups:

- `stat`: logistic, GAM, and Cox source/target/pooled baselines.
- `ml`: random survival forest and IPCW XGBoost source/target/pooled baselines.
- `rsf`: RSF baselines only.
- `xgboost`: IPCW XGBoost baselines only.
- `transfer`: fixed transfer estimators.
- `gam_gam_basis`: GAM-GAM source/target basis ablations, including the main GG-SrcTP estimator.
- `auto`: automatic selection variants.
- `oracle`: simulation-only oracle comparator.
- `all`: all implemented simulation methods.

## Slurm

An example array job is provided in `scripts/submit_simulation_array.sh`. Submit with defaults:

```bash
sbatch scripts/submit_simulation_array.sh
```

Override settings at submission:

```bash
sbatch --export=ALL,MODEL=aftlogistic,N_T=150,METHODS=stat,transfer,gam_gam_basis,oracle scripts/submit_simulation_array.sh
```

The array script uses `--rep-start` and `--rep-end` to split Monte Carlo replications into chunks.

Merge chunk outputs after the array finishes:

```bash
Rscript scripts/merge_simulation_chunks.R \
  --chunk-root results/weibull_300/chunks \
  --out-dir results/weibull_300/combined
```

## Optional TransCox Comparator

The code contains a wrapper for the TransCox comparator, but the external TransCox source code is not bundled here. To run `--methods transcox`, place `TransCox-master/` next to this repository or set:

```r
options(translogistic.transcox_dir = "/path/to/TransCox-master")
```

The proposed Trans-GAM estimators and the simulation data-generating mechanisms do not require TransCox.

## Notes

- The simulation code uses only synthetic data generated at run time.
- No local absolute paths are required.
- Output directories are created under the user-specified `--out-root` or `--out-dir`.
- The real-data analysis scripts are omitted because individual-level clinical data cannot be redistributed.
