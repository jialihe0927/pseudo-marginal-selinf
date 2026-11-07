# pseudo-marginal-selinf

Code for reproducing the experiments in *Pseudo-marginal Sampling for Selective Inference*.

## Installation

Install R and the following packages:

```sh
Rscript -e "install.packages(c('posterior', 'glmnet', 'ggplot2', 'patchwork'))"
```

`posterior` is used for diagnostics, `glmnet` for CV-Lasso preparation and
sampling, and `ggplot2` and `patchwork` for figures. PDF output requires R with
Cairo support. File Drawer calculations and Changepoint sampling use base R.
Tested with R 4.4.3, posterior 1.7.0, glmnet 4.1-9, ggplot2 4.0.3 and
patchwork 1.3.0; these versions are not enforced by the code.

Install Git and Git LFS, then obtain the code and data:

```sh
git clone https://github.com/jialihe927-maker/pseudo-marginal-selinf.git
cd pseudo-marginal-selinf
git lfs install --local
git lfs pull
```

The two large CV-Lasso RDS files, `calibration_lasso.rds` and `power_lasso.rds`,
use Git LFS. Download their actual data and keep them at the supplied paths;
GitHub's Download ZIP may omit the LFS data.

## Repository structure

| Path | Purpose |
| --- | --- |
| [reproduce.R](reproduce.R) | Unified entry for paper tables, diagnostics and figures |
| [experiments/file_drawer/](experiments/file_drawer/) | File Drawer sampling complexity |
| [experiments/cv_lasso/](experiments/cv_lasso/) | CV-Lasso experiments |
| [experiments/changepoint/](experiments/changepoint/) | Changepoint experiments |

Both applications contain four experiment directories:

| Directory | Purpose | CV-Lasso | Changepoint |
| --- | --- | --- | --- |
| `efficiency_comparison` | Sampling efficiency | [code](experiments/cv_lasso/efficiency_comparison/) | [code](experiments/changepoint/efficiency_comparison/) |
| `global_null_calibration` | Global-null calibration | [code](experiments/cv_lasso/global_null_calibration/) | [code](experiments/changepoint/global_null_calibration/) |
| `power_comparison` | Power comparison | [code](experiments/cv_lasso/power_comparison/) | [code](experiments/changepoint/power_comparison/) |
| `sampling_validation` | Sampling validation | [code](experiments/cv_lasso/sampling_validation/) | [code](experiments/changepoint/sampling_validation/) |

## Reproduce figures and tables

To reproduce figures and tables from the supplied CV-Lasso and Changepoint
results, run from the repository root:

```sh
Rscript reproduce.R ../paper-output
```

This summarizes the supplied results, recomputes coordinate and inference
indicator diagnostics from saved samples, and generates their figures.
**File Drawer is run from scratch** before plotting; the command does not rerun
the CV-Lasso or Changepoint experiments.

Tables and diagnostics are printed to the console. All nine PDFs are saved in
the specified directory: File Drawer, two calibration plots, two power plots,
two sampling-validation plots, and one combined inference figure per application.
The default output is `../paper-output`, beside the repository. Relative output
paths use the working directory; the root entry requires an output directory
outside the repository. It reads the supplied RDS paths, not the new files in
the examples below.

## Run experiments

Start each code block below from the repository root; its final `cd` returns
there. New RDS and PDF files are written in the corresponding experiment
directory with a `new_` prefix. Choose unused names to keep the supplied results
and figures intact. Plotting a new RDS writes a PDF with the same basename;
efficiency table commands print to the console.

Complete sampling can take substantial time. Seeds and stream assignments are
set in the code or supplied configurations; reruns can differ numerically across
R, package and numerical-library environments. Experiment settings are defined
in the preparation functions and the configurations loaded by the runners.

### File Drawer

Run the complete experiment and save its PDF to a new path:

```sh
cd experiments/file_drawer
Rscript run_fd.R new_sampling_complexity_fd.pdf
cd ../..
```

This runs 21 simulated configurations and computes 14 analytic rejection results.
Results stay in memory; no result RDS is saved. Without an output argument,
`run_fd.R` replaces `sampling_complexity_fd.pdf` in this directory.

### CV-Lasso

In inference figures, Proposed is the model-only method (`model_only`).

#### Sampling efficiency

Run all three methods on the supplied fixed benchmark in `efficiency_lasso.rds`,
then print the new table. This does not regenerate the benchmark; rejection uses
3,000,000 proposals.

```sh
cd experiments/cv_lasso/efficiency_comparison
Rscript sampling_efficiency_lasso.R all new_efficiency_lasso.rds
Rscript table_efficiency_lasso.R new_efficiency_lasso.rds
cd ../../..
```

#### Global-null calibration

Prepare a new bank, run inference on that bank, then plot its results:

```sh
cd experiments/cv_lasso/global_null_calibration
Rscript prepare_calibration_lasso.R new_bank_calibration_lasso.rds
Rscript calibration_lasso.R new_calibration_lasso.rds new_bank_calibration_lasso.rds
Rscript plot_calibration_lasso.R new_calibration_lasso.rds
cd ../../..
```

#### Power comparison

```sh
cd experiments/cv_lasso/power_comparison
Rscript prepare_power_lasso.R new_bank_power_lasso.rds
Rscript power_lasso.R new_power_lasso.rds new_bank_power_lasso.rds
Rscript plot_power_lasso.R new_power_lasso.rds
cd ../../..
```

#### Sampling validation

Regenerate samples using the supplied validation configuration, analyze them,
and plot the new samples:

```sh
cd experiments/cv_lasso/sampling_validation
Rscript sampling_validation_lasso.R --generate new_validation_lasso.rds
Rscript plot_sampling_validation_lasso.R new_validation_lasso.rds
cd ../../..
```

### Changepoint

In inference figures, Proposed is the full adaptive-selector method (`full_selector`).

#### Sampling efficiency

Run all three methods on the supplied fixed benchmark in `efficiency_cp.rds`,
then print the new table. This does not regenerate the benchmark; rejection uses
3,000,000 proposals.

```sh
cd experiments/changepoint/efficiency_comparison
Rscript sampling_efficiency_cp.R all new_efficiency_cp.rds
Rscript table_efficiency_cp.R new_efficiency_cp.rds
cd ../../..
```

#### Global-null calibration

```sh
cd experiments/changepoint/global_null_calibration
Rscript prepare_calibration_cp.R new_bank_calibration_cp.rds
Rscript calibration_cp.R new_calibration_cp.rds new_bank_calibration_cp.rds
Rscript plot_calibration_cp.R new_calibration_cp.rds
cd ../../..
```

#### Power comparison

```sh
cd experiments/changepoint/power_comparison
Rscript prepare_power_cp.R new_bank_power_cp.rds
Rscript power_cp.R new_power_cp.rds new_bank_power_cp.rds
Rscript plot_power_cp.R new_power_cp.rds
cd ../../..
```

#### Sampling validation

Regenerate samples using the supplied validation configuration, analyze them,
and plot the new samples.

```sh
cd experiments/changepoint/sampling_validation
Rscript sampling_validation_cp.R --generate new_validation_cp.rds
Rscript plot_sampling_validation_cp.R new_validation_cp.rds
cd ../../..
```
