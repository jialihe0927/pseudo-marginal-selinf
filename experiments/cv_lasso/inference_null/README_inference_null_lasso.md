# CV-Lasso global-null calibration

This experiment evaluates selective p-value calibration after randomized,
fixed-fold CV-Lasso selection under the global null.

## Files

- `functions_null_lasso.R`: data generation, selectors, selection events,
  conditional reconstruction, data-bank construction, and validation.
- `samplers_null_lasso.R`: the two formal PMMH kernels.
- `prepare_null_lasso.R`: creates a new eligible data bank from raw attempts;
  it does not run PMMH.
- `inference_null_lasso.R`: full experiment entry point.
- `plot_inference_null_lasso.R`: ECDF entry point; it runs no simulation.
- `inference_null_lasso.rds`: frozen datasets, random streams, and paper p-values.
- `inference_null_lasso.pdf`: paper ECDF figure.

## Obtain the data

`inference_null_lasso.rds` uses Git LFS and must contain the full 34,423,726-byte
RDS (about 32.83 MiB), not a pointer. Git and Git LFS are required. To obtain a
new clone and download its LFS data:

```sh
git clone "<repository-url>" pmmh-selective-inference
cd pmmh-selective-inference
git lfs install --local
git lfs pull
```

Replace `<repository-url>` with the published URL. Keep this RDS at
`experiments/cv_lasso/inference_null/inference_null_lasso.rds` in the clone;
then use the commands below from this experiment directory. If the file is small
and starts with `version https://git-lfs.github.com/spec/v1`, run `git lfs pull`
from the clone to replace the pointer with the complete data. Download ZIP is
not guaranteed to contain full LFS data. See the [repository data instructions](../../../README.md#data-and-git-lfs).
Remote LFS upload and clean-clone download have not yet been verified.

## Experiment

Each raw dataset has `n = 200` observations and `p = 100` Gaussian design
columns with adjacent correlation 0.5. The response is `N(0, I)` and independent
Gaussian randomization has standard deviation 0.15. The design is column
centered and scaled to column norm `sqrt(n)`. Ten balanced folds and a fixed
32-point decreasing lambda grid are used.

The frozen bank contains the first 200 eligible raw attempts, in generation order, whose
observed selected model has dimension 1 through 10. Reaching this bank required
708 raw attempts. The RDS retains raw IDs, generation streams, fixed folds,
data, observed selection results, and conditional targets.

The three paired methods are:

- **Model-only**: rerun fixed-fold CV and the full-data refit, then require the
  selected model to equal the observed model. The selected lambda is not
  additionally conditioned on.
- **Fixed-lambda**: hold the observed lambda fixed, rerun the full-data Lasso,
  and require the selected model to equal the observed model.
- **Naive**: use the analytic chi-square upper tail without selection adjustment.

For PMMH, `Q = T' Sigma_T^{-1} T`. The p-value is the number of 5000 production
draws with `Q >= Q_obs`, divided by 5000, with no `+1` correction. Rejection is
`p <= 0.05`. The naive value is `P(chi-square_d >= Q_obs)` and is marked analytic.

Each chain uses 200 warm-up transitions, no stabilization, and 5000 production
transitions. There is one Bernoulli auxiliary draw, correlation 0.8, initial
scalar scale 0.8, target acceptance 0.20, and scale bounds `[0.05, 5]`. The
proposal is `T' = T + h_prop / sqrt(d) * G * xi`, where `G G' = Sigma_T`.
Therefore `Sigma_T` fixes the covariance shape and the proposal covariance is
`(h_prop^2 / d) * Sigma_T`. The scalar adapts from the MH acceptance probability
during warm-up only and is fixed during production. Rejection retains both the
target and auxiliary states.

## Environment

The recorded formal environment used R 4.4.3, glmnet 4.1.9, and Matrix 1.7.2.
Release validation used R 4.4.3, glmnet 4.1.9, and ggplot2 4.0.3. Packages are
not installed automatically.

## Run

Run from this directory.

```sh
Rscript plot_inference_null_lasso.R
Rscript inference_null_lasso.R rerun_inference_null_lasso.rds
Rscript plot_inference_null_lasso.R rerun_inference_null_lasso.rds
Rscript prepare_null_lasso.R bank_null_lasso.rds
Rscript inference_null_lasso.R new_inference_null_lasso.rds bank_null_lasso.rds
Rscript plot_inference_null_lasso.R new_inference_null_lasso.rds
```

The first command redraws the supplied paper result. The second reruns inference
on the frozen data bank contained in `inference_null_lasso.rds`. The final three
commands instead build a fresh bank from the recorded raw-attempt streams, run
inference on that bank, and plot its result. New output paths must not already
exist.

The ECDF shading is the pointwise across-dataset plus/minus one standard-error
range, `sqrt(F_hat(x) * (1 - F_hat(x)) / N)`, truncated to `[0, 1]`. It is not a
simultaneous confidence band or a within-chain Monte Carlo standard error.

The supplied object contains the paper results but not the large production
chains, so chain diagnostics cannot be recomputed from it. Assigned formal
random streams are retained for a new run; fresh-bank or PMMH reruns are not
claimed to reproduce the supplied paper values bitwise, because agreement also
depends on the recorded R, glmnet, and numerical environment.
