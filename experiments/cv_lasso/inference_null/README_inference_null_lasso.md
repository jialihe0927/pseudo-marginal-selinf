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

`inference_null_lasso.rds` uses Git LFS. Follow the
[repository download instructions](../../../README.md#get-the-code-and-lfs-data)
and keep the complete RDS beside these scripts. A small LFS pointer is not usable
as an RDS. Plotting and inference then use the commands below.

## Experiment

Each raw dataset has `n = 200` observations and `p = 100` Gaussian design
columns with adjacent correlation 0.5. The response is `N(0, I)` and independent
Gaussian randomization has standard deviation 0.15. The design is column
centered and scaled to column norm `sqrt(n)`. Ten balanced folds and a fixed
32-point decreasing lambda grid are used. The implementation requires unit
response noise (`response_sd = 1`); unsupported values are rejected before computation.
Lambda follows [glmnet's normalization](../README_lasso.md), including the
training-sample-size scaling in CV.

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

New inference requires posterior 1.7.0 for tail-indicator diagnostics; R 4.4.3 is
the verified environment, with a warning on other R versions. Bank preparation
and plotting have no posterior version requirement. No packages are installed automatically.

Run from this directory. Redraw the frozen results without sampling:

```sh
Rscript plot_inference_null_lasso.R
```

Run complete inference on the frozen bank, then plot the new results:

```sh
Rscript inference_null_lasso.R rerun_inference_null_lasso.rds
Rscript plot_inference_null_lasso.R rerun_inference_null_lasso.rds
```

Generate a complete bank from raw attempts, run inference, and plot:

```sh
Rscript prepare_null_lasso.R bank_null_lasso.rds
Rscript inference_null_lasso.R new_inference_null_lasso.rds bank_null_lasso.rds
Rscript plot_inference_null_lasso.R new_inference_null_lasso.rds
```

Preparation and inference start full computations, not shortened examples.
The default bank is contained in `inference_null_lasso.rds`; an explicit bank
supplies its own configuration, datasets and streams. New output paths must
not already exist.

The ECDF shading is the pointwise across-dataset plus/minus one standard-error
range, `sqrt(F_hat(x) * (1 - F_hat(x)) / N)`, truncated to `[0, 1]`. It is not a
simultaneous confidence band or a within-chain Monte Carlo standard error.

The supplied object contains the paper results but not the large production
target/auxiliary chains. The root `inference_trajectories.rds` supplies complete
ordered Q statistics and the reused directional statistics under shared physical
trajectory keys. `Rscript ../../../reproduce.R` recomputes their actual tail
indicator diagnostics offline using posterior 1.7.0 and R 4.4.3. The same Proposed
production trajectory supplies two different tests at delta zero, not two runs.
Each indicator separately uses `ess_basic(..., split=TRUE)` and `mcse_mean` on
an iterations-by-one-chain matrix. Missing ESS and MCSE remain missing and are
excluded separately from their own medians; no p-value is omitted. Assigned formal
random streams are retained for a new run; fresh-bank or PMMH reruns are not
claimed to reproduce the supplied paper values bitwise, because agreement also
depends on the recorded R, glmnet, and numerical environment.
