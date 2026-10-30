# Changepoint global-null calibration

This experiment evaluates selective p-value calibration after randomized
adaptive changepoint selection under the global null.

## Files

- `functions_null_cp.R`: data generation, selectors, selection events,
  conditional reconstruction, data-bank construction, and validation.
- `samplers_null_cp.R`: the two formal PMMH kernels.
- `prepare_null_cp.R`: creates a new eligible data bank from raw attempts; it
  does not run PMMH.
- `inference_null_cp.R`: full experiment entry point.
- `plot_inference_null_cp.R`: ECDF entry point; it runs no simulation.
- `inference_null_cp.rds`: frozen datasets, random streams, and paper p-values.
- `inference_null_cp.pdf`: paper ECDF figure.

## Experiment

Each raw series has length 200 with `Y ~ N(0, I_200)`. Independent Gaussian
randomization has standard deviation 0.1 and the selector receives `Y + W`.
Bounded binary segmentation uses minimum segment length 10 and a four-split
path. The adaptive selector chooses `K` from 0 through 4 by minimizing
`RSS_K + 1.25 K log(n)`; ties within `1e-12` use the smallest K.

The frozen bank contains the first 200 eligible raw attempts, in generation order, with
observed `K >= 1`. Reaching this bank required 881 raw attempts. The RDS retains
raw IDs, generation provenance, data, observed selection results, and targets.

The three paired methods are:

- **Full-selector**: rerun the complete adaptive-K selector and require its
  final unordered changepoint set to equal the observed set.
- **Fixed-K**: fix observed K, take exactly K binary-segmentation additions,
  and require the final unordered set to equal the observed set.
- **Naive**: use the analytic chi-square upper tail without selection adjustment.

For PMMH, `Q = T' Sigma_T^{-1} T`. The p-value is the number of 5000 production
draws with `Q >= Q_obs`, divided by 5000, with no `+1` correction. Rejection is
`p <= 0.05`. The naive value is `P(chi-square_K >= Q_obs)` and is marked analytic.

Each chain uses 200 warm-up transitions, no stabilization, and 5000 production
transitions. There is one Bernoulli auxiliary draw, correlation 0.8, initial
scalar scale 0.8, target acceptance 0.20, and scale bounds `[0.05, 5]`. The
proposal is `T' = T + h_prop / sqrt(d) * G * xi`, where `G G' = Sigma_T`.
Therefore `Sigma_T` fixes the covariance shape and the proposal covariance is
`(h_prop^2 / d) * Sigma_T`. The scalar adapts from the MH acceptance probability
during warm-up only and is fixed during production. Rejection retains both the
target and auxiliary states.

## Environment

R 4.4.3 is the verified environment; indicator diagnostics require posterior
1.7.0 and warn on other R versions. Bank preparation uses base R and parallel
without a diagnostic version lock. Plotting uses ggplot2 4.0.3 and Cairo.
Packages are not installed automatically.

## Run

Commands below run from this directory. Default inputs and dependencies are
script-relative; explicit relative paths use the caller's directory. Sourcing
an entry point does not execute its main program. Redraw without sampling:

```sh
Rscript plot_inference_null_cp.R
```

Run complete inference on the frozen bank, then plot the new results:

```sh
Rscript inference_null_cp.R rerun_inference_null_cp.rds
Rscript plot_inference_null_cp.R rerun_inference_null_cp.rds
```

Generate a complete bank from raw attempts, run inference, and plot:

```sh
Rscript prepare_null_cp.R bank_null_cp.rds
Rscript inference_null_cp.R new_inference_null_cp.rds bank_null_cp.rds
Rscript plot_inference_null_cp.R new_inference_null_cp.rds
```

Preparation and inference start full computations, not shortened examples.
The default bank is contained in `inference_null_cp.rds`; an explicit bank
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
depends on the recorded R and numerical environment.
