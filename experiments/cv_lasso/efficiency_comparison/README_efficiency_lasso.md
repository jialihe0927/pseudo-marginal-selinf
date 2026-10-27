# CV-Lasso sampling efficiency comparison

This experiment compares rejection sampling, Crude PMMH (independent auxiliary
refreshment), and Correlated PMMH for a selective Gaussian target induced by a
fixed-fold CV-Lasso selector.

## Files and dependencies

- `functions_lasso.R`: frozen-fold selector, selection event, and conditional reconstruction.
- `samplers_lasso.R`: rejection and the two PMMH kernels.
- `sampling_efficiency_lasso.R`: formal simulation entry point.
- `table_efficiency_lasso.R`: reads the RDS and prints the frozen paper table.
- `efficiency_lasso.rds`: frozen benchmark and paper estimates.

Packages are not installed automatically. The original formal-run versions of
R, `glmnet`, and `posterior` were not recorded in the available run artifacts
and are therefore unknown. This release was checked with R 4.4.3, glmnet 4.1-9,
and posterior 1.7.0. These current versions are validation information, not a
claim about the historical environment.

## Frozen setting

The benchmark has 200 observations, 100 predictors, 10 fixed balanced folds,
and a frozen 32-value lambda grid. The observed six-variable model is the
conditioning event; lambda index and signs are not conditioned on. The rare
Gaussian center is displaced by 3.7 Mahalanobis standard deviations. Both PMMH
chains use one auxiliary draw and propose
`T' = T + h_prop / sqrt(d) * G * xi`, where `G G^T = Sigma_T` and
`xi ~ N(0, I_d)`. Thus `Sigma_T` fixes the covariance shape and the actual
proposal covariance is `(h_prop^2 / d) * Sigma_T`. The scalar `h_prop` starts at
0.8, adapts only during the 500 warm-up transitions toward acceptance 0.20, and
is fixed for the 19,500 production transitions. Correlated PMMH uses auxiliary
correlation 0.8; crude PMMH uses independent refreshment (`rho = 0`). Each
transition costs one selector evaluation,
and efficiency is the minimum coordinate bulk ESS divided by all 20,000 calls,
times 1,000. Rejection uses exactly 3,000,000 candidate attempts in 300 fixed
10,000-attempt RNG streams; the original run used seven workers. Three million
is not the number of accepted samples.

The bundled summaries are the paper results. The fixed benchmark is included,
but a new run has ordinary Monte Carlo variation and is not promised to restore
the historical random stream bit for bit.

The paper table reports `(h_final, efficiency, gain, acceptance, zero rate)` as
`(--, 0.043, 1.0, --, --)`, `(0.524, 3.260, 75.8, 0.193, 0.717)`, and
`(0.718, 8.215, 191.1, 0.220, 0.638)` for rejection, crude PMMH, and correlated
PMMH, respectively.

Here `h_final` is the proposal scale at the end of warm-up, held fixed during
production. `Acc.` is the actual production acceptance fraction. `Zero rate` is
the fraction of production proposals with a zero Bernoulli selection estimate;
it is not the holding fraction, since a positive estimate can still be rejected
by MH. Each formal PM/CPM budget is 20,000 selector evaluations, including
warm-up but excluding the extra initial selection verification.

Run from this directory:

```sh
Rscript table_efficiency_lasso.R
Rscript sampling_efficiency_lasso.R pmmh rerun_pmmh_lasso.rds
Rscript sampling_efficiency_lasso.R rejection rerun_rejection_lasso.rds
Rscript sampling_efficiency_lasso.R all rerun_lasso.rds
Rscript table_efficiency_lasso.R rerun_lasso.rds
```

The entry point refuses to overwrite an existing result. Only `all` produces a
complete three-method object that can reproduce the paper table. Original
wall-clock time is unknown. Before simulation, the entry point checks the
frozen transition accounting, the 3,000,000-attempt rejection budget, the
Bernoulli auxiliary-sample requirement, proposal-scale bounds, and correlation.
The full simulation was not executed during release validation.

