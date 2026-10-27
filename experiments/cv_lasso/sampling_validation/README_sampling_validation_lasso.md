# CV-Lasso: sampling accuracy validation

Compare independent rejection draws with one PM and one CPM production trajectory
for a fixed selective Gaussian target. The observed active set is
`{1, 7, 15, 30, 81, 100}`; `T1` through `T6` are its regression coefficients in
that order. The fixed design has n = 200, p = 100, randomization SD 0.15, ten fixed
folds and a 32-point lambda grid. The event is that the complete CV-and-refit
selector returns this active set, without conditioning on lambda or signs.
The frozen nuisance residual, covariance and rare Gaussian center (r = 3.7)
are included in `config$benchmark` and `config$target`.

## Files and use

- `functions_validation_lasso.R`: sample checks, statistics, numerical diagnostics
  and explicit regeneration using the sibling Efficiency samplers.
- `sampling_validation_lasso.R`: offline analysis entry; sampling is opt-in only.
- `plot_sampling_validation_lasso.R`: full-sample coordinate ECDFs, two rows by
  three columns; writes only PDF and does not recompute diagnostics.
- `sampling_validation_lasso.rds`: configuration, full samples, sample metadata,
  diagnostics and summary. `sampling_validation_lasso.pdf` is the paper figure.

Run these commands from this directory (default inputs and function sources are
also resolved relative to the script when launched elsewhere):

```sh
Rscript plot_sampling_validation_lasso.R
Rscript sampling_validation_lasso.R recomputed_lasso.rds
Rscript plot_sampling_validation_lasso.R recomputed_lasso.rds
```

The analysis accepts an optional second positional input RDS. It refuses to
overwrite any existing output. Plotting an alternate RDS writes the same basename
with `.pdf`, leaving the supplied paper figure untouched.

To actively generate new samples, using the released selector and kernels in
`../efficiency_comparison/functions_lasso.R` and `samplers_lasso.R`:

```sh
Rscript sampling_validation_lasso.R --generate regenerated_lasso.rds
Rscript plot_sampling_validation_lasso.R regenerated_lasso.rds
```

This is a full 3,000,000-proposal rejection run and two 20,000-transition PMMH
runs, not an offline diagnostic command. It saves all accepted rejection vectors
and complete production trajectories. It was inspected but **not executed**
during release preparation. No command installs packages automatically.

## Samples and numerical analysis

The reference contains 129 accepted vectors from the saved 3M rejection replay,
not a new reference run. The replay audit matched all 300 chunk hit counts and
saved start/end RNG states to the original run. Original accepted vectors were
not stored, so this does not establish samplewise equality with an unavailable
original sample file. The replay vectors are exactly those used by both the
authoritative numerical analysis and the paper coordinate plot.

PM/CPM come from the completed formal checkpoints. Each has 500 adaptive
transitions, no stabilization, and 19,500 production states (iterations
501--20,000). All holding states are retained; warm-up is removed once. Each
transition costs one selector evaluation; the formal cost is 20,000, including
adaptation but excluding initializer verification. PM has independent auxiliary
proposals; CPM has rho = 0.8; the Bernoulli estimate uses one auxiliary draw.
Both start at the selected observation and its observed auxiliary state.
The proposal covariance is `(h_prop^2 / d) * Sigma_T`, initially h_prop = 0.8;
only warm-up adapts toward 0.20 using acceptance probabilities, with bounds
0.05--5. Saved L'Ecuyer-CMRG states and chunk budgets are retained. These allow
new sampling in the same design; exact historical replay across numerical
environments has not been established by executing the released sampler.

`maximum_holding_run` is the largest number of saved production states in a
consecutive run of exactly identical full vectors: 328 for PM and 57 for CPM.
All coordinates are compared directly, without rounding or tolerance. A run
includes its first state; no warm-up or preceding state is included. This is
not a count of consecutive rejected transitions or a run of acceptance flags.

Diagnostics include all coordinates, uncentered `Q = T' Sigma_T^-1 T`, and
uncentered directional `S_proj`; their definitions and the paper plot are fixed.
Both applications use `draws <- matrix(as.numeric(x), ncol = 1L)`, with iterations
in rows and the single chain in one column, then
`posterior::ess_basic(draws, split = TRUE)` and `posterior::mcse_mean(draws)`.
`ess_mean` is the raw-statistic mean ESS; valid results satisfy
`mcse_mean = sd(x) / sqrt(ess_mean)` and `tau = length(x) / ess_mean`.
Coordinate `ess_bulk` is retained separately for the Efficiency comparison and
never used in this mean-MCSE conversion. No additional ESS cap or alternative
estimator fallback is applied; package numerical protections and warnings remain
visible and are recorded in `package_warnings`.

Upper-tail (`>=`) diagnostics apply these same public functions to each specified
indicator's own 0/1 trajectory, not to the original continuous variable. Successes,
failures and state transitions accompany indicator ESS/MCSE. Constant, too-short
or otherwise unestimable sequences have NA ESS/MCSE and an `undefined_reason`,
never a fabricated zero MCSE. Nonfinite input is an error identifying the method
and statistic. All eight saved MCMC tail indicators have finite ESS/MCSE.

The independent reference mean SE is SD/sqrt(129). Conditional on the fixed
target, separate assigned RNG streams with no saved starting-state reuse are
the basis for treating reference and chain errors as independent and adding
their variances for mean comparisons. This is an RNG-design assumption, not an
independence test. Rejection draws never enter the MCMC ESS helper. Tail
thresholds are rejection empirical quantiles (0.75 and 0.90, type 7), so chain
indicator MCSE is conditional on the realized threshold and excludes threshold
uncertainty. No independent binomial reference interval or Z test is attached
to these data-dependent thresholds. Mixing, half-chain means, moments, quantiles,
covariance comparisons and descriptive maximum ECDF discrepancies are retained.

The reference has 129 draws and each kernel contributes one production trajectory.
Minimum coordinate mean ESS is 58.080104 (PM) and 152.968611 (CPM). Splitting a
trajectory does not provide independently initialized chains or multi-chain
convergence evidence. Coordinate ECDFs are descriptive marginal comparisons,
not proof of joint agreement or accurate tail inference.

## Environment and verification

Recorded original numerical-analysis environment: R 4.4.3, posterior 1.7.0.
Original sampling R/glmnet versions are unknown. Offline recalculation used
R 4.4.3, posterior
1.7.0, ggplot2 4.0.3 and glmnet 4.1-9. Offline analysis needs posterior; plotting
needs ggplot2 and Cairo PDF support; explicit regeneration also needs glmnet and
the sibling Efficiency code. Historical run time is unknown.

Offline checks cover syntax, function-only loading, exact-vector holding-run
examples, equality of retained diagnostics, and sample/target/metadata invariance.
The posterior wrappers were also checked against direct calls and the ESS/MCSE
identity, including degenerate and nonfinite inputs. The paper PDF and coordinate
ECDFs remain unchanged. No MCMC, rejection sampling, selector evaluation or smoke
simulation was run. Frozen Efficiency results and their accounting were not changed.

