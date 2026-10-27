# Changepoint: sampling accuracy validation

Compare independent rejection draws with one PM and one CPM production trajectory
for a fixed selective Gaussian target. The observed unordered changepoint set is
`{71, 99, 140}`. `T1`, `T2`, `T3` are the successive right-minus-left segment-mean
contrasts in that order. The event is that the complete adaptive selector returns
this set: bounded binary segmentation with h = 10, K = 0,...,4 and penalty
`1.25 * K * log(n)`. Randomization SD is 0.1. The conditional residual, contrast
matrix, covariance and rare Gaussian center (r = 4) are included in the RDS.

## Files and use

- `functions_validation_cp.R`: sample checks, statistics, numerical diagnostics
  and explicit regeneration using the sibling Efficiency samplers.
- `sampling_validation_cp.R`: offline analysis entry; sampling is opt-in only.
- `plot_sampling_validation_cp.R`: full-sample coordinate ECDFs, one row by three
  columns; writes only PDF without recomputing diagnostics.
- `sampling_validation_cp.rds`: configuration, full samples, sample metadata,
  diagnostics and summary. `sampling_validation_cp.pdf` is the paper figure.

Run from this directory (default inputs and function sources are also resolved
relative to the script when launched elsewhere):

```sh
Rscript plot_sampling_validation_cp.R
Rscript sampling_validation_cp.R recomputed_cp.rds
Rscript plot_sampling_validation_cp.R recomputed_cp.rds
```

The analysis accepts an optional second positional input RDS and refuses existing
output files. Plotting an alternate RDS writes a PDF with its basename, leaving
the supplied paper figure untouched.

To actively regenerate samples using the released selector and kernels in
`../efficiency_comparison/functions_cp.R` and `samplers_cp.R`:

```sh
Rscript sampling_validation_cp.R --generate regenerated_cp.rds
Rscript plot_sampling_validation_cp.R regenerated_cp.rds
```

This runs 3,000,000 rejection proposals and two 20,000-transition PMMH chains,
and saves accepted rejection vectors and complete production trajectories.
The regeneration path was inspected but **not executed** during preparation.
No packages are installed automatically.

## Samples and diagnostics

The reference is the 78 accepted vectors from the recorded primary 3M rejection
run, with original candidate identities and 70 chunk RNG states. It is exactly
the reference read by the authoritative analysis and coordinate plot; no batches
were added. PM/CPM are the completed formal chains: 500 adaptive transitions,
zero stabilization and 19,500 production states at iterations 501--20,000.
All holding states are retained, with warm-up removed once. One selector call is
charged per transition; the 20,000-call cost includes adaptation but not initial
selection verification. PM uses independent auxiliary proposals; CPM uses
rho = 0.8; one auxiliary draw gives the Bernoulli estimate.

Initialization is the selected observation and its observed auxiliary state.
Proposal covariance is `(h_prop^2 / d) * Sigma_T`, initially h_prop = 0.8. Only
warm-up adapts toward acceptance 0.20 using acceptance probabilities, with bounds
0.05--5. L'Ecuyer-CMRG chain seeds are 2026100412/413, auxiliary seeds
2026101412/413 and proposal seeds 2026102412/413; rejection master seed is
2026100501. Explicit generation uses the saved chunk states and these separate
streams, not the Efficiency entry's default seed offsets. Exact cross-environment
historical trajectory reproduction has not been verified by rerunning sampling.

`maximum_holding_run` is the largest number of saved production states in a
consecutive run of exactly identical full vectors: 84 for PM and 53 for CPM.
All coordinates are compared directly, without rounding or tolerance. A run
includes its first state and excludes warm-up and any preceding state. It is
not the number of consecutive rejected transitions or a run of held flags.

Diagnostics retain coordinates, uncentered `Q = T' Sigma_T^-1 T`,
`S_shift = a' Sigma_T^-1 (T - theta) / sqrt(a' Sigma_T^-1 a)`, and supplementary
centered Q. Here theta is the frozen Gaussian center and a = sign(T_obs[1]) e1.
Row-wise centering is essential and has not changed. Both applications use
`draws <- matrix(as.numeric(x), ncol = 1L)`, with iterations in rows and the single
chain in one column, then `posterior::ess_basic(draws, split = TRUE)` and
`posterior::mcse_mean(draws)`. The resulting `ess_mean` is mean ESS on the raw
statistic, not rank-normalized bulk ESS. Valid results satisfy
`mcse_mean = sd(x) / sqrt(ess_mean)`; `tau = length(x) / ess_mean` is also saved.
There is no handwritten ESS estimator, extra ESS cap, or batch-means fallback.
The separate descriptive ACF retains its existing autocovariance calculation;
it does not supply ESS. Efficiency's frozen bulk ESS is not modified or used
to convert mean MCSE.

Each tail indicator is analyzed as its own 0/1 trajectory by these same public
interfaces. Successes, failures and state transitions are reported without
inventing an adequacy threshold. All eight MCMC tail indicators are estimable,
but one has only 18 state transitions; finite ESS is not a precision guarantee.
Constant/too-short/unestimable sequences retain NA ESS/MCSE and a reason, not a
zero MCSE or an estimator switch. Nonfinite input identifies the method/statistic
and stops analysis. Reasons are saved in `undefined_reason`, and original package
warnings remain visible and are saved in `package_warnings`. Splitting a
trajectory does not provide independently initialized chains or multi-chain
convergence evidence.

The independent reference mean SE is SD/sqrt(78). Conditional on the fixed target,
separately assigned RNG streams without saved starting-state reuse are the basis
for adding reference and chain variances. This is an independence assumption of
the RNG design, not an independence test. Rejection samples never enter the MCMC
ESS helper.

Tail thresholds are the pre-selection noncentral chi-square Q and standard-normal
S_shift quantiles at 0.75/0.90. Independent-reference boundary counts retain
Wilson intervals and the descriptive `reference_boundary` flag; boundary counts
are not treated as exact probabilities. Moments, quantiles, ECDF discrepancies, covariance/correlation,
autocorrelations, rolling means and mixing summaries are retained. The existing
covariance comparison uses 500 descriptive bootstrap replicates: IID reference
resampling and circular blocks of length 200 for chains (seed 20261022).
Whitening uses the Gaussian target covariance, not a claim that this equals the
selective covariance. No new diagnostic or interval was introduced.

The reference has 78 draws and each kernel contributes one production trajectory.
Minimum coordinate mean ESS is 238.955843 (PM) and 309.020906 (CPM). At the
S_shift 0.75 threshold, reference exceedances are 78/78, with a Wilson interval
of approximately [0.953, 1]. One MCMC indicator has only 18 state transitions.
Coordinate ECDFs compare marginals descriptively; they do not prove joint
agreement or uniformly accurate tail inference.

## Environment and verification

Historical R/package versions are unknown. Offline recalculation used R 4.4.3,
posterior 1.7.0; plotting was checked with ggplot2 4.0.3. Offline diagnostics need
posterior; plotting needs ggplot2 and Cairo PDF support. Explicit generation also
requires the sibling Efficiency code. Historical run time is unknown.

Offline checks cover syntax, function-only loading, exact-vector holding-run
examples, equality of retained diagnostics, and sample/target/metadata invariance.
The posterior wrappers were also checked against direct calls and the ESS/MCSE
identity, including degenerate and nonfinite inputs. Existing covariance bootstrap
estimates and intervals were retained without recomputation. Descriptive ACF,
rolling outputs, coordinate ECDFs and the paper PDF remain unchanged. No MCMC,
rejection sampling, selector evaluation or smoke simulation was run. Efficiency
results and accounting were not changed.

