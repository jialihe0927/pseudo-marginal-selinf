# CV-Lasso directional power

Proposed, Extra conditioning, and matched Data splitting are compared on the same
200 eligible datasets at each delta = 0, 0.2, 0.3, 0.4, 0.5. Internal method IDs
are `model_only`, `full_conditioning`, and `data_fission`, respectively.

## Obtain the data

`inference_nonnull_lasso.rds` uses Git LFS. Follow the
[repository download instructions](../../../README.md#get-the-code-and-lfs-data)
and keep the complete RDS beside these scripts. A small LFS pointer is not usable
as an RDS. Plotting and inference then use the commands below.

## Model, selection, and target

Each raw dataset has n=200, p=100 and a Gaussian AR(1) design with correlation
0.5, column-centered and scaled to norm sqrt(n). Y=X(delta v)+epsilon with
epsilon~N(0,I) and fixed v=(1,-0.8,0.6,-0.5,0.4,0,...). Independent
W~N(0,0.15^2 I) is added for selection. Balanced ten-fold assignments are fixed
within each dataset. The selector minimizes CV deviance on a fixed 32-lambda
grid (first minimum), then refits on all data, without intercept or additional
standardization. The first 200 raw-order eligible datasets with dimension 1--10
and a valid target are retained. dataset_id is local to delta; raw_id identifies
the generating attempt. The selector is replayed during preparation.
The implementation requires unit response noise (`response_sd = 1`).
Lambda uses [glmnet's normalization](../README_lasso.md).

Model-only matches the observed model. Full-conditioning also matches the
observed lambda index and coefficient signs. Both rerun complete CV and refit.
For the selected model M, H=(X_M'X_M)^(-1)X_M', T=HY, Sigma=(X_M'X_M)^(-1).
The nuisance residual uses C=H'(HH')^(-1), Z=Y-CT; its formal numerical path is
retained. The fixed direction is a=(X_M'X_M)^(-1)X_M'Xv. The statistic is
S(t)=a'Sigma^(-1)t/sqrt(a'Sigma^(-1)a). The multivariate null reference is
N(0,Sigma), conditional on Z and the method's selection event; the directional
test uses the upper tail at projected null mean zero. The nonzero DGP delta is
not inserted into the null reference mean.

Fission selects on Y+W and infers from Y_inf=Y-W/0.15^2. Conditional on X these
Gaussian views are independent, with covariances 1.0225 I and
(1+0.15^(-2))I. It uses the same observed model and direction, with no reselection
or new random draws. Its analytic p-value is the normal upper tail of
S(HY_inf)/sqrt(1+0.15^(-2)); the formal direct normal-equation solve is retained.

## PMMH and power

Correlated PMMH starts at the observed feasible T and W/0.15. It uses 200 adaptive
warm-up transitions, zero stabilization, and 5000 fixed-scale production draws,
one Bernoulli auxiliary estimate, rho=0.8, initial scalar h=0.8, bounds [0.05,5],
and target acceptance 0.20. The warm-up update is
log(h) += 0.5/(i+10)^0.6 * (MH acceptance probability - 0.20).
Proposal covariance is (h^2/d)Sigma. Rejection retains target and auxiliary;
the positive-signal and Extra-conditioning uniform draw occurs only when the
proposed event is positive. Delta-zero Proposed reuses the quadratic-null
trajectory, whose runner draws its uniform before the selector; the rerun branch
uses that runner and RNG assignment rather than the positive-signal kernel.
Selector cost is 5200 transition calls plus one initializer check per chain.

The inclusive upper-tail p-value is tail_count/5000, without +1 correction;
reject iff p<=0.05. Power is the mean paired rejection indicator. Error bars are
power +/- sqrt(power*(1-power)/N), N=200 per signal/method. Endpoints are not
numerically clipped; the figure retains its [0,1] scale. These are across-dataset
one-standard-error bars, not chain-level precision.

## Files and environment

- functions_nonnull_lasso.R: configuration, DGP, selectors/events, eligibility,
  target/direction, fission, bank/result checks, and power summaries.
- samplers_nonnull_lasso.R: the formal PMMH transitions.
- prepare_nonnull_lasso.R: generates an eligible bank without PMMH.
- inference_nonnull_lasso.R: paired inference on a new or frozen bank.
- plot_inference_nonnull_lasso.R: reads RDS and uses ../../../plot_style.R.
- inference_nonnull_lasso.rds: frozen data, targets/directions, generation/task
  streams, individual p-values/counts, power, and saved finite-chain summaries.
- inference_nonnull_lasso.pdf: the paper figure.

Preparation/inference requires R and glmnet; plotting requires ggplot2 and
patchwork with Cairo. Current checks used R 4.4.3, glmnet 4.1.9, Matrix 1.7.2,
ggplot2 4.0.3 and patchwork 1.3.0. Original execution versions are unknown.
Dependencies are not installed automatically. Diagnostics require posterior 1.7.0.
New computations warn if R differs from the verified R 4.4.3; bank preparation
and plotting have no diagnostic version gate.

To load the functions from this directory without running an experiment:

```r
source("../inference_null/functions_null_lasso.R")
source("functions_nonnull_lasso.R")
```

## Reproduction

Run from this directory. Redraw supplied paper results:

```sh
Rscript plot_inference_nonnull_lasso.R
```

Rerun inference on the frozen bank:

```sh
Rscript inference_nonnull_lasso.R rerun_inference_nonnull_lasso.rds
Rscript plot_inference_nonnull_lasso.R rerun_inference_nonnull_lasso.rds
```

Build a new bank, run inference, and plot:

```sh
Rscript prepare_nonnull_lasso.R bank_nonnull_lasso.rds
Rscript inference_nonnull_lasso.R new_inference_nonnull_lasso.rds bank_nonnull_lasso.rds
Rscript plot_inference_nonnull_lasso.R new_inference_nonnull_lasso.rds
```

New bank/result paths must not exist. An explicit bank supplies all configuration,
data and task streams. Each generation stream yields four component seeds;
Mersenne-Twister then generates design, noise, randomization and folds separately.
Recorded inference stream positions are preserved independently of scheduling.
Analytic fission requires no task stream.

Delta-zero reuse checks n, p, design correlation, response/randomization scales,
fold count, the exact lambda grid, coefficient and glmnet thresholds, dimension
range and RNG kinds against the null implementation. Other combinations are
rejected before preparation or sampling. PMMH lengths, correlation, acceptance
target and scalar bounds are transferred through explicit field mappings.

The supplied results are the paper values. Full bank reconstruction and PMMH
reruns have not been runtime-tested; bitwise agreement also depends on glmnet
and the numerical environment. Full target/auxiliary chains are not bundled.
Complete ordered production statistics are stored once per physical
trajectory in the root `inference_trajectories.rds`; `Rscript ../../../reproduce.R`
recomputes their diagnostics without the research workspace. Saved finite-chain
summaries include constant tails, threshold overlap and prefix decision changes.
These limit precision claims and do not establish convergence. Full-run time is unknown.

At delta zero, Proposed and Data splitting use the same 200 global-null datasets
as quadratic calibration (last retained raw ID 708). Proposed reuses its existing
production trajectories; Extra conditioning uses 200 separately completed chains
matching the final active set, CV lambda index and active refit signs. The
nonzero direction template remains `a=H_M X v` even when delta=0. Rejection
rates are 0.060, 0.085 and 0.065. Preparation uses the quadratic-null raw-attempt
streams for delta zero. Extra-conditioning streams start from L'Ecuyer-CMRG seed
20261028; positive-signal streams retain their assigned identities.

Diagnostics use each inclusive upper-tail indicator separately with
`posterior::ess_basic(matrix(as.numeric(I), ncol=1L), split=TRUE)` and
`posterior::mcse_mean(...)`, posterior 1.7.0 / R 4.4.3. Nonfinite estimates
remain missing and are excluded separately from the corresponding pooled median;
no p-value is excluded. Splitting one trajectory is not independent multi-chain
initialization. The figure includes delta zero and a horizontal 0.05 reference.
`summary$finite_chain` contains five-signal finite-count and prefix summaries.

