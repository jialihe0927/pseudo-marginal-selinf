# CV-Lasso directional power

Model-only, Full-conditioning, and matched Data Fission are compared on the same
200 eligible datasets at each delta = 0.2, 0.3, 0.4, 0.5.

## Obtain the data

`inference_nonnull_lasso.rds` uses Git LFS and must contain the full
140,486,247-byte RDS (about 133.98 MiB), not a pointer. Git and Git LFS are
required. To obtain a new clone and download its LFS data:

```sh
git clone "<repository-url>" pmmh-selective-inference
cd pmmh-selective-inference
git lfs install --local
git lfs pull
```

Replace `<repository-url>` with the published URL. Keep this RDS at
`experiments/cv_lasso/inference_nonnull/inference_nonnull_lasso.rds` in the clone;
then use the plotting and inference commands below from this experiment
directory. If the file is small and starts with
`version https://git-lfs.github.com/spec/v1`, run `git lfs pull` from the clone
to obtain the complete data. Download ZIP is not guaranteed to contain full
LFS data. See the [repository data instructions](../../../README.md#data-and-git-lfs).
Remote LFS upload and clean-clone download have not yet been verified.

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
the formal uniform draw occurs only when the proposed event is positive.
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
Dependencies are not installed automatically.

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

The supplied results are the paper values. Small deterministic generation/target
checks and complete saved-result checks were performed; full bank reconstruction
and PMMH reruns were not executed for this release. Bitwise agreement also depends
on glmnet and the numerical environment. Large production chains are not bundled,
so chain diagnostics cannot be recomputed from this RDS. Saved finite-chain
summaries include constant tails, threshold overlap and prefix decision changes.
These limit precision claims; integrity checks do not establish convergence.
No dataset or decision was changed. Full-run time is unknown.

