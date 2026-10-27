# Changepoint directional power

Full-selector, Full-path and matched Data Fission are compared at delta=0.4,
0.6, 0.8, 1.0. Each signal has 200 paired eligible datasets. The supplied results
combine two complementary sets of 100 from the same frozen bank: IDs 1--100 and
101--200. Public preparation and inference use one bank and one task loop;
source_batch and the actual task seeds retain their identities.

## Model, events, and directional target

Y=delta*b+epsilon, epsilon~N(0,I_200), where b is one on indices 101:140 and zero
elsewhere. Independent W~N(0,0.1^2 I) gives selection input Y+W. Bounded binary
segmentation makes four sequential splits, with minimum segment length 10.
Candidate ties use segment start, split location, then segment end. Adaptive K
minimizes RSS_K+1.25*K*log(200) over K=0,...,4; ties within 1e-12 choose smallest K.
Each delta retains the first 200 raw-order datasets with K_obs>=1. No dataset is
filtered by p-value or chain behavior. dataset_id is local to delta; raw_id is the
original attempt.

Full-selector matches the final unordered adaptive changepoint set. Full-path
matches observed K, all four ordered split locations and all four signed-CUSUM
signs, including splits after K_obs. It does not condition on the IC ranking or
additional pairwise IC constraints. Both rerun the complete adaptive selector.

V contains adjacent segment-mean contrasts: negative reciprocal segment length
on the left, positive reciprocal length on the right. T=V'Y, Sigma=V'V and
R=Y-V*Sigma^(-1)T. Fix R and a=V'b, with g=a'Sigma^(-1)a. The statistic is
S_dir(t)=a'Sigma^(-1)t/sqrt(g). Its observed direction is fixed throughout the
chain and uses the known signal shape b, not an estimated unknown signal.
The multivariate null reference is N(0,Sigma), with the relevant event and R
fixed; the directional test is upper-tailed at projected null mean zero. Delta
generates alternatives and is not inserted into the null reference mean.
The recorded structural-zero convention is p=1 without sampling; no supplied
dataset has a structural-zero direction.

Fission selects on Y+W and infers from Y_fiss=Y-100W. The Gaussian views are
independent, with covariances 1.01 I and 101 I. T_fiss=V'Y_fiss has covariance
101 Sigma. The analytic p-value is the normal upper tail of
a'Sigma^(-1)T_fiss/sqrt(101*g). The same observed changepoints and direction
are used, without reselection, extra random draws or a fictitious chain.

## PMMH and power

Correlated PMMH initializes at observed T and W/0.1, with 200 adaptive warm-up
transitions, zero stabilization, 5000 production draws, one Bernoulli auxiliary
estimate, rho=0.8, initial h=0.8, target acceptance 0.20 and bounds [0.05,5].
Warm-up uses log(h) += 0.5/(i+10)^0.6*(MH acceptance probability - 0.20);
production freezes h. Rejection retains both target and auxiliary. Cost is one
selector call per transition (5200), with the observed initializer taken as given.

The frozen runner uses T'=T+(h/sqrt(K))*R_chol*xi with R_chol=chol(Sigma),
without transposing the upper factor. Its actual proposal covariance is
(h^2/K)*R_chol*R_chol', where R_chol'*R_chol=Sigma. This recorded implementation
is preserved; its symmetric proposal still uses the Gaussian null target in MH.
It should not be described as proposal covariance (h^2/K)*Sigma.

P-values are inclusive upper-tail counts divided by 5000, without +1 correction;
reject iff p<=0.05. Power is recomputed from the 200 paired rejection indicators.
Error bars are power +/- sqrt(power*(1-power)/200), not single-batch SE or MCMC
precision. Endpoints are not numerically clipped; the plot keeps its [0,1] scale.

## Files and environment

- functions_nonnull_cp.R: configuration, DGP, selector/events, eligibility,
  target/direction, fission, RNG assignments, checks and summaries.
- samplers_nonnull_cp.R: the formal correlated PMMH transitions.
- prepare_nonnull_cp.R: constructs the complete eligible bank without PMMH.
- inference_nonnull_cp.R: paired inference on a supplied or frozen bank.
- plot_inference_nonnull_cp.R: reads RDS and uses ../../../plot_style.R.
- inference_nonnull_cp.rds: frozen data/targets, generation/task seeds and states,
  dataset results, power and saved finite-chain summaries.
- inference_nonnull_cp.pdf: the paper figure.

Sampling and preparation use base R and parallel. Plotting requires ggplot2,
patchwork and Cairo. Original sampler software versions are unknown; the recorded
results-analysis environment used R 4.4.3 and ggplot2 4.0.3. Current checks used
R 4.4.3, ggplot2 4.0.3 and patchwork 1.3.0. No packages are installed automatically.

## Reproduction

Run from this directory. Redraw paper results:

```sh
Rscript plot_inference_nonnull_cp.R
```

Rerun inference on the frozen bank:

```sh
Rscript inference_nonnull_cp.R rerun_inference_nonnull_cp.rds
Rscript plot_inference_nonnull_cp.R rerun_inference_nonnull_cp.rds
```

Build a new bank, run inference, and plot:

```sh
Rscript prepare_nonnull_cp.R bank_nonnull_cp.rds
Rscript inference_nonnull_cp.R new_inference_nonnull_cp.rds bank_nonnull_cp.rds
Rscript plot_inference_nonnull_cp.R new_inference_nonnull_cp.rds
```

New bank/result paths must not exist. An explicit bank supplies its own config,
datasets and task streams. Generation uses the recorded per-delta master and
raw-attempt seed formula; inference retains the two complementary seed ranges.
Fission needs no additional random stream.

The supplied p-values and decisions are the paper results. Small deterministic
generation/target/event checks and complete saved-result/count checks were
performed. Full bank reconstruction and PMMH reruns were not executed for release;
bitwise recovery also depends on the R/numerical environment. Production chains
are not bundled, so chain diagnostics cannot be recomputed from this RDS. Saved
summaries include constant tail indicators, threshold overlap and prefix changes;
they limit finite-chain precision claims. Integrity does not establish convergence.
Full-run time is unknown.

