# Changepoint directional inference

Compare Proposed (`full_selector`) and matched Data splitting (`data_fission`)
at delta = 0, 0.4, 0.6, 0.8, 1, using 200 paired eligible datasets per signal.

## Model and target

`Y = delta*b + epsilon`, with `epsilon ~ N(0,I_200)` and `b = 1{101<=i<=140}`.
Independent `W ~ N(0,0.1^2 I)` gives selection input `Y+W`. The selector builds
four bounded segmentation splits, minimum segment length 10; candidates are
ordered by absolute CUSUM, segment start, split, then segment end. It minimizes
`RSS_K + 1.25*K*log(n)` over K = 0,...,4, choosing the smallest K within tie
tolerance 1e-12. Proposed matches only the final unordered changepoint set.
The complete four-step path remains necessary for choosing K.

Each signal retains the first 200 raw-order datasets with K_obs >= 1, without
filtering on p-values or chain behavior. `dataset_id` is local to delta and
`raw_id` records the attempt. New preparation has a 10,000-attempt limit per
signal and reports delta, retained count and attempts if it cannot fill a bank.
Old frozen configurations without `raw_attempt_limit` use that same default
without rewriting their files. Recorded attempts are 881, 435, 253, 212 and 200.

V contains successive right-minus-left segment-mean contrasts. Write `T=V'Y`,
`Sigma=V'V`, `R=Y-V*Sigma^(-1)T`, `a=V'b` and `g=a'Sigma^(-1)a`.
The upper-tail statistic is `a'Sigma^(-1)t/sqrt(g)`, with fixed direction and
residual and null reference `N(0,Sigma)`. Delta generates alternatives; it does
not shift this null center. Structural-zero directions return p=1 before either
sampling branch, without fabricated tail counts or diagnostics. None occurs in
the supplied bank.

Data splitting selects on Y+W and infers from `Y_fiss=Y-100W`. These independent
Gaussian views have covariances 1.01 I and 101 I. Its p-value is the normal upper
tail of `a'Sigma^(-1)T_fiss/sqrt(101*g)` using the same selected set and direction,
with no reselection or random draws. The factors 100 and 101 follow from SD 0.1.

## Dependencies, files and commands

- `functions_nonnull_cp.R`: configuration, selector/event, target/statistic,
  bank generation, paired results and validation.
- `samplers_nonnull_cp.R`: directional PMMH and zero-signal reuse.
- `prepare_nonnull_cp.R`: bank preparation without PMMH.
- `inference_nonnull_cp.R`: paired inference on a frozen or supplied bank.
- `plot_inference_nonnull_cp.R`: PDF via the root `plot_style.R`.
- `inference_nonnull_cp.rds` / `.pdf`: frozen results and figure.

R 4.4.3 is the verified environment; inference diagnostics require posterior
1.7.0. Preparation uses base R and parallel without a diagnostic version lock.
Plotting requires ggplot2 4.0.3, patchwork 1.3.0 and Cairo. Positive-signal
sampling and the offline diagnostics use R 4.4.3 and posterior 1.7.0;
the original zero-signal sampling environment is not fully recorded.
Packages are not installed automatically.
From this directory:

```sh
Rscript plot_inference_nonnull_cp.R
Rscript inference_nonnull_cp.R rerun_inference_nonnull_cp.rds
Rscript plot_inference_nonnull_cp.R rerun_inference_nonnull_cp.rds
Rscript prepare_nonnull_cp.R bank_nonnull_cp.rds
Rscript inference_nonnull_cp.R new_inference_nonnull_cp.rds bank_nonnull_cp.rds
```

Default inputs and dependencies are script-relative. Explicit relative paths use
the caller's directory. Preparation/inference refuse existing outputs and start
full computations. A supplied bank controls configuration, data and task streams.
Source an entry point to expose its functions without executing the main program.

## Kernels, provenance and results

PMMH uses 200 warm-up transitions, zero stabilization, 5,000 production draws,
one auxiliary sample, rho=0.8, initial h=0.8 and bounds [0.05,5]. Warm-up adapts
`log(h) += 0.5/(i+10)^0.6*(acceptance probability - 0.20)`; production freezes h.
Rejection holds both states. At positive signals, the implementation is

```r
L <- t(chol(Sigma))
T_proposed <- T + (h / sqrt(K)) * drop(L %*% xi)
```

Here `tcrossprod(L) = Sigma` and K is the observed target dimension, so the
proposal covariance is `(h^2/K)*Sigma`. This symmetric random walk uses the
Gaussian MH ratio and an unconditional uniform draw. Positive-signal chains use
5,200 selector calls with their observed initializer taken as given.

Delta zero instead reuses the quadratic-null runner with a directional score:
it uses the lower factor `t(chol(Sigma))`, covariance `(h^2/K)*Sigma`, and draws
a uniform only for positive acceptance probability. Its recorded 5,201 calls
include initializer verification. Shared n, randomization SD, segment bounds,
penalty, tie tolerance, test level and PMMH settings are explicitly passed to the
null configuration. Supported bounds remain SD=0.1, minimum segment=10,
maximum changepoints=4 and penalty=1.25; incompatible settings fail before work.
The fixed positive-signal adaptation schedule cannot be changed through optional
null-kernel fields.

The 200 delta-zero datasets are the global-null bank (last raw ID 881), with
original generation states and the Proposed streams at positions `2*id-1` in
the quadratic-null assignment. The signal template remains nonzero. New dataset
and task `source_batch=0` means global-null origin. In the supplied task table,
delta-zero `source_batch=1/2` only groups IDs; its original streams are retained.
Positive signals retain complementary bank batches 1/2 (IDs 1--100 / 101--200).
Their Proposed task streams use master seed 20261030 and L'Ecuyer-CMRG with
Inversion/Rejection, allocating successive `parallel::nextRNGStream()` states
in sorted `(delta, dataset_id)` order. The saved `rng_state` is the task's actual
initial state; `seed` identifies its stream master, not an independent task seed.
Bank generation seeds are unchanged. Preparation and inference use these same
identity-based assignments, independent of worker scheduling.

P-values are inclusive tail counts / 5,000 with no +1 correction; reject at
p <= 0.05. Power error bars are +/- sqrt(power*(1-power)/200), untruncated
numerically, on a [0,1] plot. Delta-zero rejection rates are 0.025 for Proposed
and 0.040 for Data splitting. No Data-splitting chain or MCSE is invented.
At positive signals, Proposed rejects 26, 64, 115 and 153 of 200 datasets;
Data splitting rejects 8, 20, 19 and 27, respectively.

The root `inference_trajectories.rds` provides ordered production statistics,
including shared keys for quadratic-null and delta-zero tests, rather than full
target/auxiliary chains. `Rscript ../../../reproduce.R` computes each actual
upper-tail indicator's mean ESS/MCSE using posterior 1.7.0. Constant, too-short or
unestimable series retain NA with objective reasons; nonfinite input errors.
Finite ESS and MCSE are pooled independently without dropping p-values. A split
trajectory is not independent multi-chain initialization. Saved prefix and
finite-tail summaries describe precision limitations, not convergence evidence.
Full bank reconstruction and chain reruns have not been runtime-tested, and
bitwise recovery also depends on the R/numerical environment.
