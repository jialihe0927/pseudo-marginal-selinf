# Changepoint sampling efficiency

Compare rejection, PM (independent auxiliary refreshment) and CPM for one fixed
selective Gaussian target. The 200-point benchmark was generated at signal 1.50;
the adaptive selector returns the unordered set `{71, 99, 140}`. It constructs
four bounded segmentation splits with minimum segment length 10 and chooses
K = 0,...,4 by minimizing `RSS_K + 1.25 K log(n)` (tie tolerance 1e-12).
Randomization SD is 0.1.

The Gaussian center is displaced only in the first coordinate:

```r
mu_eff <- T_obs
mu_eff[1] <- T_obs[1] - 4 * sign(T_obs[1]) * sqrt(Sigma_T[1, 1])
```

Every other coordinate remains at `T_obs`. This is the fixed rare-event benchmark
center shared by PM, CPM and rejection; inference uses a zero-mean reference.
The covariance, residual and adaptive final-set event are also shared.

## Files, dependencies and commands

- `functions_cp.R`: configuration checks, selector/event and Gaussian target.
- `samplers_cp.R`: separate rejection and PMMH kernels.
- `sampling_efficiency_cp.R`: explicit simulation entry point.
- `table_efficiency_cp.R`: prints the three-method table without writing files.
- `efficiency_cp.rds`: frozen benchmark and estimates.

R 4.4.3 and posterior 1.7.0 were used for offline diagnostic verification.
Historical sampling versions and wall-clock times are unknown. No packages are
installed automatically. From this directory:

```sh
Rscript table_efficiency_cp.R
Rscript sampling_efficiency_cp.R pmmh rerun_pmmh_cp.rds
Rscript sampling_efficiency_cp.R rejection rerun_rejection_cp.rds
Rscript sampling_efficiency_cp.R all rerun_cp.rds
Rscript table_efficiency_cp.R rerun_cp.rds
```

Default inputs/dependencies are script-relative; explicit paths are relative to
the caller. Simulation refuses existing output paths. Only `all` produces a
complete three-method object. The root `reproduce.R` recomputes efficiency from
the frozen Sampling Validation trajectories without sampling or selector calls.
The full simulation has not been rerun to verify historical stream recovery.

## Kernels and efficiency

Both PMMH chains use one Bernoulli auxiliary draw and propose
`T' = T + h/sqrt(d) * G * xi`, where `G G' = Sigma_T` and `xi ~ N(0, I_d)`.
Proposal covariance is `(h^2/d) * Sigma_T`. The initial scale is 0.8; only the
500 warm-up transitions adapt, toward acceptance 0.20 with bounds [0.05, 5].
The remaining 19,500 production states include all holding states. PM uses
rho = 0 and CPM rho = 0.8. Chain uniforms, auxiliary draws and proposal draws
have separate RNG streams. Initializer verification is an extra selector call.

Efficiency is 1,000 times the minimum raw-coordinate mean ESS divided by all
20,000 transition selector calls, including warm-up. Mean ESS is
`posterior::ess_basic(matrix(as.numeric(x), ncol=1L), split=TRUE)`, using the
shared `diagnostics.R` implementation, with no cap or replacement estimator.
Gain is computed before rounding. Rejection uses 3,000,000 candidates in 70 fixed
streams (69 chunks of 42,857 and one of 42,867), historically on seven workers;
78 candidates were accepted.

The paper table `(h_final, efficiency, gain, acceptance, zero rate)` is
`(--, 0.026, 1.0, --, --)`, `(0.921, 11.948, 459.5, 0.214, 0.466)`, and
`(1.206, 15.451, 594.3, 0.184, 0.432)` for rejection, PM and CPM.
Acceptance and zero rate refer to production. Zero rate counts proposals with a
zero selection estimate; a nonzero estimate can still be rejected by MH.
