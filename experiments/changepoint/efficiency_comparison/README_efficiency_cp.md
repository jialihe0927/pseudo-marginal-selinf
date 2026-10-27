# Changepoint sampling efficiency comparison

This experiment compares rejection sampling, Crude PMMH (independent auxiliary
refreshment), and Correlated PMMH for a selective Gaussian changepoint target.

## Files and dependencies

- `functions_cp.R`: sequence reconstruction, bounded binary segmentation, and the selection event.
- `samplers_cp.R`: rejection and the two PMMH kernels.
- `sampling_efficiency_cp.R`: formal simulation entry point.
- `table_efficiency_cp.R`: reads the RDS and prints the frozen paper table.
- `efficiency_cp.rds`: frozen benchmark and paper estimates.

Packages are not installed automatically. The original formal-run versions of
R and `posterior` were not recorded in the available run artifacts and are
therefore unknown. This release was checked with R 4.4.3 and posterior 1.7.0.
These current versions are validation information, not a claim about the
historical environment.

## Frozen setting

The 200-point benchmark was generated with signal parameter 1.50 and conditions
on the unordered final changepoint set `{71, 99, 140}`. The selector follows
four bounded binary-segmentation splits,
allows K = 0,...,4, uses minimum segment length 10 and information-criterion
penalty coefficient 1.25, and selects K before testing equality of final sets.
The auxiliary scale is 0.10 and the rare-center displacement is 4 marginal
standard deviations in the first target coordinate.

Both PMMH chains use one auxiliary draw and propose
`T' = T + h_prop / sqrt(d) * G * xi`, where `G G^T = Sigma_T` and
`xi ~ N(0, I_d)`. Thus `Sigma_T` fixes the covariance shape and the actual
proposal covariance is `(h_prop^2 / d) * Sigma_T`. The scalar `h_prop` starts at
0.8, adapts only during the 500 warm-up transitions toward acceptance 0.20, and
is fixed for the 19,500 production transitions. Correlated PMMH uses auxiliary
correlation 0.8; crude PMMH uses independent refreshment (`rho = 0`). Each
transition costs one selector
evaluation. Efficiency is the minimum coordinate bulk ESS divided by all
20,000 calls, times 1,000. Rejection uses exactly 3,000,000 candidate attempts
in 70 fixed RNG streams (69 of size 42,857 and one of size 42,867); the original
run used seven workers. This is not 3,000,000 accepted samples.

The bundled summaries are the paper results. The exact fixed benchmark is
included; new simulation has ordinary Monte Carlo variation and is not promised
to restore the historical random stream bit for bit.

The paper table reports `(h_final, efficiency, gain, acceptance, zero rate)` as
`(--, 0.026, 1.0, --, --)`, `(0.921, 11.793, 453.6, 0.214, 0.466)`, and
`(1.206, 15.539, 597.6, 0.184, 0.432)` for rejection, crude PMMH, and correlated
PMMH, respectively.

Here `h_final` is the proposal scale at the end of warm-up, held fixed during
production. `Acc.` is the actual production acceptance fraction. `Zero rate` is
the fraction of production proposals with a zero Bernoulli selection estimate;
it is not the holding fraction, since a positive estimate can still be rejected
by MH. Each formal PM/CPM budget is 20,000 selector evaluations, including
warm-up but excluding the extra initial selection verification.

Run from this directory:

```sh
Rscript table_efficiency_cp.R
Rscript sampling_efficiency_cp.R pmmh rerun_pmmh_cp.rds
Rscript sampling_efficiency_cp.R rejection rerun_rejection_cp.rds
Rscript sampling_efficiency_cp.R all rerun_cp.rds
Rscript table_efficiency_cp.R rerun_cp.rds
```

The entry point refuses to overwrite an existing result. Only `all` produces a
complete three-method object that can reproduce the paper table. Original
wall-clock time is unknown. Before simulation, the entry point checks the
frozen transition accounting, the 3,000,000-attempt rejection budget, the
Bernoulli auxiliary-sample requirement, proposal-scale bounds, and correlation.
The full simulation was not executed during release validation.

