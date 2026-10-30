# Gaussian file-drawer experiment

This experiment compares rejection, Independent Pseudo-Marginal (PM), and
Correlated Pseudo-Marginal (CPM) sampling by the per-chain observation count
needed to reach a fixed accuracy target. The figure abbreviates Independent PM
as PM.

## Files

- `samplers_fd.R`: selection event, exact selected-state initialization, and
  Independent PM and CPM transitions.
- `analysis_fd.R`: selected-distribution calculations, stable crossing, and
  the analytic rejection calculation.
- `run_fd.R`: experiment settings, simulation, and result saving.
- `plot_fd.R`: the two-panel paper figure, without simulation.
- `results_fd.rds`: 35 paper estimates and their experiment configuration.
- `sampling_complexity_fd.pdf`: the paper figure.

## Dependencies

Numerical calculations use base R; plotting requires `ggplot2` and Cairo
support (`capabilities("cairo")`). Checked with R 4.4.3 and ggplot2 4.0.3.
Packages are not installed automatically. Sourcing the R files does not start
an experiment, read results, or write files.

## Plot the supplied results

Run these commands from this directory:

```sh
Rscript plot_fd.R
Rscript plot_fd.R rerun_fd.rds
```

The first command reads `results_fd.rds` and writes
`sampling_complexity_fd.pdf`. The second writes `rerun_fd.pdf` beside its input.
Only PDF is produced. When invoked by script path from another directory,
default inputs and dependencies are located beside the script; explicitly
supplied relative input and output paths remain relative to the caller.

## Run the complete experiment

```sh
Rscript run_fd.R rerun_fd.rds
```

This starts the full simulation. The default output is `rerun_fd.rds`; an
existing output is never overwritten. The original wall-clock time is unknown.

The supplied RDS contains the paper results. Historical random streams were
not retained; new simulations use L'Ecuyer-CMRG with seed
`9100000 + configuration_id` and do not promise bitwise reproduction.

## Settings and metrics

The grid is `b = 2, ..., 8`. The fixed regime uses `a = 0.9` and compares
Rejection, Independent PM, and CPM (21 displayed records). The adaptive regime
uses `a_b = 1 - b^-2` and compares Rejection and Independent PM (14 records).
The proposal scale is `2 / (b * sqrt(a))`. Independent PM uses `rho = 0`;
fixed-regime CPM uses `rho = 1 - b^-2`.

Each simulated configuration uses 1,000 chains initialized exactly from the
selected distribution, with no warm-up. The observable is
`I(X <= selected-distribution median)`. At checkpoints spaced approximately
5% apart, the variance ratio is four times the across-chain variance of the
chain averages. With `epsilon = 0.10`, the raw stable-crossing estimate is the
start of the earliest eight consecutive checkpoints whose variance ratios
are all at most epsilon. A missing crossing remains `NA`, not the run cap.

The result fields distinguish accuracy from simulation cost:

- `complexity_n`: per-chain observation count at the stable-crossing start,
  including the initial state; not the final stopping point.
- `observations_simulated`: per-chain observations actually retained through
  confirmation and the configured extra checkpoints.
- `transitions_per_chain`: `observations_simulated - 1`.
- `selector_evaluations`: `n_chains * (observations_simulated - 1)`, the total
  proposal event evaluations across chains under the experiment's billing
  convention. Initialization validation and checkpoint state-integrity
  checks are excluded; this is not every call to `selection_event()`.

The run configuration retains distinct observation caps and extra-checkpoint
counts. Fixed-regime Independent PM has no extra checkpoints after the
eight-checkpoint confirmation; fixed CPM and adaptive Independent PM have
eight extra checkpoints, subject to their caps. Neither this total cost nor
the confirmation endpoint replaces the plotted per-chain `complexity_n`.

New runs also save `checkpoint_summary`, with `regime`, `b`, `method`,
`observations`, and `variance_ratio`. These small summaries allow the crossing
start to be recomputed without storing chains. The frozen paper RDS lacks
these checkpoint variance sequences, so it alone cannot independently
recompute the historical simulated crossings.

Rejection values are computed from an analytic formula, not simulation.
The calculation retains stable small-probability expressions. At `b = 8`,
the crossing exceeds the range in which double precision represents every
integer; the numerical search does not resolve it one integer at a time.

All plotted logarithms are natural logarithms. Reference curves use the
theoretical growth forms. Their vertical positions are chosen using empirical
anchors and fixed display offsets; they are not estimates of theoretical
constants or fitted growth rates.
