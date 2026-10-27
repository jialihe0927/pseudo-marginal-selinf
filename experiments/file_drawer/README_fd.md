# Gaussian file-drawer experiment

This directory contains the self-contained R implementation and paper results
for the Gaussian file-drawer experiment. It compares rejection sampling,
Independent Pseudo-Marginal (PM), and Correlated Pseudo-Marginal (CPM) sampling
by the number of observations needed to reach a fixed accuracy target.

## Files

- `samplers_fd.R`: selection event, exact selected-state initialization, and the
  Independent PM and CPM transitions.
- `analysis_fd.R`: selected-distribution calculations, fixed-accuracy estimation,
  and the analytic rejection calculation.
- `run_fd.R`: the complete simulation entry point and final experiment settings.
- `plot_fd.R`: reads one RDS result object and creates the two-panel paper figure.
- `results_fd.rds`: supplied paper estimates and the configuration needed
  to interpret them.
- `figure1_fd.pdf`: the paper figure.

Run all commands below from this directory.

## Environment

The code requires R and `ggplot2`; all numerical work uses base R. The release
was checked with R 4.4.3 and ggplot2 4.0.3. Packages are not installed
automatically.

## Reproduce the figure

```sh
Rscript plot_fd.R
```

This reads `results_fd.rds` and writes `figure1_fd.pdf`. A different result file
can be plotted with, for example, `Rscript plot_fd.R rerun_fd.rds`; its figure
is written beside it as `rerun_fd.pdf`.

## Rerun the experiment

```sh
Rscript run_fd.R rerun_fd.rds
```

The default output is `rerun_fd.rds`. The command stops rather than
overwriting an existing file. A full run is computationally expensive; its
original wall-clock time is unknown.

The supplied RDS contains the paper results. The historical random streams are
not available in this compact release, so `run_fd.R` uses a documented,
deterministic seed scheme and a rerun will have ordinary Monte Carlo variation;
bitwise reproduction is not expected.

## Experiment and metric

The displayed grid is `b = 2, ..., 8`. The fixed regime uses `a = 0.9`; the
adaptive regime uses `a_b = 1 - b^-2`. The proposal scale is
`2 / (b * sqrt(a))`. Independent PM uses `rho = 0`, while fixed-regime CPM uses
`rho = 1 - b^-2`. The paper figure includes rejection, Independent PM, and CPM
in the fixed regime, and rejection plus Independent PM in the adaptive regime.
Adaptive-regime CPM is implemented but is not part of the default run.

Each of 1,000 chains starts exactly from the selected distribution, with no
warm-up. The initial state is observation 1, so a run with `n` observations has
`n - 1` transitions. Each transition is charged one selector evaluation. The
observable is the indicator that `X` is below the exact selected-distribution
median. At checkpoints spaced approximately 5% apart, accuracy is measured as
four times the across-chain variance of the chain averages. The reported
fixed-accuracy time is the first of eight consecutive checkpoints at or below
0.10. An unresolved crossing remains `NA`; it is never replaced by the run cap.

The plotted point estimate is the raw stable-crossing start. Rejection values
are analytic rather than simulated. All plotted logarithms are natural
logarithms. The dashed and dotted reference curves are formula-based growth
guides with fixed vertical display offsets; those offsets are not theoretical
constants or fitted simulation estimates.
