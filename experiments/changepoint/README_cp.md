# Changepoint experiments

Four experiments use randomized bounded binary segmentation:

- [Sampling efficiency](efficiency_comparison/README_efficiency_cp.md): PM, CPM and rejection for a fixed rare-event Gaussian target.
- [Global-null calibration](inference_null/README_inference_null_cp.md): paired selective and naive p-values.
- [Directional inference](inference_nonnull/README_inference_nonnull_cp.md): Proposed and Data splitting at five signal strengths, including zero.
- [Sampling validation](sampling_validation/README_sampling_validation_cp.md): saved production samples, offline diagnostics and coordinate ECDFs.

R 4.4.3 is the verified environment. New mean ESS/MCSE diagnostics require
posterior 1.7.0; plotting uses ggplot2 4.0.3, patchwork 1.3.0 where needed, and
Cairo PDF support. Bank preparation and pure plotting do not require posterior.
No package is installed automatically.

From the repository root, `Rscript reproduce.R` recomputes tables and inference
diagnostics offline. Experiment entry points resolve default inputs and function
sources relative to their scripts; explicit relative paths use the caller's
working directory. Sourcing an entry point does not execute its main program.
See the linked commands before starting bank generation or sampling. Complete
bank reconstruction and chain reruns have not been runtime-tested.
