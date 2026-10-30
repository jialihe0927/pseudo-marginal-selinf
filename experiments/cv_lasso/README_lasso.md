# CV-Lasso experiments

The CV-Lasso application is organized into four focused experiments:

- [Sampling efficiency comparison](efficiency_comparison/README_efficiency_lasso.md): organized and reproducible from the supplied result object.
- [Global-null calibration](inference_null/README_inference_null_lasso.md): organized code, frozen paired p-values, and the paper ECDF figure.
- [Non-null inference and power](inference_nonnull/README_inference_nonnull_lasso.md): paired directional tests, frozen power results, bank preparation and the paper figure.
- [Sampling accuracy validation](sampling_validation/README_sampling_validation_lasso.md): full saved samples, offline diagnostics and the paper coordinate ECDFs.

Each experiment separates reusable functions from preparation, computation and
output entry points. Sourcing an entry defines `main()` without running it.
Default inputs are located beside the script; explicit relative paths are
resolved from the caller's working directory.

The inference implementations support unit response noise (`response_sd = 1`).
The directional delta-zero branch checks shared design, selector and eligibility
settings against the global-null implementation before computation.

The recorded lambda grid uses glmnet's Gaussian objective
`sum((y - X beta)^2) / (2 * n_fit) + lambda * sum(abs(beta))`.
For a paper objective with unnormalized half squared error, the corresponding
penalty is `n_fit * lambda`; CV uses each training fold's sample size.
The supplied lambda values are retained at their recorded precision.

