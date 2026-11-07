main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  generate <- length(arguments) == 2L && identical(arguments[[1L]], "--generate")
  if (!length(arguments) %in% c(1L, 2L) ||
      (!generate && startsWith(arguments[[1L]], "--"))) {
    stop("Usage: Rscript sampling_validation_lasso.R output.rds [input.rds]; or --generate output.rds")
  }
  output_path <- if (generate) arguments[[2L]] else arguments[[1L]]
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)
  input_path <- if (!generate && length(arguments) == 2L) {
    arguments[[2L]]
  } else "sampling_validation_lasso.rds"

  for (package in if (generate) c("posterior", "glmnet") else "posterior") {
    if (!requireNamespace(package, quietly = TRUE)) stop("Package '", package, "' is required.")
  }
  source("../../../diagnostics.R", local = TRUE)
  source("functions_validation_lasso.R", local = TRUE)
  object <- readRDS(input_path)
  result <- if (generate) {
    regenerate_validation_lasso(object, "../efficiency_comparison")
  } else analyze_validation_lasso(object)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output_path, version = 3)
  cat("Undefined continuous ESS/MCSE:",
      sum(!is.finite(result$diagnostics$ess_mcse$ess_mean)), "\n")
  cat("Undefined MCMC tail indicator ESS/MCSE:",
      sum(result$diagnostics$tails$method != "Rejection" &
            !is.finite(result$diagnostics$tails$ess_indicator)), "\n")
  cat("Retained samples by method:\n")
  method_names <- c(Rejection = "Rejection", "Crude PMMH" = "PM", "Correlated PMMH" = "CPM")
  counts <- vapply(result$samples, nrow, integer(1L))
  names(counts) <- unname(method_names[names(counts)])
  print(counts)
  mixing <- result$diagnostics$mixing
  mixing$method_label <- unname(method_names[mixing$method])
  print(mixing[setdiff(names(mixing), "method")], row.names = FALSE)
  message("Saved ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
