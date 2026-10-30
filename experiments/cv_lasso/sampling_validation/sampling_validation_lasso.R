script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(script_file)

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  options(stringsAsFactors = FALSE)

  source(file.path(script_directory, "../../../diagnostics.R"), local = TRUE)
  require_diagnostic_version()
  source(file.path(script_directory, "functions_validation_lasso.R"), local = TRUE)

  # Offline analysis is the default; sampling requires an explicit --generate.
  generate <- length(arguments) == 2L && identical(arguments[[1L]], "--generate")
  if (length(arguments) && startsWith(arguments[[1L]], "--") && !generate) {
    stop("Explicit sampling requires --generate output.rds.")
  }
  if (!generate && !length(arguments) %in% c(1L, 2L)) {
    stop("Usage: Rscript sampling_validation_lasso.R output.rds [input.rds]; or --generate output.rds")
  }
  output_path <- if (generate) arguments[[2L]] else arguments[[1L]]
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)
  input_path <- if (!generate && length(arguments) == 2L) {
    arguments[[2L]]
  } else file.path(script_directory, "sampling_validation_lasso.rds")
  if (!file.exists(input_path)) stop("Validation input not found: ", input_path)
  # Validate complete production trajectories before analyzing or regenerating.
  object <- readRDS(input_path)
  validate_samples_lasso(object)
  result <- if (generate) {
    regenerate_validation_lasso(
      object,
      file.path(script_directory, "..", "efficiency_comparison")
    )
  } else analyze_validation_lasso(object)

  # Save to the protected output path and report numerical summaries.
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output_path, version = 3)
  cat("Mean diagnostics: posterior::ess_basic(split = TRUE), posterior::mcse_mean;",
      "R", result$config$software$diagnostics[["R"]],
      "posterior", result$config$software$diagnostics[["posterior"]], "\n")
  cat("Undefined continuous ESS/MCSE:",
      sum(!is.finite(result$diagnostics$ess_mcse$ess_mean)), "\n")
  cat("Undefined MCMC tail indicator ESS/MCSE:",
      sum(result$diagnostics$tails$method != "Rejection" &
            !is.finite(result$diagnostics$tails$ess_indicator)), "\n")
  cat("Retained samples by method:\n")
  print(vapply(result$samples, nrow, integer(1L)))
  print(result$diagnostics$mixing, row.names = FALSE)
  message("Saved ", normalizePath(output_path, winslash = "/"))
}

if (sys.nframe() == 0L) {
  main()
}
