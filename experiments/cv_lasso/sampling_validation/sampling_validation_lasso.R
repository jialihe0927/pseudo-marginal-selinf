options(stringsAsFactors = FALSE)

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_argument)) stop("Run this entry point with Rscript.")
script_directory <- dirname(sub("^--file=", "", script_argument[[1L]]))
source(file.path(script_directory, "functions_validation_lasso.R"), local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
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
object <- readRDS(input_path)
validate_samples_lasso(object)
result <- if (generate) {
  regenerate_validation_lasso(object, file.path(script_directory, "..", "efficiency_comparison"))
} else analyze_validation_lasso(object)

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
