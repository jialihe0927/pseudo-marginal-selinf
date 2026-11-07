main <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("../global_null_calibration/functions_calibration_lasso.R", local = TRUE)

  source("functions_power_lasso.R", local = TRUE)

  if (length(args) != 1L) stop("Usage: Rscript prepare_power_lasso.R bank_power_lasso.rds")
  output_path <- args[[1L]]
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

  bank <- prepare_bank_nonnull_lasso(default_config_nonnull_lasso())
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bank, output_path, version = 3)
  message("Saved data bank to ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
