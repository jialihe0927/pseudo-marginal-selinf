main <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("functions_calibration_lasso.R", local = TRUE)

  if (length(args) != 1L) stop("Usage: Rscript prepare_calibration_lasso.R bank_calibration_lasso.rds")
  output_path <- args[[1L]]
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

  bank <- prepare_bank_null_lasso(default_config_null_lasso())
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bank, output_path, version = 3)
  message("Saved data bank to ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
