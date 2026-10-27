options(stringsAsFactors = FALSE)

source("functions_nonnull_lasso.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript prepare_nonnull_lasso.R bank_nonnull_lasso.rds")
output_path <- args[[1L]]
if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

bank <- prepare_bank_nonnull_lasso(default_config_nonnull_lasso())
saveRDS(bank, output_path, version = 3)
message("Saved data bank to ", normalizePath(output_path, winslash = "/"))
