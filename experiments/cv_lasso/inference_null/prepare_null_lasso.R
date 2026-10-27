options(stringsAsFactors = FALSE)

source("functions_null_lasso.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript prepare_null_lasso.R bank_null_lasso.rds")
output_path <- args[[1L]]
if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

bank <- prepare_bank_null_lasso(default_config_null_lasso())
validate_bank_null_lasso(bank)
saveRDS(bank, output_path, version = 3)
message("Saved data bank to ", normalizePath(output_path, winslash = "/"))
