options(stringsAsFactors = FALSE)

source("functions_null_cp.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript prepare_null_cp.R bank_null_cp.rds")
output_path <- args[[1L]]
if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

bank <- prepare_bank_null_cp(default_config_null_cp())
validate_bank_null_cp(bank)
saveRDS(bank, output_path, version = 3)
message("Saved data bank to ", normalizePath(output_path, winslash = "/"))
