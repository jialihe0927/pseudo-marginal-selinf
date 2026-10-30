script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(script_file)

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options(stringsAsFactors = FALSE)

  source(file.path(script_directory, "functions_null_lasso.R"), local = TRUE)

  if (length(args) != 1L) stop("Usage: Rscript prepare_null_lasso.R bank_null_lasso.rds")
  output_path <- args[[1L]]
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

  bank <- prepare_bank_null_lasso(default_config_null_lasso())
  validate_bank_null_lasso(bank)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bank, output_path, version = 3)
  message("Saved data bank to ", normalizePath(output_path, winslash = "/"))
}

if (sys.nframe() == 0L) {
  main()
}
