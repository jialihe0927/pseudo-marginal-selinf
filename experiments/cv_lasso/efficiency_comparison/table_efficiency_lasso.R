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
  if (length(arguments) > 1L) stop("Usage: Rscript table_efficiency_lasso.R [efficiency_lasso.rds]")
  input <- if (length(arguments)) arguments[[1L]] else file.path(script_directory, "efficiency_lasso.rds")
  if (!file.exists(input)) stop("Results file not found: ", input)

  table <- efficiency_table(readRDS(input))
  print(table, row.names = FALSE, digits = 15)
}

if (sys.nframe() == 0L) {
  main()
}

