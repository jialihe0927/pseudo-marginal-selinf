script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))

main <- function() {
  options(stringsAsFactors = FALSE)

  source(file.path(script_directory, "../../../diagnostics.R"), local = TRUE)
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) > 1L) stop("Usage: Rscript table_efficiency_cp.R [efficiency_cp.rds]")
  input <- if (length(arguments)) arguments[[1L]] else file.path(script_directory, "efficiency_cp.rds")
  if (!file.exists(input)) stop("Results file not found: ", input)

  table <- efficiency_table(readRDS(input))
  print(table, row.names = FALSE, digits = 15)
}

if (sys.nframe() == 0L) main()
