suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")
if (!capabilities("cairo")) stop("Cairo PDF support is required for the paper figure.")
panel_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
panel_directory <- dirname(panel_file)
source(file.path(panel_directory, "../inference_null/functions_null_lasso.R"), local = TRUE)
source(file.path(panel_directory, "functions_nonnull_lasso.R"), local = TRUE)
source(file.path(panel_directory, "../../../plot_style.R"), local = TRUE)

make_directional_figure <- function(result) {
  validate_result_nonnull_lasso(result)
  power <- result$summary$power
  make_power_figure(power, c("Proposed", "Extra cond.", "Data splitting"))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) > 1L) stop("Usage: Rscript plot_inference_nonnull_lasso.R [result.rds]")
  input_path <- if (length(args)) args[[1L]] else file.path(panel_directory, "inference_nonnull_lasso.rds")
  if (!file.exists(input_path)) stop("Results file not found: ", input_path)
  result <- readRDS(input_path)

  figure <- make_directional_figure(result)
  output_path <- paste0(tools::file_path_sans_ext(input_path), ".pdf")
  ggsave(
    output_path, figure, width = 4.62, height = 3.6, units = "in",
    device = cairo_pdf, bg = "white"
  )
}

if (sys.nframe() == 0L) {
  main()
}
