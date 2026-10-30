script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))

suppressPackageStartupMessages(library(ggplot2))
source(file.path(script_directory, "functions_nonnull_cp.R"), local = TRUE)
source(file.path(script_directory, "../../../plot_style.R"), local = TRUE)

make_directional_figure <- function(result) {
  validate_result_nonnull_cp(result)
  power <- result$summary$power
  display_ids <- c(full_selector = "model_only", data_fission = "data_fission")
  power$method_id <- unname(display_ids[power$method_id])
  if (anyNA(power$method_id)) stop("Unsupported current-paper method.")
  make_power_figure(power, c("Proposed", "Data splitting"))
}

main <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")
  if (!capabilities("cairo")) stop("Cairo PDF support is required for the paper figure.")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript plot_inference_nonnull_cp.R [result.rds]")
  input_path <- if (length(args)) args[[1L]] else file.path(script_directory, "inference_nonnull_cp.rds")
  if (!file.exists(input_path)) stop("Results file not found: ", input_path)
  result <- readRDS(input_path)

  figure <- make_directional_figure(result)
  output_path <- paste0(tools::file_path_sans_ext(input_path), ".pdf")
  ggsave(
    output_path, figure, width = 4.62, height = 3.6, units = "in",
    device = cairo_pdf, bg = "white"
  )
}

if (sys.nframe() == 0L) main()
