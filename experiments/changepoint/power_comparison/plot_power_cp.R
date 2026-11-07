script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(script_file)

suppressPackageStartupMessages(library(ggplot2))
source(file.path(script_directory, "../../../plot_style.R"), local = TRUE)

make_directional_figure <- function(result) {
  make_power_figure(result$summary$power)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Usage: Rscript plot_power_cp.R [result.rds]")
  input_path <- if (length(args)) args[[1L]] else file.path(script_directory, "power_cp.rds")
  result <- readRDS(input_path)

  figure <- make_directional_figure(result)
  output_path <- paste0(tools::file_path_sans_ext(input_path), ".pdf")
  ggsave(
    output_path, figure, width = 4.62, height = 3.6, units = "in",
    device = cairo_pdf, bg = "white"
  )
}

if (sys.nframe() == 0L) main()
