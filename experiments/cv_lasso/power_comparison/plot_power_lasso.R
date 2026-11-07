suppressPackageStartupMessages(library(ggplot2))

make_directional_figure <- function(result) {
  power <- result$summary$power
  make_power_figure(power)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("../../../plot_style.R", local = environment(make_directional_figure))
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")
  if (!capabilities("cairo")) stop("Cairo PDF support is required for the paper figure.")
  if (length(args) > 1L) stop("Usage: Rscript plot_power_lasso.R [result.rds]")
  input_path <- if (length(args)) args[[1L]] else "power_lasso.rds"
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
