options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")
if (!capabilities("cairo")) stop("Cairo PDF support is required for the paper figure.")
source("functions_nonnull_cp.R", local = TRUE)
source("../../../plot_style.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript plot_inference_nonnull_cp.R [result.rds]")
input_path <- if (length(args)) args[[1L]] else "inference_nonnull_cp.rds"
if (!file.exists(input_path)) stop("Results file not found: ", input_path)
result <- readRDS(input_path)
validate_result_nonnull_cp(result)

power <- result$summary$power
display_ids <- c(full_selector = "model_only", full_path = "full_conditioning", data_fission = "data_fission")
power$method_id <- unname(display_ids[power$method_id])
figure <- make_power_figure(power, c("Full-selector", "Full-path", "Data fission"))
output_path <- paste0(tools::file_path_sans_ext(input_path), ".pdf")
ggsave(
  output_path, figure, width = 6.4, height = 4.4, units = "in",
  device = cairo_pdf, bg = "white"
)
