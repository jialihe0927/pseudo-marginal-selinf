options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(ggplot2))
source("functions_null_cp.R", local = TRUE)

step_ribbon_data <- function(data) {
  data <- data[order(data$x), , drop = FALSE]
  if (nrow(data) < 2L) return(data)
  index <- rep(seq_len(nrow(data)), each = 2L)[-c(1L, 2L * nrow(data))]
  output <- data[index, , drop = FALSE]
  output$x <- c(rbind(data$x[-nrow(data)], data$x[-1L]))
  output
}

marker_data <- function(data, method_ids, marker_x) {
  rows <- list()
  for (method_id in method_ids) {
    method_data <- data[data$method_id == method_id, , drop = FALSE]
    method_data <- method_data[order(method_data$x), , drop = FALSE]
    for (x_value in marker_x) {
      index <- max(which(method_data$x <= x_value))
      rows[[length(rows) + 1L]] <- data.frame(
        method_id = method_id,
        x = x_value,
        ecdf = method_data$ecdf[[index]]
      )
    }
  }
  do.call(rbind, rows)
}

make_ecdf_data <- function(p_values, method_id) {
  n <- length(p_values)
  x <- sort(unique(c(0, p_values, 0.05, 0.20, 0.40, 0.60, 0.80, 1)))
  empirical_cdf <- vapply(x, function(value) mean(p_values <= value), numeric(1L))
  standard_error <- sqrt(empirical_cdf * (1 - empirical_cdf) / n)

  data.frame(
    method_id = method_id,
    x = x,
    ecdf = empirical_cdf,
    lower = pmax(0, empirical_cdf - standard_error),
    upper = pmin(1, empirical_cdf + standard_error),
    stringsAsFactors = FALSE
  )
}

make_inference_figure <- function(result) {
  required <- c("config", "dataset_results")
  if (!all(required %in% names(result))) stop("The result object is incomplete.")

  validate_result_null_cp(result)
  data <- result$dataset_results
  method_ids <- as.character(result$config$methods$method_id)

  plot_data <- do.call(
    rbind,
    lapply(method_ids, function(method_id) {
      make_ecdf_data(data$p_value[data$method_id == method_id], method_id)
    })
  )
  plot_data$method_id <- factor(plot_data$method_id, levels = method_ids)

  ribbon <- do.call(
    rbind,
    lapply(method_ids, function(method_id) {
      step_ribbon_data(plot_data[plot_data$method_id == method_id, , drop = FALSE])
    })
  )
  ribbon$method_id <- factor(ribbon$method_id, levels = method_ids)
  markers <- marker_data(plot_data, method_ids, c(0.20, 0.40, 0.60, 0.80))
  markers$method_id <- factor(markers$method_id, levels = method_ids)

  colors <- c(
    full_selector = "#4878A8",
    fixed_k = "#94B786",
    naive_gaussian = "#9A86B2"
  )
  shapes <- c(full_selector = 0, fixed_k = 2, naive_gaussian = 4)
  labels <- c(
    full_selector = "Full-selector",
    fixed_k = "Fixed-K",
    naive_gaussian = "Naive"
  )

  ggplot(plot_data, aes(x = x, y = ecdf, group = method_id)) +
    geom_ribbon(
      data = ribbon,
      aes(ymin = lower, ymax = upper, fill = method_id),
      alpha = 0.30,
      color = NA,
      show.legend = FALSE
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "#222222",
      linetype = "dashed",
      linewidth = 0.65
    ) +
    geom_vline(
      xintercept = 0.05,
      color = "black",
      linetype = "dashed",
      linewidth = 0.65
    ) +
    geom_step(
      aes(color = method_id),
      direction = "hv",
      linewidth = 0.95,
      lineend = "butt"
    ) +
    geom_point(
      data = markers,
      aes(color = method_id, shape = method_id),
      size = 2.5,
      stroke = 0.70
    ) +
    scale_color_manual(
      name = NULL,
      values = colors,
      breaks = method_ids,
      labels = labels,
      drop = FALSE
    ) +
    scale_fill_manual(values = colors, breaks = method_ids, guide = "none") +
    scale_shape_manual(
      name = NULL,
      values = shapes,
      breaks = method_ids,
      labels = labels,
      drop = FALSE
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(x = "P-value", y = "Empirical CDF") +
    theme_classic(base_size = 13, base_family = "sans") +
    theme(
      axis.title = element_text(size = 13.5, color = "#222222"),
      axis.text = element_text(size = 11.5, color = "#222222"),
      axis.line = element_line(linewidth = 0.7, color = "#3F3F3F"),
      axis.ticks = element_line(linewidth = 0.65, color = "#3F3F3F"),
      axis.ticks.length = grid::unit(4.5, "pt"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.ontop = FALSE,
      legend.position = c(0.975, 0.045),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.45),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.width = grid::unit(28, "pt"),
      legend.key.height = grid::unit(18, "pt"),
      legend.margin = margin(4, 5, 4, 5, unit = "pt"),
      legend.text = element_text(size = 12, color = "#222222"),
      plot.margin = margin(10, 12, 10, 12, unit = "pt")
    )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Usage: Rscript plot_inference_null_cp.R [results.rds]")
}
input_path <- if (length(args)) args[[1L]] else "inference_null_cp.rds"
if (!file.exists(input_path)) stop("Results file not found: ", input_path)

output_path <- if (identical(input_path, "inference_null_cp.rds")) {
  "inference_null_cp.pdf"
} else {
  paste0(tools::file_path_sans_ext(input_path), ".pdf")
}

figure <- make_inference_figure(readRDS(input_path))
ggsave(
  output_path,
  figure,
  width = 6.4,
  height = 4.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
