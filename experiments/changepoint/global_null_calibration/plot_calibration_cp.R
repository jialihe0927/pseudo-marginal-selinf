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

# Pointwise across-dataset ECDF +/- one SE, clipped to [0, 1].
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
  shapes <- c(full_selector = 1, fixed_k = 2, naive_gaussian = 4)
  labels <- unname(c(full_selector = "Proposed", fixed_k = "Fixed-K",
                     naive_gaussian = "Naive")[method_ids])

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
    scale_x_continuous(
      breaks = seq(0, 1, 0.2), labels = inference_tick_labels,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2), labels = inference_tick_labels,
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(x = "P-value", y = "Empirical CDF") +
    labs(title = expression(paste("(a) Global-null ", italic(p), "-values"))) +
    theme_inference_panel() +
    theme_null_legend()
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) {
    stop("Usage: Rscript plot_calibration_cp.R [results.rds]")
  }
  input_path <- if (length(args)) args[[1L]] else file.path(script_directory, "calibration_cp.rds")

  output_path <- paste0(tools::file_path_sans_ext(input_path), ".pdf")

  figure <- make_inference_figure(readRDS(input_path))
  ggsave(
    output_path,
    figure,
    width = 3.78,
    height = 3.6,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}

if (sys.nframe() == 0L) main()
