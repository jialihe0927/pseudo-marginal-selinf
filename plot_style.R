# Shared geometry and typography for the inference panels.
inference_tick_labels <- function(x) {
  # Plot labels must not inherit full-precision numerical export options.
  formatC(x, format = "f", digits = 1, decimal.mark = ".")
}

plot_style <- list(
  base_size = 13, axis_title_size = 10, axis_text_size = 8.5,
  legend_text_size = 11.5, bar_border_linewidth = 0.55, errorbar_linewidth = 0.85
)

theme_inference_panel <- function() {
  theme_classic(base_size = plot_style$base_size, base_family = "sans") +
    theme(
      axis.title = element_text(size = plot_style$axis_title_size, color = "#222222"),
      axis.text = element_text(size = plot_style$axis_text_size, color = "#222222"),
      axis.line = element_line(linewidth = 0.7, color = "#3F3F3F"),
      axis.ticks = element_line(linewidth = 0.65, color = "#3F3F3F"),
      axis.ticks.length = grid::unit(4.5, "pt"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = plot_style$legend_text_size, color = "#222222"),
      legend.key.width = grid::unit(35, "pt"),
      legend.key.height = grid::unit(14, "pt"),
      legend.spacing.x = grid::unit(8, "pt"),
      legend.margin = margin(0, 0, 1, 0, unit = "pt"),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.ontop = FALSE,
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        size = 10, face = "plain", hjust = 0.5, vjust = 0.5, color = "#222222",
        margin = margin(b = 5, unit = "pt")
      ),
      plot.margin = margin(5, 7, 6, 7, unit = "pt")
    )
}


bar_geometry <- function(data, method_ids) {
  span <- if (length(method_ids) == 2L) 0.13 else 0.25
  offsets <- setNames(seq(-span, span, length.out = length(method_ids)), method_ids)
  data$group_x <- match(data$delta, sort(unique(data$delta)))
  data$x <- data$group_x + offsets[data$method_id]
  data$xmin <- data$x - 0.105
  data$xmax <- data$x + 0.105
  data
}

diagonal_hatch <- function(row, spacing = 0.18) {
  constants <- seq(-1, 1, by = spacing)
  pieces <- lapply(constants, function(cc) {
    u0 <- max(0, -cc)
    u1 <- min(1, 1 - cc)
    if (u0 >= u1) return(NULL)
    data.frame(
      x = row$xmin + u0 * (row$xmax - row$xmin),
      xend = row$xmin + u1 * (row$xmax - row$xmin),
      y = row$power * (u0 + cc),
      yend = row$power * (u1 + cc)
    )
  })
  do.call(rbind, pieces)
}

horizontal_hatch <- function(row) {
  vv <- seq(0.16, 0.88, by = 0.18)
  data.frame(
    x = row$xmin + 0.07 * (row$xmax - row$xmin),
    xend = row$xmax - 0.07 * (row$xmax - row$xmin),
    y = row$power * vv,
    yend = row$power * vv
  )
}

make_bar_hatches <- function(bars) {
  out <- list()
  for (i in seq_len(nrow(bars))) {
    row <- bars[i, , drop = FALSE]
    seg <- if (row$method_id == "full_conditioning") {
      diagonal_hatch(row)
    } else if (row$method_id == "data_fission") {
      horizontal_hatch(row)
    } else NULL
    if (!is.null(seg)) {
      seg$method_id <- row$method_id
      out[[length(out) + 1L]] <- seg
    }
  }
  do.call(rbind, out)
}

legend_hatches <- function(method_id, xmin, xmax, ymin, ymax) {
  row <- data.frame(xmin = xmin, xmax = xmax, power = ymax - ymin, method_id = method_id)
  if (method_id == "full_conditioning") {
    seg <- diagonal_hatch(row, spacing = 0.25)
  } else if (method_id == "data_fission") {
    seg <- horizontal_hatch(row)
  } else return(NULL)
  seg$y <- seg$y + ymin
  seg$yend <- seg$yend + ymin
  seg
}


make_power_figure <- function(power, legend_labels) {
  method_ids <- c("model_only", "full_conditioning", "data_fission")
  method_ids <- method_ids[method_ids %in% power$method_id]
  delta_grid <- sort(unique(power$delta))
  fills <- c(
    model_only = "#4878A8",
    full_conditioning = "#D8AA69",
    data_fission = "#BC6C65"
  )
  bars <- bar_geometry(power, method_ids)
  hatches <- make_bar_hatches(bars)
  edge <- "#4A4A4A"

  main <- ggplot() +
    geom_hline(yintercept = 0.05, color = "#222222", linetype = "dashed", linewidth = 0.65) +
    geom_rect(
      data = bars,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = power, fill = method_id),
      color = edge, linewidth = plot_style$bar_border_linewidth
    ) +
    geom_segment(
      data = hatches,
      aes(x = x, xend = xend, y = y, yend = yend),
      color = edge, linewidth = 0.38, lineend = "butt"
    ) +
    geom_errorbar(
      data = bars,
      aes(x = x, ymin = lower, ymax = upper),
      width = 0.075, color = edge, linewidth = plot_style$errorbar_linewidth
    ) +
    scale_fill_manual(values = fills, guide = "none") +
    scale_x_continuous(
      breaks = seq_along(delta_grid), labels = inference_tick_labels(delta_grid),
      expand = expansion(add = c(0.45, 0.45))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2), labels = inference_tick_labels, limits = c(0, 1),
      expand = expansion(mult = 0)
    ) +
    labs(x = expression(paste("Signal strength ", delta)), y = "Empirical rejection rate",
         title = "(b) Null rejection and power") +
    theme_inference_panel() +
    theme(legend.position = "none")

  two_methods <- length(method_ids) == 2L
  legend_y <- if (two_methods) c(0.7025, 0.2975) else c(0.77, 0.50, 0.23)
  half_height <- if (two_methods) 0.105 else 0.070
  legend_frame <- data.frame(xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  legend_boxes <- data.frame(
    method_id = method_ids,
    xmin = 0.10, xmax = 0.29,
    ymin = legend_y - half_height, ymax = legend_y + half_height
  )
  legend_text <- data.frame(x = 0.37, y = legend_y, label = legend_labels)
  legend_hatch_list <- lapply(seq_len(nrow(legend_boxes)), function(i) {
    legend_hatches(
      legend_boxes$method_id[[i]], legend_boxes$xmin[[i]], legend_boxes$xmax[[i]],
      legend_boxes$ymin[[i]], legend_boxes$ymax[[i]]
    )
  })
  legend_hatch <- do.call(rbind, Filter(Negate(is.null), legend_hatch_list))

  legend_plot <- ggplot() +
    geom_rect(
      data = legend_frame,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "white", color = "black", linewidth = 0.45
    ) +
    geom_rect(
      data = legend_boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = method_id),
      color = edge, linewidth = plot_style$bar_border_linewidth
    ) +
    geom_segment(
      data = legend_hatch,
      aes(x = x, xend = xend, y = y, yend = yend),
      color = edge, linewidth = 0.38, lineend = "butt"
    ) +
    geom_text(
      data = legend_text,
      aes(x = x, y = y, label = label), hjust = 0, vjust = 0.5,
      family = "sans", color = "#222222", size = 8.5 / ggplot2::.pt
    ) +
    scale_fill_manual(values = fills, guide = "none") +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork package is required")
  p <- main + patchwork::inset_element(
    legend_plot, left = 0.035, bottom = if (two_methods) 0.785 else 0.700,
    right = 0.400, top = 0.955,
    align_to = "panel", clip = FALSE
  )


  p
}

theme_null_legend <- function() {
  theme(
    legend.position = c(0.970, 0.060), legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.45),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = grid::unit(20, "pt"), legend.key.height = grid::unit(15, "pt"),
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = margin(3, 4, 3, 4, unit = "pt"),
    legend.text = element_text(size = 8.5, color = "#222222")
  )
}
