make_reference_guides <- function(estimates) {
  b_values <- sort(unique(estimates$b))
  fixed_pm_anchor <- estimates$complexity_n[
    estimates$regime == "fixed_a" &
      estimates$method == "Independent PM" & estimates$b == max(b_values)
  ]
  fixed_cpm_anchor <- estimates$complexity_n[
    estimates$regime == "fixed_a" &
      estimates$method == "CPM" & estimates$b == max(b_values)
  ]
  adaptive_pm_anchor <- median(estimates$complexity_n[
    estimates$regime == "adaptive_a" &
      estimates$method == "Independent PM"
  ])

  make_guide <- function(regime, method, values) {
    data.frame(
      regime = regime,
      b = b_values,
      method = method,
      log_value = values,
      stringsAsFactors = FALSE
    )
  }

  # Empirical anchors and fixed offsets position the theoretical growth forms;
  # their vertical positions are not theoretical constants or fitted rates.
  rbind(
    make_guide("fixed_a", "Rejection reference", b_values^2 / 2 + 4.2),
    make_guide("adaptive_a", "Rejection reference", b_values^2 / 2 + 4.2),
    make_guide(
      "fixed_a",
      "Independent PM reference",
      log(fixed_pm_anchor) + 0.05 * (b_values^2 - max(b_values)^2) + 0.65
    ),
    make_guide(
      "fixed_a",
      "CPM reference",
      log(fixed_cpm_anchor) + log(b_values^2 / max(b_values)^2) + 0.65
    ),
    make_guide(
      "adaptive_a",
      "Independent PM reference",
      rep(log(adaptive_pm_anchor) + 0.55, length(b_values))
    )
  )
}

make_figure <- function(result) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
  })

  estimates <- result$estimates

  method_colors <- c(
    "Rejection" = "#4d4d4d",
    "Independent PM" = "#bc6c65",
    "CPM" = "#4878a8",
    "Rejection reference" = "#4d4d4d",
    "Independent PM reference" = "#d8a29b",
    "CPM reference" = "#9ab7d3"
  )
  legend_order <- names(method_colors)
  shapes <- setNames(c(16, 17, 15, NA, NA, NA), legend_order)
  linetypes <- setNames(c(1, 1, 1, 4, 2, 3), legend_order)
  linewidths <- setNames(c(1.3, 1.3, 1.3, 0.75, 0.75, 0.75), legend_order)
  panel_levels <- c("fixed_a", "adaptive_a")

  # Preserve the paper's marker overplotting order where curves meet at small b.
  estimates <- estimates[order(match(estimates$regime, c("fixed_a", "adaptive_a")),
                               match(estimates$method, c("Rejection", "Independent PM", "CPM")),
                               estimates$b), ]
  empirical <- transform(
    estimates,
    panel = factor(regime, levels = panel_levels),
    b_squared = b^2,
    log_value = log(complexity_n),
    color = unname(method_colors[method]),
    shape = unname(shapes[method]),
    linetype = 1,
    linewidth = unname(linewidths[method])
  )
  guides <- make_reference_guides(estimates)
  guides <- transform(
    guides,
    panel = factor(regime, levels = panel_levels),
    b_squared = b^2,
    color = unname(method_colors[method]),
    linetype = unname(linetypes[method]),
    linewidth = unname(linewidths[method])
  )

  make_legend <- function(panel, keys, labels) {
    data.frame(
      panel = factor(panel, levels = panel_levels),
      key = keys,
      y = 38.5 - 2.30 * (seq_along(keys) - 1),
      label = labels,
      color = unname(method_colors[keys]),
      linetype = unname(linetypes[keys]),
      linewidth = unname(linewidths[keys]),
      shape = unname(shapes[keys]),
      stringsAsFactors = FALSE
    )
  }
  method_names <- c(Rejection = "Rejection", "Independent PM" = "PM", CPM = "CPM")
  fixed_legend <- make_legend(
    "fixed_a",
    legend_order,
    c(encodeString(unname(method_names[legend_order[1:3]]), quote = "'"),
      "slope~0.5", "slope~0.05", "log(b^2)")
  )
  adaptive_keys <- c(
    "Rejection", "Independent PM", "Rejection reference",
    "Independent PM reference"
  )
  adaptive_legend <- make_legend(
    "adaptive_a",
    adaptive_keys,
    c(encodeString(unname(method_names[adaptive_keys[1:2]]), quote = "'"),
      "slope~0.5", "slope~0")
  )
  custom_legend <- rbind(fixed_legend, adaptive_legend)
  legend_boxes <- data.frame(
    panel = factor(panel_levels, levels = panel_levels),
    xmin = 5.2,
    xmax = 22,
    ymin = c(min(fixed_legend$y), min(adaptive_legend$y)) - 1.25,
    ymax = 39.75
  )

  panel_labels <- c(
    "fixed_a" = "'(a)'~Fixed~italic(a)==0.9",
    "adaptive_a" = "'(b)'~Adaptive~italic(a)[b]==1~'-'~italic(b)^'-2'"
  )

  ggplot() +
    geom_line(
      data = guides,
      aes(
        x = b_squared, y = log_value,
        group = interaction(panel, method),
        color = color, linetype = linetype, linewidth = linewidth
      )
    ) +
    geom_line(
      data = empirical,
      aes(
        x = b_squared, y = log_value,
        group = interaction(panel, factor(method, levels = c("Rejection", "Independent PM", "CPM"))),
        color = color, linetype = linetype, linewidth = linewidth
      )
    ) +
    geom_point(
      data = empirical,
      aes(x = b_squared, y = log_value, color = color, shape = shape),
      size = 2.8,
      stroke = 0.85
    ) +
    geom_rect(
      data = legend_boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white",
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = custom_legend,
      aes(
        x = 6.2, xend = 10.3, y = y, yend = y,
        color = color, linetype = linetype, linewidth = linewidth
      )
    ) +
    geom_point(
      data = custom_legend[!is.na(custom_legend$shape), ],
      aes(x = 8.25, y = y, color = color, shape = shape),
      size = 2.8,
      stroke = 0.85
    ) +
    geom_text(
      data = custom_legend,
      aes(x = 12, y = y, label = label),
      parse = TRUE,
      hjust = 0,
      family = "sans",
      color = "#222222",
      size = 11.5 / ggplot2::.pt
    ) +
    facet_grid(
      . ~ panel,
      labeller = labeller(panel = as_labeller(panel_labels, label_parsed))
    ) +
    scale_color_identity() +
    scale_shape_identity() +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    scale_x_continuous(
      breaks = (2:8)^2,
      limits = c(4, 64),
      expand = expansion(mult = c(0.025, 0.025))
    ) +
    scale_y_continuous(
      breaks = c(5, 10, 20, 30, 40),
      limits = c(4.5, 40),
      expand = expansion(mult = c(0.012, 0.012))
    ) +
    labs(x = expression(b^2), y = expression(log~N(f[b]*";"~0.10))) +
    theme_classic(base_size = 12.5, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.42),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.7, color = "#3F3F3F"),
      axis.ticks = element_line(linewidth = 0.65, color = "#3F3F3F"),
      axis.ticks.length = unit(4.5, "pt"),
      axis.title = element_text(size = 15.5, color = "#222222"),
      axis.text = element_text(size = 12.5, color = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(
        size = 15.5,
        color = "#222222",
        face = "plain",
        margin = margin(b = 5)
      ),
      panel.spacing.x = unit(12, "pt"),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8, unit = "pt")
    )
}
