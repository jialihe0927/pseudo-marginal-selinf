#!/usr/bin/env Rscript

script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
# Windows can leave UTF-8 command-line paths unmarked under the C locale.
if (.Platform$OS.type == "windows" &&
    isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
  Encoding(script_file) <- "UTF-8"
}
script_directory <- dirname(script_file)

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

  required <- c("regime", "b", "method", "complexity_n", "resolved")
  if (!is.list(result) || !is.data.frame(result$estimates)) {
    stop("The RDS file must contain a list with an estimates data frame.")
  }
  if (!all(required %in% names(result$estimates))) {
    stop("The estimates data frame is missing required columns.")
  }

  estimates <- result$estimates
  expected <- rbind(
    expand.grid(
      regime = "fixed_a",
      b = 2:8,
      method = c("Rejection", "Independent PM", "CPM"),
      stringsAsFactors = FALSE
    ),
    expand.grid(
      regime = "adaptive_a",
      b = 2:8,
      method = c("Rejection", "Independent PM"),
      stringsAsFactors = FALSE
    )
  )
  key <- function(x) paste(x$regime, x$b, x$method, sep = "|")
  estimates <- estimates[key(estimates) %in% key(expected), ]
  if (nrow(estimates) != nrow(expected) || anyDuplicated(key(estimates)) ||
      !setequal(key(estimates), key(expected))) {
    stop("The RDS file does not contain exactly the displayed configurations.")
  }
  if (any(!estimates$resolved) || any(!is.finite(estimates$complexity_n)) ||
      any(estimates$complexity_n <= 0)) {
    stop("Every displayed complexity estimate must be resolved and positive.")
  }

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
  fixed_legend <- make_legend(
    "fixed_a",
    legend_order,
    c("'Rejection'", "'PM'", "'CPM'", "slope~0.5", "slope~0.05", "log(b^2)")
  )
  adaptive_keys <- c(
    "Rejection", "Independent PM", "Rejection reference",
    "Independent PM reference"
  )
  adaptive_legend <- make_legend(
    "adaptive_a",
    adaptive_keys,
    c("'Rejection'", "'PM'", "slope~0.5", "slope~0")
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

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  # R file connections need a UTF-8 locale for non-ASCII Windows paths.
  if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C" &&
      isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
    previous_locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", ".UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
  }

  if (length(args) > 1L) stop("Usage: Rscript plot_fd.R [results.rds]")
  input_path <- if (length(args)) {
    args[[1L]]
  } else {
    file.path(script_directory, "results_fd.rds")
  }
  if (!file.exists(input_path)) stop("Results file not found: ", input_path)
  if (!capabilities("cairo")) stop("PDF plotting requires R with Cairo support.")

  output_stem <- if (!length(args)) {
    file.path(script_directory, "sampling_complexity_fd")
  } else if (identical(input_path, "results_fd.rds")) {
    "sampling_complexity_fd"
  } else {
    tools::file_path_sans_ext(input_path)
  }
  figure <- make_figure(readRDS(input_path))
  ggplot2::ggsave(
    paste0(output_stem, ".pdf"),
    figure,
    width = 12,
    height = 5,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

if (sys.nframe() == 0L) {
  main()
}
