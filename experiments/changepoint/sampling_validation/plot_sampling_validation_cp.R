suppressPackageStartupMessages(library(ggplot2))

# Descriptive coordinate ECDFs compare marginals, without asserting joint agreement or convergence.
make_coordinate_figure_cp <- function(object) {
  methods <- c("Rejection", "Crude PMMH", "Correlated PMMH")
  coordinates <- object$config$coordinates
  colors <- setNames(c("#4d4d4d", "#bc6c65", "#4878a8"), methods)
  shapes <- setNames(c(16L, 17L, 15L), methods)
  lines <- list()
  markers <- list()
  bounds <- list()
  for (coordinate in coordinates) {
    pooled <- unlist(lapply(object$samples, function(x) x[, coordinate]), use.names = FALSE)
    limits <- range(pooled)
    padding <- diff(limits) * 0.035
    lower <- limits[1L] - padding
    upper <- limits[2L] + padding
    bounds[[coordinate]] <- c(lower, upper)
    for (method in methods) {
      x <- object$samples[[method]][, coordinate]
      distribution <- stats::ecdf(x)
      grid <- c(lower, sort(unique(x)), upper)
      lines[[length(lines) + 1L]] <- data.frame(
        variable = coordinate, method = method, x = grid, ecdf = distribution(grid)
      )
      # Markers use quantiles; every state, including repeats, enters the ECDF.
      positions <- unname(stats::quantile(x, seq(0.1, 0.9, length.out = 7L), type = 1))
      markers[[length(markers) + 1L]] <- data.frame(
        variable = coordinate, method = method, x = positions, ecdf = distribution(positions)
      )
    }
  }
  lines <- do.call(rbind, lines)
  markers <- do.call(rbind, markers)
  limits <- unname(bounds[["T1"]])
  span <- diff(limits)
  legend <- data.frame(
    variable = "T1", method = methods,
    x = limits[1L] + 0.055 * span,
    xend = limits[1L] + 0.155 * span,
    point_x = limits[1L] + 0.105 * span,
    text_x = limits[1L] + 0.185 * span,
    y = c(0.952, 0.887, 0.822),
    label = unname(c(Rejection = "Rejection", "Crude PMMH" = "PM", "Correlated PMMH" = "CPM")[methods])
  )
  box <- data.frame(
    variable = "T1", xmin = limits[1L] + 0.025 * span,
    xmax = limits[1L] + 0.43 * span, ymin = 0.788, ymax = 0.986
  )
  lines$variable <- factor(lines$variable, levels = coordinates)
  markers$variable <- factor(markers$variable, levels = coordinates)
  legend$variable <- factor(legend$variable, levels = coordinates)
  box$variable <- factor(box$variable, levels = coordinates)
  lines$method <- factor(lines$method, levels = methods)
  markers$method <- factor(markers$method, levels = methods)
  labels <- setNames(paste0("'(", letters[seq_along(coordinates)], ")'~t[",
                            seq_along(coordinates), "]"), coordinates)
  ggplot() +
    geom_step(data = lines, aes(x = x, y = ecdf, color = method, group = method),
              direction = "hv", linewidth = 1.3, lineend = "butt") +
    geom_point(data = markers, aes(x = x, y = ecdf, color = method, shape = method),
               size = 2.8, stroke = 0.85) +
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", color = "black", linewidth = 0.45) +
    geom_segment(data = legend, aes(x = x, xend = xend, y = y, yend = y, color = method),
                 linewidth = 1.3) +
    geom_point(data = legend, aes(x = point_x, y = y, color = method, shape = method),
               size = 2.8, stroke = 0.85) +
    geom_text(data = legend, aes(x = text_x, y = y, label = label),
              family = "sans", hjust = 0, size = 11.5 / ggplot2::.pt, color = "#222222") +
    facet_grid(. ~ variable, scales = "free_x",
               labeller = labeller(variable = as_labeller(labels, label_parsed))) +
    scale_color_manual(values = colors, breaks = methods) +
    scale_shape_manual(values = shapes, breaks = methods) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 4L),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.012, 0.012))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Empirical CDF") +
    theme_classic(base_size = 12.5, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.42),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.7, color = "#3F3F3F"),
      axis.ticks = element_line(linewidth = 0.65, color = "#3F3F3F"),
      axis.ticks.length = grid::unit(4.5, "pt"),
      axis.title = element_text(size = 15.5, color = "#222222"),
      axis.text = element_text(size = 12.5, color = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(size = 15.5, color = "#222222",
                                face = "plain", margin = margin(b = 5)),
      panel.spacing.x = grid::unit(12, "pt"),
      legend.position = "none", plot.margin = margin(8, 8, 8, 8, unit = "pt")
    )
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (!capabilities("cairo")) stop("Cairo PDF support is required for the paper figure.")
  if (length(arguments) > 1L) stop("Usage: Rscript plot_sampling_validation_cp.R [input.rds]")
  input_path <- if (length(arguments)) arguments[[1L]] else "sampling_validation_cp.rds"
  object <- readRDS(input_path)
  figure <- make_coordinate_figure_cp(object)
  ggsave(
    paste0(tools::file_path_sans_ext(input_path), ".pdf"), figure,
    width = 12, height = 4, units = "in",
    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
}

if (sys.nframe() == 0L) main()
