script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
if (.Platform$OS.type == "windows" &&
    isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
  Encoding(script_file) <- "UTF-8"
}
script_directory <- dirname(script_file)

# Compose the same independently runnable panel objects used by the experiments.
# The offline summary entry can supply its already-loaded inference results.
inference_panels <- function(repository = ".", results = NULL) {
  panels <- list()
  for (application in c("cv_lasso", "changepoint")) {
    suffix <- if (application == "cv_lasso") "lasso" else "cp"
    base <- file.path(repository, "experiments", application)
    for (experiment in c("null", "nonnull")) {
      directory <- file.path(base, paste0("inference_", experiment))
      environment <- new.env(parent = globalenv())
      source(file.path(directory, paste0("plot_inference_", experiment, "_", suffix, ".R")),
             local = environment)
      result <- if (is.null(results)) {
        readRDS(file.path(directory, paste0("inference_", experiment, "_", suffix, ".rds")))
      } else {
        results[[application]][[experiment]]
      }
      make_panel <- if (experiment == "null") environment$make_inference_figure else environment$make_directional_figure
      panels[[paste(application, experiment, sep = ":")]] <- make_panel(result)
    }
  }
  panels
}

save_inference_figures <- function(panels, output_directory) {
  if (!capabilities("cairo")) stop("Cairo PDF support is required.")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  for (application in c("cv_lasso", "changepoint")) {
    suffix <- if (application == "cv_lasso") "lasso" else "cp"
    left <- panels[[paste0(application, ":null")]]
    right <- panels[[paste0(application, ":nonnull")]]
    ggplot2::ggsave(file.path(output_directory, paste0("inference_null_", suffix, ".pdf")),
                    left, width = 3.78, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
    ggplot2::ggsave(file.path(output_directory, paste0("inference_nonnull_", suffix, ".pdf")),
                    right, width = 4.62, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
    combined <- patchwork::wrap_plots(left, right, ncol = 2, widths = c(4.5, 5.5))
    label <- if (application == "cv_lasso") "CV" else "CP"
    ggplot2::ggsave(file.path(output_directory, paste0("Figure_", label, "_Inference_combined.pdf")),
                    combined, width = 8.4, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
  }
  invisible(panels)
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (length(arguments) != 1L) stop("Usage: Rscript plot_inference.R output-directory")
  if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C" &&
      isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
    previous_locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", ".UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
  }
  output <- arguments[[1L]]
  save_inference_figures(inference_panels(script_directory), output)
}

if (sys.nframe() == 0L) {
  main()
}
