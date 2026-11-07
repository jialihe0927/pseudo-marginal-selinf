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

read_result <- function(path) {
  signature <- charToRaw("version https://git-lfs.github.com/spec/v1")
  if (identical(readBin(path, "raw", n = length(signature)), signature)) {
    stop("Git LFS pointer instead of results: ", path, ". Run git lfs pull.")
  }
  readRDS(path)
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  output <- if (length(arguments)) arguments[[1L]] else "../paper-output"
  repository <- script_directory

  # File connections need a UTF-8 locale for non-ASCII Windows paths.
  if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C" &&
      isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
    previous_locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", ".UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
  }
  source(file.path(repository, "diagnostics.R"), local = TRUE)
  for (package in c("posterior", "ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Package '", package, "' is required.")
  }
  if (!capabilities("cairo")) stop("Cairo PDF support is required.")

  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  repository <- normalizePath(repository, winslash = "/", mustWork = TRUE)
  if (tolower(output) == tolower(repository) ||
      startsWith(tolower(output), paste0(tolower(repository), "/"))) {
    stop("Choose an output directory outside the repository to protect the supplied results.")
  }

  # Read the supplied CV-Lasso and Changepoint results and construct their panels.
  applications <- c(cv_lasso = "lasso", changepoint = "cp")
  validation_inputs <- list()
  efficiency_inputs <- list()
  inference_results <- list()
  panels <- list()
  validation_figures <- list()
  experiments <- c(null = "global_null_calibration", nonnull = "power_comparison")
  stems <- c(null = "calibration", nonnull = "power")
  for (application in names(applications)) {
    suffix <- applications[[application]]
    base <- file.path(repository, "experiments", application)
    validation_inputs[[application]] <- read_result(file.path(
      base, "sampling_validation", paste0("sampling_validation_", suffix, ".rds")
    ))
    efficiency_inputs[[application]] <- read_result(file.path(
      base, "efficiency_comparison", paste0("efficiency_", suffix, ".rds")
    ))
    helpers <- new.env(parent = globalenv())
    source(file.path(base, "sampling_validation",
                     paste0("plot_sampling_validation_", suffix, ".R")), local = helpers)
    validation_figures[[application]] <- helpers[[paste0("make_coordinate_figure_", suffix)]](
      validation_inputs[[application]]
    )
    for (experiment in names(experiments)) {
      directory <- file.path(base, experiments[[experiment]])
      stem <- paste0(stems[[experiment]], "_", suffix)
      helpers <- new.env(parent = globalenv())
      if (application == "cv_lasso") {
        source(file.path(repository, "plot_style.R"), local = helpers)
      }
      source(file.path(directory, paste0("plot_", stem, ".R")), local = helpers)
      result <- read_result(file.path(directory, paste0(stem, ".rds")))
      make_panel <- if (experiment == "null") helpers$make_inference_figure else helpers$make_directional_figure
      panels[[application]][[experiment]] <- make_panel(result)
      inference_results[[application]][[experiment]] <- result
    }
  }
  # File Drawer is computed here; it has no saved result input.
  helpers <- new.env(parent = globalenv())
  for (file in c("samplers_fd.R", "analysis_fd.R", "run_fd.R", "plot_fd.R")) {
    source(file.path(repository, "experiments/file_drawer", file), local = helpers)
  }
  file_drawer_result <- helpers$run_experiment_fd()
  file_drawer <- helpers$make_figure(file_drawer_result)
  bank <- read_result(file.path(repository, "inference_trajectories.rds"))

  # Coordinate diagnostics serve both efficiency and sampling validation.
  efficiency_tables <- list()
  coordinates <- list()
  validation_rows <- list()
  power_rows <- list()
  null_rows <- list()
  method_labels <- c("Crude PMMH" = "PM", "Correlated PMMH" = "CPM")
  for (application in names(applications)) {
    suffix <- applications[[application]]
    efficiency <- recompute_efficiency(
      efficiency_inputs[[application]], validation_inputs[[application]]
    )
    table <- efficiency_table(efficiency)
    efficiency_tables[[paste0("efficiency_", suffix)]] <- table
    coordinates[[application]] <- efficiency$coordinate_ess
    for (method in names(method_labels)) {
      rows <- coordinates[[application]][coordinates[[application]]$method == method, ]
      validation_rows[[length(validation_rows) + 1L]] <- data.frame(
        Application = application,
        Method = unname(method_labels[method]),
        minimum_mean_ESS = min(rows$ess_mean),
        maximum_mean_MCSE = max(rows$mcse_mean)
      )
    }

    null <- inference_results[[application]]$null$summary
    names(null)[names(null) == "type_i_error"] <- "rejection_rate"
    null_rows[[application]] <- data.frame(application = application, null)
    power <- inference_results[[application]]$nonnull$summary$power
    power$application <- application
    power_rows[[application]] <- power
  }

  # Diagnose the actual upper-tail indicators from the saved inference trajectories.
  diagnostics <- recompute_inference_diagnostics(bank)
  pooled <- summarize_inference_diagnostics(diagnostics)
  row_order <- order(
    match(pooled$application, c("CV-Lasso", "Changepoint")),
    match(pooled$experiment, c("quadratic_null", "directional")),
    match(pooled$method_id, c(
      "model_only", "fixed_lambda", "full_selector", "fixed_k", "full_conditioning"
    ))
  )
  pooled <- pooled[row_order, ]
  paper_table <- pooled[c(
    "application", "experiment", "method_id", "test_records",
    "median_indicator_ess", "median_p_value_mcse"
  )]
  labels <- c(
    model_only = "Proposed", full_selector = "Proposed", fixed_lambda = "Fixed-lambda",
    fixed_k = "Fixed-K", full_conditioning = "Extra conditioning"
  )
  paper_table$method_id <- unname(labels[paper_table$method_id])
  paper_table$experiment <- ifelse(
    paper_table$experiment == "quadratic_null", "Quadratic null", "Directional tests"
  )
  summaries <- c(efficiency_tables, list(
    sampling_validation = do.call(rbind, validation_rows),
    directional_rejection = do.call(rbind, power_rows),
    quadratic_null_rejection = do.call(rbind, null_rows),
    inference_diagnostics = paper_table
  ))

  # Tables stay on the console; only paper PDFs are written to the output directory.
  for (name in names(summaries)) {
    cat("\n", name, "\n", sep = "")
    print(summaries[[name]], row.names = FALSE, digits = 6)
  }
  ggplot2::ggsave(file.path(output, "sampling_complexity_fd.pdf"), file_drawer,
                  width = 12, height = 5, device = grDevices::cairo_pdf)
  for (application in names(applications)) {
    suffix <- applications[[application]]
    left <- panels[[application]]$null
    right <- panels[[application]]$nonnull
    ggplot2::ggsave(file.path(output, paste0("calibration_", suffix, ".pdf")),
                    left, width = 3.78, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
    ggplot2::ggsave(file.path(output, paste0("power_", suffix, ".pdf")),
                    right, width = 4.62, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
    # Reuse the same calibration and power panels in the combined figure.
    combined <- patchwork::wrap_plots(left, right, ncol = 2, widths = c(4.5, 5.5))
    label <- if (application == "cv_lasso") "CV" else "CP"
    ggplot2::ggsave(file.path(output, paste0("Figure_", label, "_Inference_combined.pdf")),
                    combined, width = 8.4, height = 3.6, device = grDevices::cairo_pdf, bg = "white")
    ggplot2::ggsave(file.path(output, paste0("sampling_validation_", suffix, ".pdf")),
                    validation_figures[[application]], width = 12,
                    height = if (application == "cv_lasso") 8 else 4,
                    device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  }
  invisible(list(
    summaries = summaries, coordinate_ess = coordinates,
    inference_diagnostics = diagnostics, file_drawer = file_drawer_result
  ))
}

if (sys.nframe() == 0L) {
  main()
}
