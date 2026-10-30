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

write_table_csv <- function(data, stem) {
  csv <- data
  for (name in names(csv)) {
    if (is.double(csv[[name]])) {
      csv[[name]] <- ifelse(
        is.finite(csv[[name]]),
        sprintf("%.17g", csv[[name]]),
        NA_character_
      )
    }
  }
  utils::write.csv(csv, paste0(stem, ".csv"), row.names = FALSE, na = "NA")
}

main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if (length(arguments) > 1L) stop("Usage: Rscript reproduce.R [output-directory]")
  output <- if (length(arguments)) arguments[[1L]] else NULL
  repository <- script_directory

  # File connections need a UTF-8 locale for non-ASCII Windows paths.
  if (.Platform$OS.type == "windows" && Sys.getlocale("LC_CTYPE") == "C" &&
      isTRUE(l10n_info()[["system.codepage"]] == 65001L)) {
    previous_locale <- Sys.getlocale("LC_CTYPE")
    Sys.setlocale("LC_CTYPE", ".UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale), add = TRUE)
  }
  source(file.path(repository, "diagnostics.R"), local = TRUE)
  require_diagnostic_version()

  # Read each frozen object once. The helper files only define functions.
  applications <- c(cv_lasso = "lasso", changepoint = "cp")
  validation_inputs <- list()
  efficiency_inputs <- list()
  inference_results <- list()
  inference_helpers <- list()
  for (application in names(applications)) {
    suffix <- applications[[application]]
    base <- file.path(repository, "experiments", application)
    validation_inputs[[application]] <- readRDS(file.path(
      base, "sampling_validation", paste0("sampling_validation_", suffix, ".rds")
    ))
    efficiency_inputs[[application]] <- readRDS(file.path(
      base, "efficiency_comparison", paste0("efficiency_", suffix, ".rds")
    ))
    helpers <- new.env(parent = globalenv())
    for (experiment in c("null", "nonnull")) {
      directory <- file.path(base, paste0("inference_", experiment))
      source(file.path(directory, paste0("functions_", experiment, "_", suffix, ".R")),
             local = helpers)
      result <- readRDS(file.path(
        directory, paste0("inference_", experiment, "_", suffix, ".rds")
      ))
      helpers[[paste0("validate_result_", experiment, "_", suffix)]](result)
      inference_results[[application]][[experiment]] <- result
    }
    inference_helpers[[application]] <- helpers
  }
  bank <- readRDS(file.path(repository, "inference_trajectories.rds"))

  # Match sampled tests, not physical trajectories: a trajectory may support
  # both a quadratic and a directional test. Analytic baselines are excluded.
  validate_inference_results(bank, inference_results)

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
    if (!isTRUE(all.equal(
      table, efficiency_table(efficiency_inputs[[application]]), tolerance = 1e-12
    ))) {
      stop("Recomputed efficiency differs from the frozen result: ", application)
    }
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

    null <- inference_results[[application]]$null
    for (method in null$config$methods$method_id) {
      rows <- null$dataset_results[null$dataset_results$method_id == method, ]
      null_rows[[length(null_rows) + 1L]] <- data.frame(
        application = application, method_id = method, datasets = nrow(rows),
        rejections = sum(rows$p_value <= null$config$test_level),
        rejection_rate = mean(rows$p_value <= null$config$test_level)
      )
    }
    nonnull <- inference_results[[application]]$nonnull
    summarize <- inference_helpers[[application]][[paste0("summarize_nonnull_", suffix)]]
    power <- summarize(nonnull$dataset_results, nonnull$config)$power
    power$application <- application
    power_rows[[application]] <- power
  }

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

  # Default execution prints paper summaries only; exports require a directory.
  for (name in names(summaries)) {
    cat("\n", name, "\n", sep = "")
    print(summaries[[name]], row.names = FALSE, digits = 6)
  }
  if (!is.null(output)) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    for (name in names(summaries)) {
      write_table_csv(summaries[[name]], file.path(output, name))
    }
    for (application in names(applications)) {
      write_table_csv(
        coordinates[[application]],
        file.path(output, paste0("validation_coordinates_", applications[[application]]))
      )
    }
    write_table_csv(diagnostics, file.path(output, "inference_diagnostics_dataset"))
    write_table_csv(pooled, file.path(output, "inference_diagnostics_pooled"))
    write_table_csv(
      summarize_inference_diagnostics(diagnostics, by_signal = TRUE),
      file.path(output, "inference_diagnostics_by_signal")
    )
    source(file.path(repository, "plot_inference.R"), local = TRUE)
    save_inference_figures(inference_panels(repository, inference_results), output)
  }
  cat(
    "Offline reproduction: ", nrow(diagnostics), " test diagnostics; ",
    length(bank$trajectories), " unique production trajectories. New selector evaluations: 0.\n",
    sep = ""
  )
  invisible(list(
    summaries = summaries, coordinate_ess = coordinates,
    inference_diagnostics = diagnostics
  ))
}

if (sys.nframe() == 0L) {
  main()
}
