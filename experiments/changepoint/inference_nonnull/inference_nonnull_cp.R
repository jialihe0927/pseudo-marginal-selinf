script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))

main <- function() {
  options(stringsAsFactors = FALSE)

  source(file.path(script_directory, "../inference_null/functions_null_cp.R"), local = TRUE)
  source(file.path(script_directory, "../inference_null/samplers_null_cp.R"), local = TRUE)
  source(file.path(script_directory, "../../../diagnostics.R"), local = TRUE)
  require_diagnostic_version()

  source(file.path(script_directory, "functions_nonnull_cp.R"), local = TRUE)
  source(file.path(script_directory, "samplers_nonnull_cp.R"), local = TRUE)

  # Select an output and a frozen or user-supplied bank.
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: Rscript inference_nonnull_cp.R [output.rds] [bank.rds]")
  output_path <- if (length(args)) args[[1L]] else "rerun_inference_nonnull_cp.rds"
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

  # All model settings and task streams come from this bank.
  bank_path <- if (length(args) >= 2L) args[[2L]] else file.path(script_directory, "inference_nonnull_cp.rds")
  if (!file.exists(bank_path)) stop("Data bank not found: ", bank_path)
  input <- readRDS(bank_path)
  bank <- input[c("config", "datasets", "rng_assignments")]
  if (!is.null(input$raw_attempts)) bank$raw_attempts <- input$raw_attempts
  validate_bank_nonnull_cp(bank)
  config <- bank$config
  datasets <- bank$datasets
  assignments <- bank$rng_assignments
  dataset_keys <- vapply(datasets, function(dataset) {
    paste(dataset$delta, dataset$dataset_id, sep = ":")
  }, character(1L))

  # Match tasks by signal and dataset identity, not by row position.
  rows <- vector("list", nrow(assignments) + length(datasets))
  for (job in seq_len(nrow(assignments))) {
    assignment <- assignments[job, , drop = FALSE]
    index <- match(
      paste(assignment$delta, assignment$dataset_id, sep = ":"),
      dataset_keys
    )
    message(
      "delta=", assignment$delta, ", dataset=", assignment$dataset_id,
      "/", config$datasets_per_delta, ", method=", assignment$method_id
    )
    rows[[job]] <- run_pmmh_nonnull_cp(
      datasets[[index]], assignment$method_id, assignment$rng_state[[1L]], config
    )
  }
  # Data splitting uses the same selected datasets and an analytic reference.
  for (index in seq_along(datasets)) {
    dataset <- datasets[[index]]
    score <- fission_statistic_nonnull_cp(dataset, config)
    p_value <- fission_pvalue_nonnull_cp(dataset, config, score)
    rows[[nrow(assignments) + index]] <- result_row_nonnull_cp(
      dataset, "data_fission", p_value, config,
      result_type = if (dataset$direction$structural_zero) "structural_zero" else "analytic",
      statistic_observed = score
    )
  }
  offline_output <- collect_inference_output(rows)
  data <- do.call(rbind, rows)
  # Restore paired signal/dataset/method order, then summarize.
  data <- data[
    order(
      match(data$delta, config$delta_grid),
      data$dataset_id,
      match(data$method_id, config$methods)
    ),
    ,
    drop = FALSE
  ]
  row.names(data) <- NULL
  result <- c(
    bank,
    list(
      dataset_results = data,
      summary = summarize_nonnull_cp(data, config)
    )
  )
  validate_result_nonnull_cp(result)
  result <- c(result, offline_output)
  saveRDS(result, output_path, version = 3)
  message("Saved results to ", normalizePath(output_path, winslash = "/"))
}

if (sys.nframe() == 0L) main()
