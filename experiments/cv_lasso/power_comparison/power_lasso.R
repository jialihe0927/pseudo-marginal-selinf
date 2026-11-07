main <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("../global_null_calibration/functions_calibration_lasso.R", local = TRUE)
  source("../global_null_calibration/samplers_calibration_lasso.R", local = TRUE)
  source("../../../diagnostics.R", local = TRUE)
  for (package in c("glmnet", "posterior")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Package '", package, "' is required.")
  }

  source("functions_power_lasso.R", local = TRUE)
  source("samplers_power_lasso.R", local = TRUE)

  if (length(args) > 2L) stop("Usage: Rscript power_lasso.R [output.rds] [bank.rds]")
  output_path <- if (length(args)) args[[1L]] else "rerun_power_lasso.rds"
  if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

  bank_path <- if (length(args) >= 2L) args[[2L]] else "power_lasso.rds"
  input <- readRDS(bank_path)
  bank <- input[c("config", "datasets", "rng_assignments")]
  if (!is.null(input$raw_attempts)) bank$raw_attempts <- input$raw_attempts
  config <- bank$config
  datasets <- bank$datasets
  assignments <- bank$rng_assignments
  dataset_keys <- vapply(datasets, function(dataset) {
    paste(dataset$delta, dataset$dataset_id, dataset$raw_id, sep = ":")
  }, character(1L))
  dataset_indices <- match(
    paste(assignments$delta, assignments$dataset_id, assignments$raw_id, sep = ":"),
    dataset_keys
  )
  if (anyDuplicated(dataset_keys) || anyNA(dataset_indices)) {
    stop("Sampling tasks do not uniquely match the bank datasets.")
  }

  # Match tasks by signal and dataset identity, not by row position.
  rows <- vector("list", nrow(assignments) + length(datasets))
  for (job in seq_len(nrow(assignments))) {
    assignment <- assignments[job, , drop = FALSE]
    rows[[job]] <- run_pmmh_nonnull_lasso(
      datasets[[dataset_indices[job]]], assignment$method_id, assignment$rng_state[[1L]], config
    )
  }
  # Data splitting uses the same selected datasets and an analytic reference.
  for (index in seq_along(datasets)) {
    dataset <- datasets[[index]]
    score <- fission_statistic_nonnull_lasso(dataset, config)
    p_value <- fission_pvalue_nonnull_lasso(dataset, config, score = score)
    rows[[nrow(assignments) + index]] <- result_row_nonnull_lasso(
      dataset, "data_fission", p_value, config, statistic_observed = score
    )
  }
  offline_output <- collect_inference_output(rows)
  data <- do.call(rbind, rows)
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
      summary = summarize_nonnull_lasso(data, config)
    )
  )
  result <- c(result, offline_output)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output_path, version = 3)
  message("Saved results to ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
