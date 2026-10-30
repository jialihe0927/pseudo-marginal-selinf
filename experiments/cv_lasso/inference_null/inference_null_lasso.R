script_file <- if (sys.nframe() > 0L) {
  files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  tail(files[nzchar(files)], 1L)
} else {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
}
script_directory <- dirname(script_file)

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options(stringsAsFactors = FALSE)

  source(file.path(script_directory, "../../../diagnostics.R"), local = TRUE)
  require_diagnostic_version()

  source(file.path(script_directory, "functions_null_lasso.R"), local = TRUE)
  source(file.path(script_directory, "samplers_null_lasso.R"), local = TRUE)

  # Choose the output and bank without overwriting saved results.
  if (length(args) > 2L) {
    stop("Usage: Rscript inference_null_lasso.R [output.rds] [bank.rds]")
  }
  output_path <- if (length(args)) args[[1L]] else "rerun_inference_null_lasso.rds"
  if (file.exists(output_path) || dir.exists(output_path)) {
    stop("Output already exists: ", output_path)
  }

  # Keep the bank's configuration, dataset order and task streams together.
  bank_path <- if (length(args) >= 2L) args[[2L]] else file.path(script_directory, "inference_null_lasso.rds")
  if (!file.exists(bank_path)) stop("Data bank not found: ", bank_path)
  bank_source <- readRDS(bank_path)
  bank <- bank_source[c("config", "datasets", "rng_assignments")]
  validate_bank_null_lasso(bank)
  config <- bank$config
  datasets <- bank$datasets
  assignments <- bank$rng_assignments
  if (length(args) >= 2L && !is.null(bank_source$dataset_results)) {
    warning("The supplied bank contains results; only config, datasets, and rng_assignments are used.", call. = FALSE)
  }

  # Run each selective test with its assigned stream.
  sampled_results <- vector("list", nrow(assignments))
  for (job in seq_len(nrow(assignments))) {
    assignment <- assignments[job, , drop = FALSE]
    message(
      sprintf(
        "Dataset %d/%d, %s",
        assignment$dataset_id,
        config$dataset_count,
        assignment$method_id
      )
    )
    sampled_results[[job]] <- run_pmmh_null_lasso(
      dataset = datasets[[assignment$dataset_id]],
      method_id = assignment$method_id,
      rng_state = assignment$rng_state[[1L]],
      config = config
    )
  }
  offline_output <- collect_inference_output(sampled_results)
  sampled_results <- do.call(rbind, sampled_results)

  # The paired naive comparison is analytic and consumes no random draws.
  naive_results <- do.call(
    rbind,
    lapply(datasets, function(dataset) {
      p_value <- pchisq(
        dataset$target$Q_obs,
        df = dataset$target_dimension,
        lower.tail = FALSE
      )
      data.frame(
        dataset_id = dataset$dataset_id,
        raw_id = dataset$raw_id,
        target_id = "selected_model",
        method_id = "naive_gaussian",
        target_dimension = dataset$target_dimension,
        statistic_observed = dataset$target$Q_obs,
        tail_count = NA_integer_,
        production_draws = NA_integer_,
        p_value = p_value,
        reject = as.integer(p_value <= config$test_level),
        result_type = "analytic",
        status = "complete",
        warmup_acceptance = NA_real_,
        production_acceptance = NA_real_,
        final_proposal_scale = NA_real_,
        selector_evaluations = 0L,
        stringsAsFactors = FALSE
      )
    })
  )

  # Restore dataset/method order before summarizing and saving.
  dataset_results <- rbind(sampled_results, naive_results)
  dataset_results <- dataset_results[
    order(
      dataset_results$dataset_id,
      match(dataset_results$method_id, config$methods$method_id)
    ),
    ,
    drop = FALSE
  ]
  row.names(dataset_results) <- NULL

  result <- list(
    config = config,
    datasets = datasets,
    rng_assignments = assignments,
    dataset_results = dataset_results,
    summary = summarize_null_results_lasso(
      dataset_results,
      config$test_level,
      config$methods$method_id
    )
  )
  validate_result_null_lasso(result)
  result <- c(result, offline_output)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output_path, version = 3)
  message("Saved results to ", normalizePath(output_path, winslash = "/"))
}

if (sys.nframe() == 0L) {
  main()
}
