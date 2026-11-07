main <- function(args = commandArgs(trailingOnly = TRUE)) {
  source("../../../diagnostics.R", local = TRUE)
  for (package in c("glmnet", "posterior")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Package '", package, "' is required.")
  }

  source("functions_calibration_lasso.R", local = TRUE)
  source("samplers_calibration_lasso.R", local = TRUE)

  if (length(args) > 2L) {
    stop("Usage: Rscript calibration_lasso.R [output.rds] [bank.rds]")
  }
  output_path <- if (length(args)) args[[1L]] else "rerun_calibration_lasso.rds"
  if (file.exists(output_path) || dir.exists(output_path)) {
    stop("Output already exists: ", output_path)
  }

  bank_path <- if (length(args) >= 2L) args[[2L]] else "calibration_lasso.rds"
  bank_source <- readRDS(bank_path)
  bank <- bank_source[c("config", "datasets", "rng_assignments")]
  config <- bank$config
  datasets <- bank$datasets
  assignments <- bank$rng_assignments
  dataset_keys <- vapply(datasets, function(dataset) {
    paste(dataset$dataset_id, dataset$raw_id, sep = ":")
  }, character(1L))
  dataset_indices <- match(
    paste(assignments$dataset_id, assignments$raw_id, sep = ":"), dataset_keys
  )
  if (anyDuplicated(dataset_keys) || anyNA(dataset_indices)) {
    stop("Sampling tasks do not uniquely match the bank datasets.")
  }

  sampled_results <- vector("list", nrow(assignments))
  for (job in seq_len(nrow(assignments))) {
    assignment <- assignments[job, , drop = FALSE]
    sampled_results[[job]] <- run_pmmh_null_lasso(
      dataset = datasets[[dataset_indices[job]]],
      method_id = assignment$method_id,
      rng_state = assignment$rng_state[[1L]],
      config = config
    )
  }
  offline_output <- collect_inference_output(sampled_results)
  sampled_results <- do.call(rbind, sampled_results)

  # Naive uses the unadjusted chi-square upper tail on the same selected datasets.
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
  result <- c(result, offline_output)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output_path, version = 3)
  message("Saved results to ", output_path)
}

if (sys.nframe() == 0L) {
  main()
}
