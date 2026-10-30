require_diagnostic_version <- function(strict_r = FALSE) {
  if (as.character(getRversion()) != "4.4.3") {
    description <- paste0(
      "Paper diagnostics were verified with R 4.4.3; found ", getRversion(), "."
    )
    if (strict_r) stop(description)
    warning(description, call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE) ||
      as.character(utils::packageVersion("posterior")) != "1.7.0") {
    stop("These diagnostics require posterior 1.7.0; no package is installed automatically.")
  }
  invisible(TRUE)
}

mean_diagnostics <- function(x, label = "statistic") {
  require_diagnostic_version()
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("Nonfinite production input: ", label)
  reason <- if (length(x) < 4L) {
    "insufficient observations"
  } else if (length(unique(x)) < 2L) {
    "constant sequence"
  } else ""
  if (nzchar(reason)) {
    return(list(ess_mean = NA_real_, mcse_mean = NA_real_, tau = NA_real_,
                undefined_reason = reason, package_warnings = ""))
  }

  # A single statistic from one ordered trajectory: iterations by one chain.
  draws <- matrix(as.numeric(x), ncol = 1L)
  messages <- character()
  estimates <- withCallingHandlers({
    list(
      ess_mean = as.numeric(posterior::ess_basic(draws, split = TRUE)),
      mcse_mean = as.numeric(posterior::mcse_mean(draws))
    )
  }, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    message(label, ": ", conditionMessage(w))
  })
  if (length(estimates$ess_mean) != 1L || length(estimates$mcse_mean) != 1L) {
    stop("Unexpected posterior result shape: ", label)
  }
  reasons <- character()
  if (!is.finite(estimates$ess_mean) || estimates$ess_mean <= 0) {
    reasons <- c(reasons, paste("undefined mean ESS:", estimates$ess_mean))
    estimates$ess_mean <- NA_real_
  }
  if (!is.finite(estimates$mcse_mean) || estimates$mcse_mean <= 0) {
    reasons <- c(reasons, paste("undefined mean MCSE:", estimates$mcse_mean))
    estimates$mcse_mean <- NA_real_
  }
  if (is.finite(estimates$ess_mean) && is.finite(estimates$mcse_mean) &&
      !isTRUE(all.equal(estimates$mcse_mean, stats::sd(x) / sqrt(estimates$ess_mean),
                       tolerance = 1e-10))) {
    stop("Posterior mean ESS/MCSE identity differs: ", label)
  }
  c(estimates, list(
    tau = length(x) / estimates$ess_mean,
    undefined_reason = paste(reasons, collapse = "; "),
    package_warnings = paste(unique(messages), collapse = " | ")
  ))
}

finite_median <- function(x) {
  if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
}

collect_inference_output <- function(rows) {
  sampled <- Filter(function(row) !is.null(attr(row, "production_statistics")), rows)
  records <- lapply(sampled, function(row) {
    identity <- row[intersect(c("delta", "dataset_id", "raw_id", "method_id"), names(row))]
    cbind(identity, as.data.frame(attr(row, "indicator_diagnostics")))
  })
  list(
    indicator_diagnostics = if (length(records)) do.call(rbind, records) else data.frame(),
    production_statistics = lapply(sampled, function(row) {
      list(identity = row[intersect(c("delta", "dataset_id", "raw_id", "method_id"), names(row))],
           observed = row$statistic_observed, values = attr(row, "production_statistics"))
    })
  )
}

coordinate_mean_diagnostics <- function(draws, method) {
  if (!is.matrix(draws) || any(!is.finite(draws))) stop("Invalid draws for ", method)
  coordinates <- colnames(draws)
  if (is.null(coordinates)) coordinates <- paste0("T", seq_len(ncol(draws)))
  do.call(rbind, lapply(seq_len(ncol(draws)), function(j) {
    estimate <- mean_diagnostics(draws[, j], paste(method, coordinates[j]))
    data.frame(method = method, coordinate = coordinates[j], n = nrow(draws),
               ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean,
               tau = estimate$tau, undefined_reason = estimate$undefined_reason,
               package_warnings = estimate$package_warnings)
  }))
}

recompute_efficiency <- function(result, validation) {
  methods <- c("Crude PMMH", "Correlated PMMH")
  required <- c("Rejection", methods)
  sample_names <- names(validation$samples)
  if (anyDuplicated(sample_names) || !setequal(sample_names, required) ||
      !all(required %in% names(validation$sample_info))) {
    stop("Efficiency requires named Rejection, Crude PMMH and Correlated PMMH samples and metadata.")
  }
  if (!isTRUE(all.equal(result$benchmark, validation$config$benchmark, tolerance = 1e-12)) ||
      !isTRUE(all.equal(result$benchmark$Sigma_T,
                        validation$config$target$covariance, tolerance = 1e-12)) ||
      !isTRUE(all.equal(result$benchmark$gaussian_center,
                        validation$config$target$gaussian_center, tolerance = 1e-12))) {
    stop("Efficiency and validation must use the same benchmark and conditional target.")
  }
  coordinate <- do.call(rbind, lapply(methods, function(method) {
    draws <- validation$samples[[method]]
    info <- validation$sample_info[[method]]
    if (!identical(colnames(draws), validation$config$coordinates) ||
        nrow(draws) != info$sample_count ||
        !identical(info$target_id, validation$config$target_id)) {
      stop("Production coordinates, sample count or target identity differ: ", method)
    }
    coordinate_mean_diagnostics(draws, method)
  }))
  for (method in methods) {
    selected <- result$estimates$method == method
    if (sum(selected) != 1L) stop("Expected one efficiency estimate for ", method)
    minimum <- min(coordinate$ess_mean[coordinate$method == method])
    info <- validation$sample_info[[method]]
    if (info$selector_evaluations != result$estimates$selector_evaluations[selected] ||
        abs(mean(info$accepted) - result$estimates$production_acceptance[selected]) > 1e-12 ||
        abs(mean(info$proposed_zero) - result$estimates$zero_estimate_rate[selected]) > 1e-12 ||
        abs(info$final_scale - result$estimates$final_proposal_scale[selected]) > 1e-12) {
      stop("Frozen sampling cost, acceptance, mismatch or proposal scale differ: ", method)
    }
    result$estimates$effective_samples[selected] <- minimum
    result$estimates$efficiency_per_1000[selected] <-
      1000 * minimum / info$selector_evaluations
  }
  result$coordinate_ess <- coordinate
  result$config$ess_definition <- "minimum coordinate posterior::ess_basic(split = TRUE)"
  result$config$ess <- result$config$ess_definition
  result$config$diagnostic_software <- c(
    R = as.character(getRversion()), posterior = as.character(utils::packageVersion("posterior"))
  )
  result
}

efficiency_table <- function(result) {
  ids <- c("Rejection", "Crude PMMH", "Correlated PMMH")
  rows <- result$estimates[match(ids, result$estimates$method), ]
  if (anyNA(rows$method)) stop("The efficiency table requires all three methods.")
  data.frame(
    Method = c("Rejection", "PM", "CPM"), Efficiency = rows$efficiency_per_1000,
    Gain = rows$efficiency_per_1000 / rows$efficiency_per_1000[1L],
    Acceptance = rows$production_acceptance, Mismatch = rows$zero_estimate_rate,
    h_final = rows$final_proposal_scale
  )
}

validate_inference_trajectories <- function(bank) {
  tests <- bank$tests
  required <- c("application", "experiment", "delta", "dataset_id", "raw_id",
                "method_id", "trajectory_key", "statistic_id", "statistic_observed",
                "official_p_value")
  if (!is.data.frame(tests) || !all(required %in% names(tests)) ||
      anyNA(tests[required]) || any(!is.finite(tests$statistic_observed))) {
    stop("Incomplete inference test identities or observed statistics.")
  }
  keys <- paste(tests$application, tests$experiment, tests$delta,
                tests$method_id, tests$dataset_id, sep = "|")
  if (anyDuplicated(keys) || anyDuplicated(names(bank$trajectories)) ||
      !setequal(tests$trajectory_key, names(bank$trajectories))) {
    stop("Incomplete or duplicate test/trajectory identities.")
  }
  if (any(!is.finite(tests$official_p_value)) ||
      any(tests$official_p_value < 0 | tests$official_p_value > 1)) {
    stop("Invalid frozen p-values.")
  }
  invisible(bank)
}

validate_inference_results <- function(bank, results, tolerance = 1e-12) {
  validate_inference_trajectories(bank)
  expected <- list()
  applications <- c(cv_lasso = "CV-Lasso", changepoint = "Changepoint")
  experiments <- c(null = "quadratic_null", nonnull = "directional")
  for (application in names(applications)) {
    for (experiment in names(experiments)) {
      result <- results[[application]][[experiment]]
      data <- result$dataset_results
      if (is.null(data) || !all(c("result_type", "method_id", "dataset_id", "raw_id",
                                  "statistic_observed", "p_value") %in% names(data))) {
        stop("Missing inference results: ", application, "/", experiment)
      }
      # Analytic baselines have no production trajectory; each sampled test does.
      sampled <- data[data$result_type == "simulation", , drop = FALSE]
      expected[[length(expected) + 1L]] <- data.frame(
        application = unname(applications[application]),
        experiment = unname(experiments[experiment]),
        delta = if (experiment == "null") 0 else sampled$delta,
        dataset_id = sampled$dataset_id,
        raw_id = sampled$raw_id,
        method_id = sampled$method_id,
        statistic_observed = sampled$statistic_observed,
        p_value = sampled$p_value
      )
      if (result$config$pmmh$production != bank$config$production ||
          result$config$test_level != bank$config$test_level) {
        stop("Inference archive and result configurations differ: ", application, "/", experiment)
      }
    }
  }
  expected <- do.call(rbind, expected)
  tests <- bank$tests
  test_key <- function(data) {
    paste(data$application, data$experiment, data$delta,
          data$dataset_id, data$method_id, sep = "|")
  }
  keys <- test_key(tests)
  expected_keys <- test_key(expected)
  if (anyDuplicated(expected_keys)) {
    stop("Duplicate sampled result identity: ", expected_keys[duplicated(expected_keys)][1L])
  }
  missing <- setdiff(expected_keys, keys)
  extra <- setdiff(keys, expected_keys)
  if (length(missing)) stop("Inference trajectory record is missing: ", missing[1L])
  if (length(extra)) stop("Inference trajectory record has no sampled result: ", extra[1L])
  matched <- expected[match(keys, expected_keys), , drop = FALSE]
  if (anyNA(matched) || any(!is.finite(matched$p_value)) ||
      any(!is.finite(matched$statistic_observed))) {
    stop("Nonfinite or missing sampled inference results.")
  }
  mismatched_raw <- which(tests$raw_id != matched$raw_id)
  if (length(mismatched_raw)) {
    stop("Inference raw_id differs for ", keys[mismatched_raw[1L]])
  }
  # P-values use absolute tolerance; statistics use the same relative precision.
  comparisons <- list(
    p_value = abs(tests$official_p_value - matched$p_value) > tolerance,
    statistic_observed = abs(tests$statistic_observed - matched$statistic_observed) >
      tolerance * pmax(1, abs(tests$statistic_observed), abs(matched$statistic_observed))
  )
  for (name in names(comparisons)) {
    mismatched <- which(comparisons[[name]])
    if (length(mismatched)) {
      stop("Inference ", name, " differs between trajectory archive and result: ",
           keys[mismatched[1L]])
    }
  }
  invisible(TRUE)
}

recompute_inference_diagnostics <- function(bank) {
  validate_inference_trajectories(bank)
  rows <- lapply(seq_len(nrow(bank$tests)), function(i) {
    test <- bank$tests[i, ]
    series <- bank$trajectories[[test$trajectory_key]][, test$statistic_id]
    if (length(series) != bank$config$production || any(!is.finite(series))) {
      stop("Invalid production statistic for ", test$trajectory_key)
    }
    indicator <- as.numeric(series >= test$statistic_observed)
    p_value <- mean(indicator)
    if (abs(p_value - test$official_p_value) > 1e-12) {
      stop("Tail mean differs from the frozen p-value: ", test$trajectory_key)
    }
    estimate <- mean_diagnostics(indicator, paste(test$trajectory_key, test$statistic_id))
    cbind(test, data.frame(
      production_draws = length(indicator), tail_count = sum(indicator),
      p_value = p_value, reject = p_value <= bank$config$test_level,
      ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean,
      tau = estimate$tau, undefined_reason = estimate$undefined_reason,
      package_warnings = estimate$package_warnings
    ))
  })
  do.call(rbind, rows)
}

summarize_inference_diagnostics <- function(data, by_signal = FALSE) {
  groups <- unique(data[c("application", "experiment", "method_id")])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, ]
    selected <- data[data$application == group$application &
                       data$experiment == group$experiment & data$method_id == group$method_id, ]
    signals <- if (by_signal) sort(unique(selected$delta)) else NA_real_
    for (delta in signals) {
      block <- if (by_signal) selected[selected$delta == delta, ] else selected
      rows[[length(rows) + 1L]] <- cbind(group, data.frame(
        delta = delta, test_records = nrow(block),
        median_indicator_ess = finite_median(block$ess_mean),
        median_p_value_mcse = finite_median(block$mcse_mean),
        finite_ess = sum(is.finite(block$ess_mean)),
        finite_mcse = sum(is.finite(block$mcse_mean))
      ))
    }
  }
  do.call(rbind, rows)
}

