# Mean ESS/MCSE for an untransformed statistic or the actual 0/1 tail indicator.
# Short or constant series retain NA estimates with an explicit reason.
mean_diagnostics <- function(x, label = "statistic") {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("Nonfinite input: ", label)
  reason <- if (length(x) < 4L) {
    "insufficient observations"
  } else if (length(unique(x)) < 2L) {
    "constant sequence"
  } else ""
  if (nzchar(reason)) {
    return(list(ess_mean = NA_real_, mcse_mean = NA_real_, tau = NA_real_,
                undefined_reason = reason))
  }

  # Rows are iterations and the sole column is one chain; splitting does not
  # create independently initialized chains.
  draws <- matrix(x, ncol = 1L)
  estimates <- list(
    ess_mean = as.numeric(posterior::ess_basic(draws, split = TRUE)),
    mcse_mean = as.numeric(posterior::mcse_mean(draws))
  )
  reasons <- character()
  if (!is.finite(estimates$ess_mean) || estimates$ess_mean <= 0) {
    reasons <- c(reasons, paste("undefined mean ESS:", estimates$ess_mean))
    estimates$ess_mean <- NA_real_
  }
  if (!is.finite(estimates$mcse_mean) || estimates$mcse_mean <= 0) {
    reasons <- c(reasons, paste("undefined mean MCSE:", estimates$mcse_mean))
    estimates$mcse_mean <- NA_real_
  }
  c(estimates, list(
    tau = length(x) / estimates$ess_mean,
    undefined_reason = paste(reasons, collapse = "; ")
  ))
}

finite_median <- function(x) {
  stats::median(x[is.finite(x)])
}

# Keep each sampled test paired with its production statistics and diagnostics.
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
  coordinates <- colnames(draws)
  do.call(rbind, lapply(seq_len(ncol(draws)), function(j) {
    estimate <- mean_diagnostics(draws[, j], paste(method, coordinates[j]))
    data.frame(method = method, coordinate = coordinates[j], n = nrow(draws),
               ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean,
               tau = estimate$tau, undefined_reason = estimate$undefined_reason)
  }))
}

# Count consecutive identical full vectors, including the first state of each run.
maximum_holding_run <- function(draws) {
  n <- nrow(draws)
  same_as_previous <- rowSums(
    draws[-1L, , drop = FALSE] != draws[-n, , drop = FALSE]
  ) == 0L
  run_starts <- c(1L, which(!same_as_previous) + 1L)
  as.integer(max(diff(c(run_starts, n + 1L))))
}

# Minimum raw-coordinate ESS per selector budget, including warm-up calls.
recompute_efficiency <- function(result, validation) {
  methods <- c("Crude PMMH", "Correlated PMMH")
  rows <- match(methods, result$estimates$method)
  if (anyNA(rows) || anyDuplicated(result$estimates$method)) {
    stop("Efficiency requires one estimate for each PM and CPM method.")
  }
  coordinate <- do.call(rbind, lapply(methods, function(method) {
    coordinate_mean_diagnostics(validation$samples[[method]], method)
  }))
  for (i in seq_along(methods)) {
    method <- methods[i]
    selected <- rows[i]
    minimum <- min(coordinate$ess_mean[coordinate$method == method])
    info <- validation$sample_info[[method]]
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
  if (anyNA(rows$method) || anyDuplicated(result$estimates$method)) {
    stop("The efficiency table requires one estimate for each method.")
  }
  data.frame(
    method_id = ids, Method = c("Rejection", "PM", "CPM"),
    Efficiency = rows$efficiency_per_1000,
    Gain = rows$efficiency_per_1000 / rows$efficiency_per_1000[1L],
    Acceptance = rows$production_acceptance, Mismatch = rows$zero_estimate_rate,
    h_final = rows$final_proposal_scale
  )
}

# Rebuild each test indicator from its named trajectory and statistic column.
recompute_inference_diagnostics <- function(bank) {
  rows <- lapply(seq_len(nrow(bank$tests)), function(i) {
    test <- bank$tests[i, ]
    trajectory <- bank$trajectories[[test$trajectory_key]]
    series <- trajectory[, test$statistic_id]
    if (length(series) != bank$config$production ||
        any(!is.finite(series)) || !is.finite(test$statistic_observed)) {
      stop("Invalid production statistic for ", test$trajectory_key)
    }
    indicator <- as.numeric(series >= test$statistic_observed)
    p_value <- mean(indicator)
    estimate <- mean_diagnostics(indicator, paste(test$trajectory_key, test$statistic_id))
    cbind(test, data.frame(
      production_draws = length(indicator), tail_count = sum(indicator),
      p_value = p_value, reject = p_value <= bank$config$test_level,
      ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean,
      tau = estimate$tau, undefined_reason = estimate$undefined_reason
    ))
  })
  do.call(rbind, rows)
}

# Pool within application, experiment and method; omit undefined ESS/MCSE separately.
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

