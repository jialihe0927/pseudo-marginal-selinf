log_target_nonnull_lasso <- function(T, covariance) {
  factor <- chol(covariance)
  standardized <- backsolve(factor, as.numeric(T), transpose = TRUE)
  -0.5 * (length(T) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(standardized^2))
}

run_pmmh_nonnull_lasso <- function(dataset, method_id, rng_state, config) {
  validate_config_nonnull_lasso(config)
  if (!method_id %in% config$sampled_methods) stop("Unknown Lasso PMMH method: ", method_id)
  # Proposed at delta zero reuses the quadratic-null kernel and its RNG order.
  if (dataset$delta == 0 && method_id == "model_only") {
    null <- list(
      dataset_id = dataset$dataset_id,
      raw_id = dataset$raw_id,
      design = dataset$X,
      response = dataset$Y,
      randomization = dataset$W,
      fold_id = dataset$fold_id,
      lambda_grid = config$lambda_grid,
      selected_model = dataset$selected_model,
      target_dimension = dataset$d,
      lambda_index = dataset$lambda_index,
      lambda_observed = dataset$lambda_observed
    )
    null_config <- null_config_nonnull_lasso(config)
    null$target <- build_target_null_lasso(null, null_config)
    score <- function(T) statistic_nonnull_lasso(T, dataset$direction)
    observed <- score(null$target$T_obs)
    sampled <- run_pmmh_null_lasso(null, method_id, rng_state, null_config, score, observed)
    row <- result_row_nonnull_lasso(
      dataset, method_id, sampled$p_value, config,
      sampled$tail_count, sampled$production_draws,
      "simulation", sampled$selector_evaluations
    )
    attr(row, "production_statistics") <- attr(sampled, "production_statistics")
    attr(row, "indicator_diagnostics") <- attr(sampled, "indicator_diagnostics")
    return(row)
  }
  with_rng_nonnull_lasso(function() {
    pmmh <- config$pmmh
    target <- dataset$target
    if (!selection_event_nonnull_lasso(
      target$T_obs, as.numeric(dataset$W / config$randomization_sd), dataset, method_id, config
    )) stop("The observed augmented initializer is infeasible.")
    T <- as.numeric(target$T_obs)
    U <- as.numeric(dataset$W / config$randomization_sd)
    proposal_factor <- t(chol(target$Sigma_M))
    scale <- pmmh$scale_initial
    production <- matrix(NA_real_, pmmh$production, dataset$d)
    for (iteration in seq_len(pmmh$total_transitions)) {
      proposed_T <- as.numeric(
        T + (scale / sqrt(dataset$d)) * proposal_factor %*% rnorm(dataset$d)
      )
      proposed_U <- as.numeric(
        pmmh$rho * U + sqrt(1 - pmmh$rho^2) * rnorm(config$n)
      )
      selected <- selection_event_nonnull_lasso(
        proposed_T, proposed_U, dataset, method_id, config
      )
      log_ratio <- if (selected) {
        log_target_nonnull_lasso(proposed_T, target$Sigma_M) -
          log_target_nonnull_lasso(T, target$Sigma_M)
      } else -Inf
      probability <- if (is.finite(log_ratio)) min(1, exp(min(0, log_ratio))) else 0
      # Draw the uniform only when the event and finite log ratio are admissible.
      accepted <- selected && is.finite(log_ratio) && log(runif(1L)) < min(0, log_ratio)
      if (accepted) {
        T <- proposed_T
        U <- proposed_U
      }
      # Rejection preserves both states; production retains every holding state.
      if (iteration <= pmmh$warmup) {
        next_scale <- exp(log(scale) + 0.5 / (iteration + 10)^0.6 *
                            (probability - pmmh$target_acceptance))
        scale <- min(pmmh$scale_max, max(pmmh$scale_min, next_scale))
      } else {
        production[iteration - pmmh$warmup, ] <- T
      }
    }
    # Keep this matrix multiplication path for the inclusive-tail comparison.
    scores <- statistic_nonnull_lasso(production, dataset$direction)
    observed <- statistic_nonnull_lasso(target$T_obs, dataset$direction)
    tail_count <- sum(scores >= observed)
    row <- result_row_nonnull_lasso(
      dataset, method_id, tail_count / pmmh$production, config,
      tail_count = as.integer(tail_count), production_draws = pmmh$production,
      result_type = "simulation", selector_evaluations = 1L + pmmh$total_transitions
    )
    attr(row, "production_statistics") <- scores
    attr(row, "indicator_diagnostics") <- mean_diagnostics(
      as.numeric(scores >= observed),
      paste(method_id, dataset$dataset_id)
    )
    row
  }, state = rng_state)
}
