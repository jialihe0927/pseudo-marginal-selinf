log_target_nonnull_lasso <- function(T, covariance) {
  factor <- chol(covariance)
  standardized <- backsolve(factor, as.numeric(T), transpose = TRUE)
  -0.5 * (length(T) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(standardized^2))
}

run_pmmh_nonnull_lasso <- function(dataset, method_id, rng_state, config) {
  if (!method_id %in% config$sampled_methods) stop("Unknown Lasso PMMH method: ", method_id)
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
      proposed_T <- as.numeric(T + (scale / sqrt(dataset$d)) * proposal_factor %*% rnorm(dataset$d))
      proposed_U <- as.numeric(pmmh$rho * U + sqrt(1 - pmmh$rho^2) * rnorm(config$n))
      selected <- selection_event_nonnull_lasso(proposed_T, proposed_U, dataset, method_id, config)
      log_ratio <- if (selected) {
        log_target_nonnull_lasso(proposed_T, target$Sigma_M) -
          log_target_nonnull_lasso(T, target$Sigma_M)
      } else -Inf
      probability <- if (is.finite(log_ratio)) min(1, exp(min(0, log_ratio))) else 0
      # The formal function lazily evaluated its uniform only for a positive event.
      accepted <- selected && is.finite(log_ratio) && log(runif(1L)) < min(0, log_ratio)
      if (accepted) {
        T <- proposed_T
        U <- proposed_U
      }
      if (iteration <= pmmh$warmup) {
        next_scale <- exp(log(scale) + 0.5 / (iteration + 10)^0.6 *
                            (probability - pmmh$target_acceptance))
        scale <- min(pmmh$scale_max, max(pmmh$scale_min, next_scale))
      } else {
        production[iteration - pmmh$warmup, ] <- T
      }
    }
    # Retain the matrix multiplication used in the formal inclusive-tail comparison.
    scores <- statistic_nonnull_lasso(production, dataset$direction)
    observed <- statistic_nonnull_lasso(target$T_obs, dataset$direction)
    tail_count <- sum(scores >= observed)
    result_row_nonnull_lasso(
      dataset, method_id, tail_count / pmmh$production, config,
      tail_count = as.integer(tail_count), production_draws = pmmh$production,
      result_type = "simulation", selector_evaluations = 1L + pmmh$total_transitions
    )
  }, state = rng_state)
}
