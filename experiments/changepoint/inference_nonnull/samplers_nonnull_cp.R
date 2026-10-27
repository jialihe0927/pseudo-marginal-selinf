run_pmmh_nonnull_cp <- function(dataset, method_id, rng_state, config) {
  if (!method_id %in% config$sampled_methods) stop("Unknown changepoint PMMH method: ", method_id)
  if (dataset$direction$structural_zero) {
    return(result_row_nonnull_cp(dataset, method_id, 1, config, result_type = "structural_zero"))
  }
  with_rng_nonnull_cp(function() {
    pmmh <- config$pmmh
    target <- dataset$target
    T <- as.numeric(target$T_obs)
    U <- as.numeric(dataset$W / config$randomization_sd)
    # The frozen runner left-multiplies by the upper factor, without transposing.
    # Its symmetric random-walk covariance is (h^2/K) R R', where R'R = Sigma_T.
    proposal_factor <- chol(target$Sigma_T)
    scale <- pmmh$scale_initial
    production <- matrix(NA_real_, pmmh$production, dataset$selected_k)
    for (iteration in seq_len(pmmh$total_transitions)) {
      proposed_T <- drop(T + (scale / sqrt(dataset$selected_k)) *
                           proposal_factor %*% rnorm(dataset$selected_k))
      proposed_U <- pmmh$rho * U + sqrt(1 - pmmh$rho^2) * rnorm(config$n)
      selected <- selection_event_nonnull_cp(proposed_T, proposed_U, dataset, method_id, config)
      probability <- 0
      if (selected) {
        log_ratio <- -0.5 * drop(crossprod(proposed_T, target$Sigma_T_inv %*% proposed_T)) -
          (-0.5 * drop(crossprod(T, target$Sigma_T_inv %*% T)))
        probability <- min(1, exp(min(0, log_ratio)))
      }
      accepted <- runif(1L) < probability
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
    scores <- statistic_nonnull_cp(production, dataset$direction)
    tail_count <- sum(scores >= dataset$direction$observed)
    result_row_nonnull_cp(
      dataset, method_id, tail_count / pmmh$production, config,
      tail_count = as.integer(tail_count), production_draws = pmmh$production,
      result_type = "simulation", selector_evaluations = pmmh$total_transitions
    )
  }, state = rng_state)
}
