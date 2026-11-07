# The DGP may have a positive signal; the selective test reference remains zero-mean Gaussian.
run_pmmh_nonnull_cp <- function(dataset, method_id, rng_state, config) {
  # A zero projected direction has no standardized score; return p=1 without sampling.
  if (dataset$direction$structural_zero) {
    return(result_row_nonnull_cp(dataset, method_id, 1, config, result_type = "structural_zero"))
  }
  # Rerun the calibration kernel and assigned stream with a directional score.
  if (dataset$delta == 0) {
    null <- list(
      dataset_id = dataset$dataset_id,
      raw_id = dataset$raw_id,
      response = dataset$Y,
      randomization = dataset$W,
      selected_k = dataset$selected_k,
      selected_changepoints = dataset$selected_changepoints
    )
    null_config <- null_config_nonnull_cp(config)
    null$target <- build_target_null_cp(null, null_config)
    score <- function(T) statistic_nonnull_cp(T, dataset$direction)
    sampled <- run_pmmh_null_cp(
      null, method_id, rng_state, null_config, score, dataset$direction$observed
    )
    row <- result_row_nonnull_cp(
      dataset, method_id, sampled$p_value, config,
      sampled$tail_count, sampled$production_draws,
      "simulation", sampled$selector_evaluations
    )
    attr(row, "production_statistics") <- attr(sampled, "production_statistics")
    attr(row, "indicator_diagnostics") <- attr(sampled, "indicator_diagnostics")
    return(row)
  }
  with_rng_nonnull_cp(function() {
    pmmh <- config$pmmh
    target <- dataset$target
    T <- as.numeric(target$T_obs)
    U <- as.numeric(dataset$W / config$randomization_sd)
    # R chol() returns an upper factor; transposing gives L L^T = Sigma_T,
    # hence proposal covariance (h^2/K) Sigma_T.
    proposal_factor <- t(chol(target$Sigma_T))
    scale <- pmmh$scale_initial
    production <- matrix(NA_real_, pmmh$production, dataset$selected_k)
    for (iteration in seq_len(pmmh$total_transitions)) {
      proposed_T <- drop(T + (scale / sqrt(dataset$selected_k)) *
                           drop(proposal_factor %*% rnorm(dataset$selected_k)))
      proposed_U <- pmmh$rho * U + sqrt(1 - pmmh$rho^2) * rnorm(config$n)
      selected <- selection_event_nonnull_cp(
        proposed_T, proposed_U, dataset, method_id, config
      )
      probability <- 0
      if (selected) {
        log_ratio <-
          -0.5 * drop(crossprod(proposed_T, target$Sigma_T_inv %*% proposed_T)) -
          (-0.5 * drop(crossprod(T, target$Sigma_T_inv %*% T)))
        probability <- min(1, exp(min(0, log_ratio)))
      }
      # This kernel draws a uniform even when the selection estimate is zero.
      accepted <- runif(1L) < probability
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
    scores <- statistic_nonnull_cp(production, dataset$direction)
    # Inclusive tail count / production length, without a +1 correction.
    tail_count <- sum(scores >= dataset$direction$observed)
    row <- result_row_nonnull_cp(
      dataset, method_id, tail_count / pmmh$production, config,
      tail_count = as.integer(tail_count), production_draws = pmmh$production,
      result_type = "simulation", selector_evaluations = pmmh$total_transitions
    )
    attr(row, "production_statistics") <- scores
    attr(row, "indicator_diagnostics") <- mean_diagnostics(
      as.numeric(scores >= dataset$direction$observed), paste(method_id, dataset$dataset_id)
    )
    row
  }, state = rng_state)
}
