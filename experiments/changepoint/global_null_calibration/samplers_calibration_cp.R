adaptation_rate_null_cp <- function(iteration, config) {
  config$pmmh$adaptation_rate /
    (iteration + config$pmmh$adaptation_offset)^config$pmmh$adaptation_power
}

# Sample the zero-mean selective Gaussian law, retaining every production statistic.
# A supplied statistic function lets the directional delta-zero test use this same kernel.
run_pmmh_null_cp <- function(dataset, method_id, rng_state, config,
                              statistic_function = NULL, statistic_observed = NULL) {
  old_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", as.integer(rng_state), envir = .GlobalEnv)

  pmmh <- config$pmmh
  target <- dataset$target
  dimension <- dataset$selected_k
  target_value <- as.numeric(target$T_obs)
  auxiliary <- as.numeric(target$U_obs)

  if (!selection_event_null_cp(
    target_value,
    auxiliary,
    dataset,
    method_id,
    config
  )) {
    stop("The observed augmented state does not satisfy the selection event.")
  }

  # chol returns an upper factor; its transpose L has L L^T = Sigma_T.
  # With dimension K, the proposal covariance is (h^2/K) Sigma_T.
  proposal_factor <- t(chol(target$Sigma_T))
  proposal_scale <- pmmh$proposal_scale_initial
  tail_count <- 0L
  production_statistics <- numeric(pmmh$production)
  if (is.null(statistic_observed)) statistic_observed <- target$Q_obs
  warmup_accepted <- 0L
  production_accepted <- 0L

  for (iteration in seq_len(pmmh$total_transitions)) {
    current_target <- target_value
    current_auxiliary <- auxiliary

    proposed_target <- as.numeric(
      current_target +
      proposal_scale / sqrt(dimension) * proposal_factor %*% rnorm(dimension)
    )
    proposed_auxiliary <-
      pmmh$rho * current_auxiliary +
      sqrt(1 - pmmh$rho^2) * rnorm(config$n)

    selected <- selection_event_null_cp(
      proposed_target,
      proposed_auxiliary,
      dataset,
      method_id,
      config
    )
    log_ratio <- if (selected) {
      -0.5 * target_statistic_null_cp(proposed_target, target) +
        0.5 * target_statistic_null_cp(current_target, target)
    } else {
      -Inf
    }
    acceptance_probability <- if (is.finite(log_ratio)) {
      min(1, exp(min(0, log_ratio)))
    } else {
      0
    }
    # This kernel consumes a uniform only when the acceptance probability is positive.
    accepted <- acceptance_probability > 0 && runif(1L) < acceptance_probability

    # A rejection keeps both the target and auxiliary state.
    if (accepted) {
      target_value <- proposed_target
      auxiliary <- proposed_auxiliary
    }

    # Adapt in warm-up only; record every subsequent state at the fixed scale.
    if (iteration <= pmmh$warmup) {
      warmup_accepted <- warmup_accepted + as.integer(accepted)
      next_scale <- exp(
        log(proposal_scale) +
        adaptation_rate_null_cp(iteration, config) *
          (acceptance_probability - pmmh$target_acceptance)
      )
      proposal_scale <- min(
        pmmh$proposal_scale_max,
        max(pmmh$proposal_scale_min, next_scale)
      )
    } else {
      production_accepted <- production_accepted + as.integer(accepted)
      statistic <- if (is.null(statistic_function)) {
        target_statistic_null_cp(target_value, target)
      } else statistic_function(target_value)
      production_statistics[iteration - pmmh$warmup] <- statistic
      tail_count <- tail_count + as.integer(statistic >= statistic_observed)
    }
  }

  # Inclusive upper-tail count over the production length, without a +1 correction.
  p_value <- tail_count / pmmh$production
  row <- data.frame(
    dataset_id = dataset$dataset_id,
    raw_id = dataset$raw_id,
    target_id = "selected_changepoints",
    method_id = method_id,
    target_dimension = dimension,
    statistic_observed = statistic_observed,
    tail_count = tail_count,
    production_draws = pmmh$production,
    p_value = p_value,
    reject = as.integer(p_value <= config$test_level),
    result_type = "simulation",
    status = "complete",
    warmup_acceptance = warmup_accepted / pmmh$warmup,
    production_acceptance = production_accepted / pmmh$production,
    final_proposal_scale = proposal_scale,
    # Include the observed-state check as well as each transition selector call.
    selector_evaluations = 1L + pmmh$total_transitions,
    stringsAsFactors = FALSE
  )
  attr(row, "production_statistics") <- production_statistics
  attr(row, "indicator_diagnostics") <- mean_diagnostics(
    as.numeric(production_statistics >= statistic_observed),
    paste(method_id, dataset$dataset_id, "quadratic or directional tail")
  )
  row
}
