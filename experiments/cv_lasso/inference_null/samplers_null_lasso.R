log_gaussian_target_null_lasso <- function(target_value, covariance) {
  factor <- chol(covariance)
  standardized <- backsolve(factor, as.numeric(target_value), transpose = TRUE)
  -0.5 * (
    length(target_value) * log(2 * pi) +
    2 * sum(log(diag(factor))) +
    sum(standardized^2)
  )
}

adaptation_rate_null_lasso <- function(iteration) {
  0.5 / (as.numeric(iteration) + 10)^0.6
}

run_pmmh_null_lasso <- function(dataset, method_id, rng_state, config,
                              statistic_function = NULL, statistic_observed = NULL) {
  validate_config_null_lasso(config)
  if (!method_id %in% c("model_only", "fixed_lambda")) {
    stop("Unknown CV-Lasso PMMH method: ", method_id)
  }

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

  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", as.integer(rng_state), envir = .GlobalEnv)

  pmmh <- config$pmmh
  target <- dataset$target
  dimension <- dataset$target_dimension
  target_value <- as.numeric(target$T_obs)
  auxiliary <- as.numeric(dataset$randomization / config$randomization_sd)

  if (!selection_event_null_lasso(
    target_value,
    auxiliary,
    dataset,
    method_id,
    config
  )) {
    stop("The observed augmented state does not satisfy the selection event.")
  }

  proposal_factor <- t(chol(target$Sigma_M))
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
    log_uniform <- log(runif(1L))

    selected <- selection_event_null_lasso(
      proposed_target,
      proposed_auxiliary,
      dataset,
      method_id,
      config
    )
    log_ratio <- if (selected) {
      log_gaussian_target_null_lasso(proposed_target, target$Sigma_M) -
        log_gaussian_target_null_lasso(current_target, target$Sigma_M)
    } else {
      -Inf
    }
    acceptance_probability <- if (is.finite(log_ratio)) {
      min(1, exp(min(0, log_ratio)))
    } else {
      0
    }
    accepted <- selected && log_uniform < min(0, log_ratio)

    if (accepted) {
      target_value <- proposed_target
      auxiliary <- proposed_auxiliary
    }

    if (iteration <= pmmh$warmup) {
      warmup_accepted <- warmup_accepted + as.integer(accepted)
      next_scale <- exp(
        log(proposal_scale) +
        adaptation_rate_null_lasso(iteration) *
          (acceptance_probability - pmmh$target_acceptance)
      )
      proposal_scale <- min(
        pmmh$proposal_scale_max,
        max(pmmh$proposal_scale_min, next_scale)
      )
    } else {
      production_accepted <- production_accepted + as.integer(accepted)
      statistic <- if (is.null(statistic_function)) {
        target_statistic_null_lasso(target_value, target)
      } else statistic_function(target_value)
      production_statistics[iteration - pmmh$warmup] <- statistic
      tail_count <- tail_count + as.integer(statistic >= statistic_observed)
    }
  }

  p_value <- tail_count / pmmh$production
  row <- data.frame(
    dataset_id = dataset$dataset_id,
    raw_id = dataset$raw_id,
    target_id = "selected_model",
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
