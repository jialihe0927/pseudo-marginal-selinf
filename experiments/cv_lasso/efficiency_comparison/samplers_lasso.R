run_pmmh_lasso <- function(benchmark, config, method) {
  if (!method %in% c("Crude PMMH", "Correlated PMMH")) {
    stop("Unknown PMMH method: ", method)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required.")
  }

  rho <- if (method == "Crude PMMH") 0 else config$rho_correlated
  seed <- config$pmmh_seeds[[method]]
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", as.integer(seed), envir = .GlobalEnv)

  dimension <- length(benchmark$T_obs)
  proposal_factor <- t(chol(benchmark$Sigma_T))
  covariance_inverse <- solve(benchmark$Sigma_T)
  log_determinant <- 2 * sum(log(diag(chol(benchmark$Sigma_T))))
  state <- list(
    target = as.numeric(benchmark$T_obs),
    auxiliary = as.numeric(benchmark$W_obs / benchmark$gamma)
  )
  if (!selection_event_lasso(state$target, state$auxiliary, benchmark)) {
    stop("The exact observed initializer is not selected.")
  }

  scale <- config$proposal_scale_initial
  draws <- matrix(NA_real_, config$total_transitions, dimension)
  accepted <- logical(config$total_transitions)
  proposed_selected <- logical(config$total_transitions)

  for (iteration in seq_len(config$total_transitions)) {
    # Sigma_T fixes the proposal shape; scale^2 / d sets its magnitude.
    proposed_target <- as.numeric(
      state$target + scale / sqrt(dimension) * proposal_factor %*% rnorm(dimension)
    )
    noise <- rnorm(length(state$auxiliary))
    proposed_auxiliary <- if (rho == 0) {
      noise
    } else {
      rho * state$auxiliary + sqrt(1 - rho^2) * noise
    }
    selected <- selection_event_lasso(
      proposed_target,
      proposed_auxiliary,
      benchmark
    )
    proposed_selected[iteration] <- selected
    log_ratio <- if (selected) {
      log_gaussian_lasso(
        proposed_target,
        benchmark$gaussian_center,
        covariance_inverse,
        log_determinant
      ) - log_gaussian_lasso(
        state$target,
        benchmark$gaussian_center,
        covariance_inverse,
        log_determinant
      )
    } else {
      -Inf
    }
    acceptance_probability <- if (is.finite(log_ratio)) {
      min(1, exp(min(0, log_ratio)))
    } else {
      0
    }
    accepted[iteration] <- selected && log(runif(1L)) < min(0, log_ratio)

    # A rejected PMMH proposal retains both target and auxiliary states.
    if (accepted[iteration]) {
      state <- list(target = proposed_target, auxiliary = proposed_auxiliary)
    }
    draws[iteration, ] <- state$target

    if (iteration <= config$adaptive_transitions) {
      rate <- 0.5 / (iteration + 10)^0.6
      scale <- exp(
        log(scale) +
          rate * (acceptance_probability - config$target_acceptance)
      )
      scale <- min(
        config$proposal_scale_max,
        max(config$proposal_scale_min, scale)
      )
    }
  }

  production_start <- config$adaptive_transitions + 1L
  production <- seq.int(production_start, config$total_transitions)
  list(
    method = method,
    rho = rho,
    draws = draws[production, , drop = FALSE],
    accepted = accepted[production],
    proposed_selected = proposed_selected[production],
    selector_evaluations = config$total_transitions,
    final_proposal_scale = scale
  )
}

run_rejection_lasso <- function(benchmark, config) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(config$rejection_seed)
  streams <- vector("list", length(config$rejection_chunk_sizes))
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  if (length(streams) > 1L) {
    for (index in 2:length(streams)) {
      streams[[index]] <- parallel::nextRNGStream(streams[[index - 1L]])
    }
  }
  if (!is.null(config$rejection_rng_states)) {
    if (length(config$rejection_rng_states) != length(streams)) {
      stop("The saved rejection streams do not match the chunk budget.")
    }
    streams <- config$rejection_rng_states
  }
  dimension <- length(benchmark$T_obs)
  factor <- t(chol(benchmark$Sigma_T))
  hits <- 0L
  accepted <- matrix(numeric(0), ncol = dimension)

  for (chunk in seq_along(streams)) {
    assign(".Random.seed", streams[[chunk]], envir = .GlobalEnv)
    for (iteration in seq_len(config$rejection_chunk_sizes[chunk])) {
      target <- as.numeric(
        benchmark$gaussian_center + factor %*% rnorm(dimension)
      )
      auxiliary <- rnorm(length(benchmark$W_obs))
      if (selection_event_lasso(target, auxiliary, benchmark)) {
        hits <- hits + 1L
        accepted <- rbind(accepted, target)
      }
    }
  }
  list(hits = hits, attempts = config$rejection_attempts, accepted = accepted)
}

