run_pmmh_cp <- function(benchmark, config, method) {
  # PM refreshes the auxiliary independently; CPM correlates it with the current state.
  rho <- if (method == "Crude PMMH") 0 else config$rho_correlated
  seed <- config$pmmh_seeds[[method]]
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  chain_state <- get(".Random.seed", envir = .GlobalEnv)
  auxiliary_seed <- if (is.null(config$pmmh_auxiliary_seeds)) {
    seed + 1000L
  } else config$pmmh_auxiliary_seeds[[method]]
  proposal_seed <- if (is.null(config$pmmh_proposal_seeds)) {
    seed + 2000L
  } else config$pmmh_proposal_seeds[[method]]
  set.seed(auxiliary_seed)
  auxiliary_state <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(proposal_seed)
  proposal_state <- get(".Random.seed", envir = .GlobalEnv)

  # Keep proposal, auxiliary and acceptance draws on their separate streams.
  draw_with_state <- function(state, expression) {
    assign(".Random.seed", state, envir = .GlobalEnv)
    value <- force(expression)
    list(value = value, state = get(".Random.seed", envir = .GlobalEnv))
  }

  dimension <- length(benchmark$T_obs)
  state <- list(target = benchmark$T_obs, auxiliary = benchmark$U_obs)
  scale <- config$proposal_scale_initial
  draws <- matrix(NA_real_, config$total_transitions, dimension)
  accepted <- logical(config$total_transitions)
  proposed_selected <- logical(config$total_transitions)

  for (iteration in seq_len(config$total_transitions)) {
    # Sigma_T fixes the proposal shape; scale^2 / d sets its magnitude.
    proposal_draw <- draw_with_state(proposal_state, rnorm(dimension))
    proposal_state <- proposal_draw$state
    proposed_target <- as.numeric(
      state$target +
        scale / sqrt(dimension) *
          benchmark$proposal_factor %*% proposal_draw$value
    )
    auxiliary_draw <- draw_with_state(auxiliary_state, rnorm(length(state$auxiliary)))
    auxiliary_state <- auxiliary_draw$state
    proposed_auxiliary <- if (rho == 0) {
      auxiliary_draw$value
    } else {
      rho * state$auxiliary +
        sqrt(1 - rho^2) * auxiliary_draw$value
    }
    selected <- selection_event_cp(proposed_target, proposed_auxiliary, benchmark)
    proposed_selected[iteration] <- selected
    log_ratio <- if (selected) {
      log_gaussian_cp(proposed_target, benchmark) -
        log_gaussian_cp(state$target, benchmark)
    } else {
      -Inf
    }
    acceptance_probability <- if (is.finite(log_ratio)) {
      min(1, exp(min(0, log_ratio)))
    } else {
      0
    }
    uniform_draw <- draw_with_state(chain_state, runif(1L))
    chain_state <- uniform_draw$state
    accepted[iteration] <- selected && log(uniform_draw$value) < min(0, log_ratio)

    # Rejection preserves both components of the extended state.
    if (accepted[iteration]) {
      state <- list(
        target = proposed_target,
        auxiliary = proposed_auxiliary
      )
    }
    draws[iteration, ] <- state$target
    # Adapt only during warm-up; the production scale stays fixed.
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

  # Remove warm-up once and retain all production states, including repeats.
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

# Draw independent Gaussian candidates at the same center and event as PM/CPM.
run_rejection_cp <- function(benchmark, config) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(config$rejection_seed)
  streams <- config$rejection_rng_states
  if (is.null(streams)) {
    streams <- vector("list", length(config$rejection_chunk_sizes))
    streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
    if (length(streams) > 1L) {
      for (index in 2:length(streams)) {
        streams[[index]] <- parallel::nextRNGStream(streams[[index - 1L]])
      }
    }
  }
  hits <- 0L
  accepted <- matrix(numeric(0), ncol = length(benchmark$T_obs))
  for (chunk in seq_along(streams)) {
    assign(".Random.seed", streams[[chunk]], envir = .GlobalEnv)
    for (iteration in seq_len(config$rejection_chunk_sizes[chunk])) {
      target <- as.numeric(
        benchmark$gaussian_center +
          benchmark$proposal_factor %*% rnorm(length(benchmark$T_obs))
      )
      auxiliary <- rnorm(length(benchmark$U_obs))
      if (selection_event_cp(target, auxiliary, benchmark)) {
        hits <- hits + 1L
        accepted <- rbind(accepted, target)
      }
    }
  }
  list(hits = hits, attempts = config$rejection_attempts, accepted = accepted)
}

