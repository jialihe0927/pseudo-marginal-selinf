selection_event <- function(x, u, a, b) {
  sqrt(a) * x + sqrt(1 - a) * u > b
}

log_selection_probability <- function(b) {
  pnorm(b, lower.tail = FALSE, log.p = TRUE)
}

draw_selected_states <- function(n, a, b) {
  if (n < 1L || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (!is.finite(a) || a <= 0 || a >= 1) {
    stop("a must lie strictly between zero and one.")
  }

  tail_uniform <- pmax(runif(n), .Machine$double.xmin)
  selected_score <- qnorm(
    log(tail_uniform) + log_selection_probability(b),
    lower.tail = FALSE,
    log.p = TRUE
  )
  orthogonal_score <- rnorm(n)
  x <- sqrt(a) * selected_score + sqrt(1 - a) * orthogonal_score
  u <- sqrt(1 - a) * selected_score - sqrt(a) * orthogonal_score

  if (any(!selection_event(x, u, a, b))) {
    stop("Exact selected-state initialization failed.")
  }
  list(x = x, u = u)
}

independent_pm_transition <- function(state, a, b, proposal_scale) {
  n_chains <- length(state$x)
  proposal_x <- state$x + proposal_scale * rnorm(n_chains)
  log_mh_ratio <- -0.5 * (proposal_x^2 - state$x^2)
  base_pass <- log(runif(n_chains)) <= pmin(0, log_mh_ratio)
  proposal_u <- rnorm(n_chains)
  event_pass <- selection_event(proposal_x, proposal_u, a, b)
  move <- base_pass & event_pass

  state$x <- ifelse(move, proposal_x, state$x)
  state$u <- ifelse(move, proposal_u, state$u)
  state
}

correlated_pm_transition <- function(state, a, b, proposal_scale, rho) {
  n_chains <- length(state$x)
  proposal_x <- state$x + proposal_scale * rnorm(n_chains)
  log_mh_ratio <- -0.5 * (proposal_x^2 - state$x^2)
  base_pass <- log(runif(n_chains)) <= pmin(0, log_mh_ratio)
  candidate_x <- ifelse(base_pass, proposal_x, state$x)

  proposal_u <- rho * state$u + sqrt((1 - rho) * (1 + rho)) * rnorm(n_chains)
  event_pass <- selection_event(candidate_x, proposal_u, a, b)

  # A base-MH rejection can still refresh the auxiliary variable. If the
  # selector rejects, both coordinates retain their previous values.
  state$x <- ifelse(event_pass, candidate_x, state$x)
  state$u <- ifelse(event_pass, proposal_u, state$u)
  state
}

initialize_chains <- function(n_chains, a, b, median_x, seed) {
  set.seed(seed)
  state <- draw_selected_states(n_chains, a, b)
  state$indicator_sum <- as.integer(state$x <= median_x)
  state$n_observations <- 1L
  state
}

advance_chains <- function(state, method, a, b, proposal_scale, rho,
                           median_x, target_observations) {
  method <- match.arg(method, c("Independent PM", "CPM"))
  if (target_observations < state$n_observations) {
    stop("target_observations cannot move backward.")
  }

  while (state$n_observations < target_observations) {
    if (method == "Independent PM") {
      state <- independent_pm_transition(state, a, b, proposal_scale)
    } else if (method == "CPM") {
      state <- correlated_pm_transition(state, a, b, proposal_scale, rho)
    } else {
      stop("Unknown sampling method: ", method)
    }

    if (any(!is.finite(state$x)) || any(!is.finite(state$u))) {
      stop("A sampler produced a non-finite state.")
    }
    if (any(!selection_event(state$x, state$u, a, b))) {
      stop("A sampler left the selected state space.")
    }

    state$indicator_sum <- state$indicator_sum + as.integer(state$x <= median_x)
    state$n_observations <- state$n_observations + 1L
  }
  state
}
