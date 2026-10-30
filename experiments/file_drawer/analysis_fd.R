posterior_cdf <- function(x, a, b) {
  integrand <- function(y) {
    log_likelihood <- pnorm(
      (b - sqrt(a) * y) / sqrt(1 - a),
      lower.tail = FALSE,
      log.p = TRUE
    )
    exp(dnorm(y, log = TRUE) + log_likelihood - log_selection_probability(b))
  }
  integrate(
    integrand,
    lower = -Inf,
    upper = x,
    rel.tol = 1e-9,
    abs.tol = 0,
    subdivisions = 1000L
  )$value
}

posterior_median <- function(a, b) {
  inverse_mills <- exp(dnorm(b, log = TRUE) - log_selection_probability(b))
  selected_variance <- 1 + inverse_mills * (b - inverse_mills)
  center <- sqrt(a) * inverse_mills
  half_width <- max(8 * sqrt(a * selected_variance + 1 - a), 2)
  objective <- function(x) posterior_cdf(x, a, b) - 0.5

  lower <- center - half_width
  upper <- center + half_width
  for (attempt in seq_len(20L)) {
    if (objective(lower) < 0 && objective(upper) > 0) {
      return(uniroot(objective, c(lower, upper), tol = 1e-9)$root)
    }
    half_width <- 2 * half_width
    lower <- center - half_width
    upper <- center + half_width
  }
  stop("Could not bracket the selected-distribution median.")
}

checkpoint_grid <- function(maximum_observations, log_step = log(1.05)) {
  if (length(maximum_observations) != 1L ||
      !is.numeric(maximum_observations) || !is.finite(maximum_observations) ||
      maximum_observations < 2 || maximum_observations != floor(maximum_observations) ||
      maximum_observations > .Machine$integer.max) {
    stop("maximum_observations must be an integer from 2 to .Machine$integer.max.")
  }
  if (length(log_step) != 1L || !is.numeric(log_step) ||
      !is.finite(log_step) || log_step <= 0) {
    stop("checkpoint_log_step must be finite and positive.")
  }
  sort(unique(c(
    1L,
    as.integer(round(exp(seq(
      log(1),
      log(maximum_observations),
      by = log_step
    )))),
    as.integer(maximum_observations)
  )))
}

stable_crossing_index <- function(pass, required = 8L) {
  runs <- rle(as.logical(pass))
  starts <- cumsum(c(1L, head(runs$lengths, -1L)))
  eligible <- which(runs$values & runs$lengths >= required)
  if (length(eligible)) starts[eligible[1L]] else NA_integer_
}

estimate_fixed_accuracy <- function(method, a, b, rho, proposal_scale,
                                    maximum_observations, n_chains,
                                    epsilon, stable_checkpoints,
                                    post_certification_checkpoints,
                                    seed,
                                    checkpoint_log_step = log(1.05)) {
  method <- match.arg(method, c("Independent PM", "CPM"))
  if (length(n_chains) != 1L || !is.numeric(n_chains) ||
      !is.finite(n_chains) || n_chains < 2 || n_chains != floor(n_chains)) {
    stop("n_chains must be an integer of at least 2.")
  }
  if (length(stable_checkpoints) != 1L || !is.numeric(stable_checkpoints) ||
      !is.finite(stable_checkpoints) || stable_checkpoints < 1 ||
      stable_checkpoints != floor(stable_checkpoints)) {
    stop("stable_checkpoints must be a positive integer.")
  }
  if (length(post_certification_checkpoints) != 1L ||
      !is.numeric(post_certification_checkpoints) ||
      !is.finite(post_certification_checkpoints) ||
      post_certification_checkpoints < 0 ||
      post_certification_checkpoints != floor(post_certification_checkpoints)) {
    stop("post_certification_checkpoints must be a nonnegative integer.")
  }
  if (length(epsilon) != 1L || !is.numeric(epsilon) ||
      !is.finite(epsilon) || epsilon <= 0 || epsilon >= 1) {
    stop("epsilon must be finite and strictly between 0 and 1.")
  }

  # Validate the grid before the median calculation or random initialization.
  checkpoints <- checkpoint_grid(
    maximum_observations,
    log_step = checkpoint_log_step
  )
  median_x <- posterior_median(a, b)
  state <- initialize_chains(n_chains, a, b, median_x, seed)
  means <- matrix(NA_real_, nrow = length(checkpoints), ncol = n_chains)
  means[1L, ] <- state$indicator_sum
  crossing <- NA_integer_
  stop_index <- length(checkpoints)
  final_index <- 1L

  for (checkpoint_index in seq_along(checkpoints)[-1L]) {
    state <- advance_chains(
      state = state,
      method = method,
      a = a,
      b = b,
      proposal_scale = proposal_scale,
      rho = rho,
      median_x = median_x,
      target_observations = checkpoints[checkpoint_index]
    )
    means[checkpoint_index, ] <- state$indicator_sum / state$n_observations
    final_index <- checkpoint_index
    variance_ratio <- 4 * apply(
      means[seq_len(checkpoint_index), , drop = FALSE],
      1L,
      var
    )

    if (is.na(crossing)) {
      crossing <- stable_crossing_index(
        variance_ratio <= epsilon,
        stable_checkpoints
      )
      if (!is.na(crossing)) {
        certification <- crossing + stable_checkpoints - 1L
        stop_index <- min(
          length(checkpoints),
          certification + post_certification_checkpoints
        )
      }
    }
    if (!is.na(crossing) && checkpoint_index >= stop_index) break
  }

  means <- means[seq_len(final_index), , drop = FALSE]
  checkpoints <- checkpoints[seq_len(final_index)]
  variance_ratio <- 4 * apply(means, 1L, var)
  crossing <- stable_crossing_index(
    variance_ratio <= epsilon,
    stable_checkpoints
  )

  # complexity_n includes the initial observation and precedes certification.
  # The budget counts proposal event evaluations across all chains, excluding
  # exact-initialization validation and checkpoint integrity checks.
  estimate <- data.frame(
    b = b,
    method = method,
    a = a,
    rho = rho,
    proposal_scale = proposal_scale,
    complexity_n = if (is.na(crossing)) NA_real_ else checkpoints[crossing],
    observations_simulated = tail(checkpoints, 1L),
    transitions_per_chain = tail(checkpoints, 1L) - 1L,
    selector_evaluations = n_chains * (tail(checkpoints, 1L) - 1L),
    resolved = !is.na(crossing),
    estimate_type = "stable crossing start",
    result_type = "simulation",
    stringsAsFactors = FALSE
  )

  list(
    estimate = estimate,
    checkpoints = data.frame(
      observations = checkpoints,
      variance_ratio = variance_ratio
    )
  )
}

expm1_minus_x <- function(x) {
  ifelse(
    abs(x) < 1e-4,
    x^2 / 2 + x^3 / 6 + x^4 / 24 + x^5 / 120 + x^6 / 720,
    expm1(x) - x
  )
}

log1m_plus_probability <- function(probability) {
  ifelse(
    probability < 1e-4,
    -(probability^2 / 2 + probability^3 / 3 + probability^4 / 4 +
        probability^5 / 5 + probability^6 / 6),
    log1p(-probability) + probability
  )
}

rejection_variance_ratio <- function(n, b) {
  selection_probability <- exp(log_selection_probability(b))
  log_persistence <- log1p(-selection_probability)
  persistence <- exp(log_persistence)

  calculate_one <- function(n_value) {
    if (n_value == 1) return(1)
    scaled_log_persistence <- n_value * log_persistence
    numerator <- if (abs(scaled_log_persistence) < 0.1) {
      expm1_minus_x(scaled_log_persistence) +
        n_value * log1m_plus_probability(selection_probability)
    } else {
      n_value * selection_probability + expm1(scaled_log_persistence)
    }
    weighted_sum <- persistence * numerator / selection_probability^2
    (n_value + 2 * weighted_sum) / n_value^2
  }
  vapply(n, calculate_one, numeric(1L))
}

rejection_fixed_accuracy_time <- function(b, epsilon = 0.10) {
  lower <- 1
  upper <- 1
  while (rejection_variance_ratio(upper, b) > epsilon) {
    lower <- upper
    upper <- upper * 2
    if (!is.finite(upper) || upper > 1e100) {
      stop("The rejection-time search exceeded its numerical range.")
    }
  }

  repeat {
    resolution <- max(1, 8 * .Machine$double.eps * upper)
    if (upper - lower <= resolution) break
    middle <- floor((lower + upper) / 2)
    if (middle <= lower || middle >= upper) break
    if (rejection_variance_ratio(middle, b) <= epsilon) {
      upper <- middle
    } else {
      lower <- middle
    }
  }
  upper
}
