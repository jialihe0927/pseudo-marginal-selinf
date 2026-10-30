# Configuration.
validate_efficiency_config_cp <- function(config) {
  chunk_sizes <- config$rejection_chunk_sizes
  valid_chunks <- is.numeric(chunk_sizes) &&
    length(chunk_sizes) > 0L &&
    all(is.finite(chunk_sizes)) &&
    all(chunk_sizes > 0) &&
    all(chunk_sizes == as.integer(chunk_sizes))
  if (!valid_chunks) {
    stop("rejection_chunk_sizes must contain positive integers.")
  }
  if (sum(chunk_sizes) != config$rejection_attempts) {
    stop("rejection_chunk_sizes must sum to rejection_attempts.")
  }
  if (config$rejection_attempts != 3000000L) {
    stop("The formal rejection budget must be 3,000,000 attempts.")
  }

  transition_total <- config$adaptive_transitions +
    config$stabilization_transitions +
    config$production_transitions
  if (transition_total != config$total_transitions) {
    stop("Adaptive, stabilization, and production transitions must sum to total_transitions.")
  }
  if (config$stabilization_transitions != 0L) {
    stop("This implementation supports stabilization_transitions = 0 only.")
  }
  if (config$auxiliary_samples != 1L) {
    stop("The Bernoulli likelihood implementation requires auxiliary_samples = 1.")
  }
  if (config$production_transitions <= 0L) {
    stop("production_transitions must be positive.")
  }

  scales <- c(
    config$proposal_scale_min,
    config$proposal_scale_initial,
    config$proposal_scale_max
  )
  if (any(!is.finite(scales)) ||
      config$proposal_scale_min <= 0 ||
      config$proposal_scale_min > config$proposal_scale_initial ||
      config$proposal_scale_initial > config$proposal_scale_max) {
    stop("Proposal scale bounds must satisfy 0 < min <= initial <= max.")
  }
  if (!is.finite(config$rho_correlated) ||
      config$rho_correlated < 0 ||
      config$rho_correlated >= 1) {
    stop("rho_correlated must lie in [0, 1).")
  }

  invisible(config)
}

# Bounded segmentation and the final-set event.
segments_cp <- function(n, changepoints) {
  boundaries <- c(0L, sort(as.integer(changepoints)), as.integer(n))
  data.frame(
    start = head(boundaries, -1L) + 1L,
    end = tail(boundaries, -1L)
  )
}

cusum_cp <- function(series, start, split, end) {
  n_left <- split - start + 1L
  n_right <- end - split
  sqrt(n_left * n_right / (n_left + n_right)) *
    (mean(series[(split + 1L):end]) - mean(series[start:split]))
}

candidate_splits_cp <- function(series, segments, minimum_segment) {
  rows <- vector("list", nrow(segments))
  for (i in seq_len(nrow(segments))) {
    start <- segments$start[i]
    end <- segments$end[i]
    if (end - start + 1L < 2L * minimum_segment) next
    splits <- seq.int(start + minimum_segment - 1L, end - minimum_segment)
    signed <- vapply(
      splits,
      function(split) cusum_cp(series, start, split, end),
      numeric(1L)
    )
    rows[[i]] <- data.frame(
      start = start,
      split = splits,
      end = end,
      signed = signed,
      absolute = abs(signed)
    )
  }
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  candidates <- do.call(rbind, rows)
  ordering <- order(
    -candidates$absolute,
    candidates$start,
    candidates$split,
    candidates$end,
    method = "radix"
  )
  candidates[ordering, , drop = FALSE]
}

piecewise_rss_cp <- function(series, changepoints) {
  segments <- segments_cp(length(series), changepoints)
  sum(vapply(seq_len(nrow(segments)), function(i) {
    values <- series[segments$start[i]:segments$end[i]]
    sum((values - mean(values))^2)
  }, numeric(1L)))
}

select_changepoints_cp <- function(
    series,
    minimum_segment = 10L,
    maximum_changepoints = 4L,
    penalty_coefficient = 1.25,
    tie_tolerance = 1e-12) {
  series <- as.numeric(series)
  if (any(!is.finite(series)) ||
      length(series) < 2L * minimum_segment) {
    stop("Invalid series.")
  }
  path <- integer(0)
  rss <- numeric(maximum_changepoints + 1L)
  rss[1L] <- piecewise_rss_cp(series, path)
  for (step in seq_len(maximum_changepoints)) {
    current_segments <- segments_cp(length(series), path)
    candidates <- candidate_splits_cp(
      series,
      current_segments,
      minimum_segment
    )
    if (is.null(candidates)) stop("No legal split remains.")
    path <- c(path, as.integer(candidates$split[1L]))
    rss[step + 1L] <- piecewise_rss_cp(series, path)
  }
  penalty <- penalty_coefficient *
    (0:maximum_changepoints) *
    log(length(series))
  criterion <- rss + penalty
  minimizers <- which(abs(criterion - min(criterion)) <= tie_tolerance)
  selected_k <- as.integer(min(minimizers) - 1L)
  selected <- if (selected_k == 0L) integer(0) else sort(path[seq_len(selected_k)])
  list(
    selected_set = selected,
    selected_k = selected_k,
    ordered_path = path,
    criterion = criterion
  )
}

selection_event_cp <- function(target, auxiliary, benchmark) {
  series <- reconstruct_series_cp(target, benchmark) +
    benchmark$tau * as.numeric(auxiliary)
  selection <- select_changepoints_cp(
    series,
    benchmark$minimum_segment,
    benchmark$maximum_changepoints,
    benchmark$penalty_coefficient,
    benchmark$tie_tolerance
  )
  selected_set <- sort(as.integer(selection$selected_set))
  identical(selected_set, benchmark$selected_changepoints)
}

# Frozen selective Gaussian benchmark.
reconstruct_series_cp <- function(target, benchmark) {
  target_coefficients <- benchmark$Sigma_T_inv %*% as.numeric(target)
  as.numeric(benchmark$R_obs + benchmark$V %*% target_coefficients)
}

log_gaussian_cp <- function(target, benchmark) {
  delta <- as.numeric(target) - benchmark$gaussian_center
  -0.5 * (
    length(delta) * log(2 * pi) +
      benchmark$log_determinant +
      drop(crossprod(delta, benchmark$Sigma_T_inv %*% delta))
  )
}

# Mean ESS on each retained coordinate.
mean_ess_cp <- function(x, method, coordinate) {
  label <- paste(method, coordinate, sep = ", ")
  estimate <- mean_diagnostics(x, label)
  if (nzchar(estimate$undefined_reason)) message(label, ": ", estimate$undefined_reason)
  estimate$ess_mean
}
