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

# Build the bounded segmentation path, then choose K by penalized RSS.
# Ties in the criterion select the smaller K.
select_changepoints_cp <- function(
    series,
    minimum_segment = 10L,
    maximum_changepoints = 4L,
    penalty_coefficient = 1.25,
    tie_tolerance = 1e-12) {
  series <- as.numeric(series)
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

# Match the final unordered set after choosing K; do not condition on the path.
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

# Vary the contrast statistic while keeping the observed residual R_obs fixed.
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
