default_config_null_cp <- function() {
  list(
    experiment = "Changepoint global-null calibration", n = 200L,
    response_distribution = "N(0, I_n)", randomization_sd = 0.1,
    minimum_segment = 10L, maximum_changepoints = 4L,
    penalty_coefficient = 1.25, tie_tolerance = 1e-12,
    dataset_count = 200L, raw_attempt_limit = 10000L,
    retention_rule = "first eligible datasets in raw-attempt order",
    methods = data.frame(
      method_id = c("full_selector", "fixed_k", "naive_gaussian"),
      display_name = c("Full-selector", "Fixed-K", "Naive"),
      event = c(
        "adaptive final unordered changepoint set equals the observed set",
        "fix observed K and match the final unordered changepoint set",
        "no selection adjustment; analytic chi-square reference"
      ),
      stringsAsFactors = FALSE
    ),
    sampled_method_ids = c("full_selector", "fixed_k"),
    statistic = "Q = T' Sigma_T^{-1} T",
    p_value = "inclusive upper-tail count / 5000; no +1 correction",
    test_level = 0.05, rejection_rule = "p <= 0.05",
    pmmh = list(
      warmup = 200L, stabilization = 0L, production = 5000L,
      total_transitions = 5200L, auxiliary_samples = 1L, rho = 0.8,
      proposal_scale_initial = 0.8, proposal_scale_min = 0.05,
      proposal_scale_max = 5, target_acceptance = 0.20,
      adaptation_rate = 0.5, adaptation_offset = 10, adaptation_power = 0.6,
      production_adaptation = FALSE,
      proposal = "T' = T + h/sqrt(d) G xi, with G G' = Sigma_T"
    ),
    rng = list(
      kind = "L'Ecuyer-CMRG", raw_master_seed = 20261013L,
      formal_master_seed = 20261014L,
      assignment = "preallocated nextRNGStream state per raw_id/job; scheduling independent"
    )
  )
}

segments_null_cp <- function(n, changepoints) {
  boundaries <- c(0L, sort(as.integer(changepoints)), as.integer(n))
  data.frame(
    start = head(boundaries, -1L) + 1L,
    end = tail(boundaries, -1L)
  )
}

cusum_null_cp <- function(series, start, split, end) {
  left_length <- split - start + 1L
  right_length <- end - split
  sqrt(left_length * right_length / (left_length + right_length)) *
    (mean(series[(split + 1L):end]) - mean(series[start:split]))
}

candidate_splits_null_cp <- function(series, segments, minimum_segment) {
  rows <- vector("list", nrow(segments))
  for (index in seq_len(nrow(segments))) {
    start <- segments$start[index]
    end <- segments$end[index]
    if (end - start + 1L < 2L * minimum_segment) next

    splits <- seq.int(start + minimum_segment - 1L, end - minimum_segment)
    signed <- vapply(
      splits,
      function(split) cusum_null_cp(series, start, split, end),
      numeric(1L)
    )
    rows[[index]] <- data.frame(
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

piecewise_rss_null_cp <- function(series, changepoints) {
  segments <- segments_null_cp(length(series), changepoints)
  sum(vapply(seq_len(nrow(segments)), function(index) {
    values <- series[segments$start[index]:segments$end[index]]
    sum((values - mean(values))^2)
  }, numeric(1L)))
}

binary_segmentation_path_null_cp <- function(series, config) {
  path <- integer(0L)
  rss <- numeric(config$maximum_changepoints + 1L)
  rss[1L] <- piecewise_rss_null_cp(series, path)

  for (step in seq_len(config$maximum_changepoints)) {
    candidates <- candidate_splits_null_cp(
      series,
      segments_null_cp(length(series), path),
      config$minimum_segment
    )
    if (is.null(candidates)) stop("No legal split remains.")
    winner <- candidates[1L, , drop = FALSE]
    path <- c(path, as.integer(winner$split))
    rss[step + 1L] <- piecewise_rss_null_cp(series, path)
  }

  list(path = path, rss = rss)
}

# Choose K from the complete bounded path by penalized RSS, breaking ties toward smaller K.
select_adaptive_changepoints_null <- function(series, config) {
  path_result <- binary_segmentation_path_null_cp(series, config)
  candidate_k <- 0:config$maximum_changepoints
  penalty <- config$penalty_coefficient * candidate_k * log(length(series))
  criterion <- path_result$rss + penalty
  minimum <- min(criterion)
  minimizers <- which(abs(criterion - minimum) <= config$tie_tolerance)
  selected_k <- as.integer(min(minimizers) - 1L)
  selected <- if (selected_k == 0L) {
    integer(0L)
  } else {
    sort(path_result$path[seq_len(selected_k)])
  }

  list(
    selected_k = selected_k,
    selected = selected,
    ordered_path = path_result$path,
    criterion = criterion
  )
}

select_fixed_k_changepoints_null <- function(series, selected_k, config) {
  selected_k <- as.integer(selected_k)

  path <- integer(0L)
  for (step in seq_len(selected_k)) {
    candidates <- candidate_splits_null_cp(
      series,
      segments_null_cp(length(series), path),
      config$minimum_segment
    )
    if (is.null(candidates)) stop("No legal fixed-K split remains.")
    path <- c(path, as.integer(candidates$split[1L]))
  }
  list(selected_k = selected_k, selected = sort(path[seq_len(selected_k)]))
}

same_changepoints_null_cp <- function(left, right) {
  identical(sort(as.integer(left)), sort(as.integer(right)))
}

# Proposed reselects K; Fixed-K takes the observed number of splits.
# Both compare only the final unordered set, not the ordered search path.
selection_event_null_cp <- function(target_value, auxiliary, dataset, method_id, config) {
  selection_input <-
    reconstruct_response_null_cp(target_value, dataset$target) +
    config$randomization_sd * as.numeric(auxiliary)

  selected <- if (identical(method_id, "full_selector")) {
    select_adaptive_changepoints_null(selection_input, config)$selected
  } else if (identical(method_id, "fixed_k")) {
    select_fixed_k_changepoints_null(
      selection_input,
      dataset$selected_k,
      config
    )$selected
  } else {
    stop("Unknown changepoint method: ", method_id)
  }

  same_changepoints_null_cp(selected, dataset$selected_changepoints)
}

# Reconstruct Y at a new contrast T while keeping the nuisance residual R_obs fixed.
reconstruct_response_null_cp <- function(target_value, target) {
  as.numeric(
    target$R_obs +
    target$V %*% (target$Sigma_T_inv %*% as.numeric(target_value))
  )
}

# Quadratic null statistic Q = T^T Sigma_T^(-1) T.
target_statistic_null_cp <- function(target_value, target) {
  drop(crossprod(target_value, target$Sigma_T_inv %*% target_value))
}

# Adjacent right-minus-left segment means define the columns of V and T = V^T Y.
build_target_null_cp <- function(dataset, config) {
  selected <- sort(as.integer(dataset$selected_changepoints))
  selected_k <- as.integer(dataset$selected_k)
  segments <- segments_null_cp(config$n, selected)
  V <- matrix(0, config$n, selected_k)
  for (index in seq_len(selected_k)) {
    V[segments$start[index + 1L]:segments$end[index + 1L], index] <- 1 / (segments$end[index + 1L] - segments$start[index + 1L] + 1L)
    V[segments$start[index]:segments$end[index], index] <- -1 / (segments$end[index] - segments$start[index] + 1L)
  }
  Sigma_T <- crossprod(V)
  Sigma_T_inv <- solve(Sigma_T)
  T_obs <- drop(crossprod(V, dataset$response))
  R_obs <- drop(dataset$response - V %*% (Sigma_T_inv %*% T_obs))
  Q_obs <- drop(crossprod(T_obs, Sigma_T_inv %*% T_obs))
  list(
    dataset_id = dataset$dataset_id, raw_id = dataset$raw_id, K_obs = selected_k,
    M_obs = selected,
    segments = lapply(seq_len(nrow(segments)), function(index) {
      seq.int(segments$start[index], segments$end[index])
    }),
    V = V, Sigma_T = Sigma_T, Sigma_T_inv = Sigma_T_inv, T_obs = T_obs,
    R_obs = R_obs, Q_obs = Q_obs,
    p_naive = stats::pchisq(Q_obs, df = selected_k, lower.tail = FALSE),
    U_obs = dataset$randomization / config$randomization_sd
  )
}

# Generate null Gaussian responses and independent selector randomization.
generate_global_null_cp <- function(config) {
  response <- rnorm(config$n)
  randomization <- rnorm(config$n, sd = config$randomization_sd)
  list(
    response = response,
    randomization = randomization,
    selection_input = response + randomization
  )
}

rng_streams_null_cp <- function(master_seed, n_streams) {
  if (n_streams < 1L) return(list())
  old_kind <- RNGkind()
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(master_seed))
  streams <- vector("list", n_streams)
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  if (n_streams > 1L) for (index in 2:n_streams) streams[[index]] <- parallel::nextRNGStream(streams[[index - 1L]])
  lapply(streams, as.integer)
}

with_rng_state_null_cp <- function(state, fun) {
  old_kind <- RNGkind()
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", as.integer(state), envir = .GlobalEnv)
  value <- fun()
  list(value = value, final_state = as.integer(get(".Random.seed", envir = .GlobalEnv)))
}

generate_raw_attempt_null_cp <- function(raw_id, rng_state, config) {
  generated <- with_rng_state_null_cp(rng_state, function() generate_global_null_cp(config))
  observed <- select_adaptive_changepoints_null(generated$value$selection_input, config)
  c(
    list(raw_id = as.integer(raw_id), generation_stream_id = as.integer(raw_id),
         generation_rng_state = as.integer(rng_state)), generated$value,
    list(selected_k = observed$selected_k, selected_changepoints = observed$selected,
         ordered_path = observed$ordered_path, criterion = observed$criterion)
  )
}

# Assign one stream per retained dataset/method pair in bank order.
assign_rng_null_cp <- function(datasets, config) {
  sampled <- if (is.null(config$sampled_method_ids)) {
    setdiff(as.character(config$methods$method_id), "naive_gaussian")
  } else as.character(config$sampled_method_ids)
  streams <- rng_streams_null_cp(config$rng$formal_master_seed, length(datasets) * length(sampled))
  records <- vector("list", length(streams)); job <- 0L
  for (dataset in datasets) for (method_id in sampled) {
    job <- job + 1L
    records[[job]] <- data.frame(
      job_id = sprintf("job_%03d", job), dataset_id = dataset$dataset_id,
      raw_id = dataset$raw_id, method_id = method_id,
      rng_stream_id = sprintf("formal_stream_%03d", job)
    )
    records[[job]]$rng_state <- I(list(streams[[job]]))
  }
  do.call(rbind, records)
}

# Retain the first raw-order datasets with at least one selected changepoint.
prepare_bank_null_cp <- function(config = default_config_null_cp()) {
  streams <- rng_streams_null_cp(config$rng$raw_master_seed, config$raw_attempt_limit)
  datasets <- list(); raw_id <- 0L
  while (length(datasets) < config$dataset_count && raw_id < config$raw_attempt_limit) {
    raw_id <- raw_id + 1L
    candidate <- generate_raw_attempt_null_cp(raw_id, streams[[raw_id]], config)
    if (candidate$selected_k < 1L) next
    candidate$dataset_id <- length(datasets) + 1L
    candidate$target <- build_target_null_cp(candidate, config)
    datasets[[candidate$dataset_id]] <- candidate
  }
  if (length(datasets) != config$dataset_count) stop("The raw-attempt limit was reached before the bank was complete.")
  config$raw_attempts <- raw_id
  bank <- list(config = config, datasets = datasets, rng_assignments = assign_rng_null_cp(datasets, config), raw_attempts = raw_id)
  bank
}

summarize_null_results_cp <- function(dataset_results, alpha = 0.05, methods = NULL) {
  if (is.null(methods)) methods <- unique(as.character(dataset_results$method_id))
  do.call(
    rbind,
    lapply(methods, function(method_id) {
      rows <- dataset_results[dataset_results$method_id == method_id, , drop = FALSE]
      data.frame(
        method_id = method_id,
        datasets = nrow(rows),
        rejections = sum(rows$p_value <= alpha),
        type_i_error = mean(rows$p_value <= alpha),
        stringsAsFactors = FALSE
      )
    })
  )
}
