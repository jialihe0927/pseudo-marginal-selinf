default_config_nonnull_cp <- function() {
  list(
    application = "changepoint", delta_grid = c(0, 0.4, 0.6, 0.8, 1),
    datasets_per_delta = 200L, raw_attempt_limit = 10000L,
    n = 200L, randomization_sd = 0.1,
    signal_direction = replace(numeric(200L), 101:140, 1),
    minimum_segment = 10L, maximum_changepoints = 4L,
    penalty_coefficient = 1.25, tie_tolerance = 1e-12,
    test_level = 0.05,
    methods = c("full_selector", "data_fission"),
    sampled_methods = "full_selector",
    pmmh = list(
      warmup = 200L, stabilization = 0L, production = 5000L,
      total_transitions = 5200L, auxiliary_samples = 1L, rho = 0.8,
      scale_initial = 0.8, scale_min = 0.05, scale_max = 5,
      target_acceptance = 0.20
    ),
    rng = list(
      kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      generation_masters = c(80421041L, 80422041L, 80423041L, 80424041L),
      generation_stride = 1009L,
      null_generation_master = 20261013L,
      null_proposed_master = 20261014L,
      positive_proposed_master = 20261030L,
      complementary_bank_boundary = 100L
    )
  )
}

raw_attempt_limit_nonnull_cp <- function(config) {
  if (is.null(config$raw_attempt_limit)) 10000L else as.integer(config$raw_attempt_limit)
}

null_config_nonnull_cp <- function(config) {
  null <- default_config_null_cp()
  # Shared scientific settings; labels and dataset-count field names differ.
  shared <- c("n", "randomization_sd", "minimum_segment", "maximum_changepoints",
              "penalty_coefficient", "tie_tolerance", "test_level")
  for (field in shared) null[[field]] <- config[[field]]
  null$dataset_count <- config$datasets_per_delta
  null$raw_attempt_limit <- raw_attempt_limit_nonnull_cp(config)
  null$rng$raw_master_seed <- config$rng$null_generation_master
  null$rng$formal_master_seed <- config$rng$null_proposed_master
  # Explicit correspondence preserves the quadratic-null kernel and its RNG order.
  fields <- c(warmup = "warmup", stabilization = "stabilization",
              production = "production", total_transitions = "total_transitions",
              auxiliary_samples = "auxiliary_samples", rho = "rho",
              proposal_scale_initial = "scale_initial", proposal_scale_min = "scale_min",
              proposal_scale_max = "scale_max", target_acceptance = "target_acceptance")
  for (field in names(fields)) null$pmmh[[field]] <- config$pmmh[[fields[[field]]]]
  null$pmmh$adaptation_rate <- 0.5
  null$pmmh$adaptation_offset <- 10
  null$pmmh$adaptation_power <- 0.6
  null$pmmh$production_adaptation <- FALSE
  null
}

segments_nonnull_cp <- function(n, changepoints) {
  boundaries <- c(0L, sort(as.integer(changepoints)), as.integer(n))
  data.frame(start = head(boundaries, -1L) + 1L, end = tail(boundaries, -1L))
}

cusum_nonnull_cp <- function(series, start, split, end) {
  left <- split - start + 1L
  right <- end - split
  sqrt(left * right / (left + right)) *
    (mean(series[(split + 1L):end]) - mean(series[start:split]))
}

candidates_nonnull_cp <- function(series, segments, minimum_segment) {
  rows <- vector("list", nrow(segments))
  for (index in seq_len(nrow(segments))) {
    start <- segments$start[index]
    end <- segments$end[index]
    if (end - start + 1L < 2L * minimum_segment) next
    splits <- seq.int(start + minimum_segment - 1L, end - minimum_segment)
    signed <- vapply(splits, function(split) {
      cusum_nonnull_cp(series, start, split, end)
    }, numeric(1L))
    rows[[index]] <- data.frame(start = start, split = splits, end = end,
                                signed = signed, absolute = abs(signed))
  }
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  candidates <- do.call(rbind, rows)
  candidates[order(-candidates$absolute, candidates$start, candidates$split,
                   candidates$end, method = "radix"), , drop = FALSE]
}

rss_nonnull_cp <- function(series, changepoints) {
  segments <- segments_nonnull_cp(length(series), changepoints)
  sum(vapply(seq_len(nrow(segments)), function(index) {
    values <- series[segments$start[index]:segments$end[index]]
    sum((values - mean(values))^2)
  }, numeric(1L)))
}

# Use the full segmentation path to select K, but condition only on the final unordered set.
select_nonnull_cp <- function(series, config) {
  path <- integer(0L)
  path_signs <- integer(config$maximum_changepoints)
  rss <- numeric(config$maximum_changepoints + 1L)
  rss[1L] <- rss_nonnull_cp(series, path)
  for (step in seq_len(config$maximum_changepoints)) {
    candidates <- candidates_nonnull_cp(
      series, segments_nonnull_cp(length(series), path), config$minimum_segment
    )
    if (is.null(candidates)) stop("No legal split remains before completing the path.")
    winner <- candidates[1L, , drop = FALSE]
    path <- c(path, as.integer(winner$split))
    path_signs[step] <- as.integer(sign(winner$signed))
    rss[step + 1L] <- rss_nonnull_cp(series, path)
  }
  candidate_k <- 0:config$maximum_changepoints
  criterion <- rss + config$penalty_coefficient * candidate_k * log(length(series))
  minimizers <- which(abs(criterion - min(criterion)) <= config$tie_tolerance)
  selected_k <- as.integer(min(minimizers) - 1L)
  selected <- if (selected_k == 0L) integer(0L) else sort(path[seq_len(selected_k)])
  list(selected_k = selected_k, selected_changepoints = selected,
       observed_path = path, observed_signs = path_signs, criterion = criterion)
}

# Vary T while holding R_obs fixed, then compare the adaptively selected final set.
selection_event_nonnull_cp <- function(T, U, dataset, method_id, config) {
  target <- dataset$target
  Y <- drop(target$R_obs + target$V %*% (target$Sigma_T_inv %*% as.numeric(T)))
  selected <- select_nonnull_cp(Y + config$randomization_sd * U, config)
  identical(selected$selected_changepoints, dataset$selected_changepoints)
}

# Conditional Gaussian target and fixed directional statistic.
build_target_nonnull_cp <- function(dataset, config) {
  segments <- segments_nonnull_cp(config$n, dataset$selected_changepoints)
  V <- matrix(0, config$n, dataset$selected_k)
  for (index in seq_len(dataset$selected_k)) {
    left <- seq.int(segments$start[index], segments$end[index])
    right <- seq.int(segments$start[index + 1L], segments$end[index + 1L])
    V[right, index] <- 1 / length(right)
    V[left, index] <- V[left, index] - 1 / length(left)
  }
  Sigma_T <- crossprod(V)
  inverse <- solve(Sigma_T)
  T_obs <- drop(crossprod(V, dataset$Y))
  R_obs <- drop(dataset$Y - V %*% (inverse %*% T_obs))
  list(V = V, Sigma_T = Sigma_T, Sigma_T_inv = inverse, T_obs = T_obs, R_obs = R_obs)
}

# Project the fixed signal template onto the selected contrasts; delta affects only the DGP.
build_direction_nonnull_cp <- function(dataset, config) {
  target <- dataset$target
  a <- drop(crossprod(target$V, config$signal_direction))
  coefficient <- as.numeric(solve(target$Sigma_T, a))
  information <- drop(crossprod(a, coefficient))
  # Exact integer overlap identities distinguish a structural zero direction.
  segments <- segments_nonnull_cp(config$n, dataset$selected_changepoints)
  lengths <- segments$end - segments$start + 1L
  overlap <- vapply(seq_len(nrow(segments)), function(index) {
    sum(config$signal_direction[seq.int(segments$start[index], segments$end[index])])
  }, numeric(1L))
  zero <- all(overlap[-1L] * head(lengths, -1L) == head(overlap, -1L) * lengths[-1L])
  if (!zero && (!is.finite(information) || information <= 0)) stop("Invalid directional information.")
  list(
    a = a, coefficient = coefficient, information = information, structural_zero = zero,
    observed = if (zero) NA_real_ else drop(crossprod(a, solve(target$Sigma_T, target$T_obs)) / sqrt(information))
  )
}

statistic_nonnull_cp <- function(T, direction) {
  if (direction$structural_zero) stop("The directional statistic is undefined for a structural zero.")
  if (is.null(dim(T))) {
    return(drop(crossprod(direction$coefficient, T) / sqrt(direction$information)))
  }
  drop(as.matrix(T) %*% direction$coefficient) / sqrt(direction$information)
}

fission_statistic_nonnull_cp <- function(dataset, config) {
  if (dataset$direction$structural_zero) return(NA_real_)
  # Selection uses Y+W; inference uses Y-tau^(-2)W, independent under Gaussian noise.
  Y_fission <- dataset$Y - 100 * dataset$W
  T_fission <- drop(crossprod(dataset$target$V, Y_fission))
  direction <- dataset$direction
  drop(crossprod(direction$a, solve(dataset$target$Sigma_T, T_fission)) /
         sqrt(direction$information))
}

fission_pvalue_nonnull_cp <- function(
    dataset, config, score = fission_statistic_nonnull_cp(dataset, config)) {
  if (dataset$direction$structural_zero) return(1)
  stats::pnorm(score / sqrt(101), lower.tail = FALSE)
}

with_rng_nonnull_cp <- function(fun, seed = NULL, state = NULL) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  if (is.null(state)) set.seed(as.integer(seed))
  else assign(".Random.seed", as.integer(state), envir = .GlobalEnv)
  fun()
}

convert_null_nonnull_cp <- function(dataset, config) {
  converted <- list(
    delta = 0, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
    Y = dataset$response, W = dataset$randomization, selected_k = dataset$selected_k,
    selected_changepoints = dataset$selected_changepoints,
    observed_path = dataset$ordered_path, observed_signs = integer(),
    criterion = dataset$criterion, generation_rng_state = dataset$generation_rng_state,
    generation_rng = dataset$generation_rng,
    source_batch = 0L
  )
  converted$target <- build_target_nonnull_cp(converted, config)
  converted$direction <- build_direction_nonnull_cp(converted, config)
  converted
}

generate_attempt_nonnull_cp <- function(delta, raw_id, config) {
  group <- match(delta, config$delta_grid[config$delta_grid > 0])
  seed <- as.integer((config$rng$generation_masters[group] +
                        raw_id * config$rng$generation_stride) %% 2147483627L)
  values <- with_rng_nonnull_cp(function() {
    epsilon <- rnorm(config$n)
    W <- rnorm(config$n, sd = config$randomization_sd)
    list(epsilon = epsilon, W = W)
  }, seed = seed)
  Y <- delta * config$signal_direction + values$epsilon
  observed <- select_nonnull_cp(Y + values$W, config)
  if (observed$selected_k < 1L) return(NULL)
  dataset <- c(list(delta = delta, raw_id = as.integer(raw_id), Y = Y,
                   W = values$W, generation_seed = seed), observed)
  dataset$target <- build_target_nonnull_cp(dataset, config)
  dataset$direction <- build_direction_nonnull_cp(dataset, config)
  dataset
}

assign_rng_nonnull_cp <- function(datasets, config) {
  rows <- list()
  null_streams <- rng_streams_null_cp(config$rng$null_proposed_master, 2L * config$datasets_per_delta)
  # Allocate positive-signal streams by scientific identity, independent of
  # worker scheduling or the order of datasets supplied to this function.
  positive_jobs <- expand.grid(
    dataset_id = seq_len(config$datasets_per_delta),
    delta = sort(config$delta_grid[config$delta_grid > 0])
  )
  positive_keys <- paste(positive_jobs$delta, positive_jobs$dataset_id, sep = ":")
  positive_streams <- with_rng_nonnull_cp(function() {
    streams <- vector("list", nrow(positive_jobs))
    streams[[1L]] <- as.integer(get(".Random.seed", envir = .GlobalEnv))
    if (length(streams) > 1L) {
      for (index in seq.int(2L, length(streams))) {
        streams[[index]] <- parallel::nextRNGStream(streams[[index - 1L]])
      }
    }
    streams
  }, seed = config$rng$positive_proposed_master)
  for (dataset in datasets) for (method_id in config$sampled_methods) {
    rng_batch <- 1L + as.integer(dataset$dataset_id > config$rng$complementary_bank_boundary)
    seed <- if (dataset$delta == 0) config$rng$null_proposed_master else config$rng$positive_proposed_master
    row <- data.frame(
      delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
      method_id = method_id,
      source_batch = if (dataset$delta == 0) 0L else rng_batch, seed = seed,
      stringsAsFactors = FALSE
    )
    state <- if (dataset$delta == 0) {
      null_streams[[2L * dataset$dataset_id - 1L]]
    } else {
      index <- match(paste(dataset$delta, dataset$dataset_id, sep = ":"), positive_keys)
      if (is.na(index)) stop("Unknown positive-signal task identity.")
      positive_streams[[index]]
    }
    row$rng_state <- I(list(state))
    rows[[length(rows) + 1L]] <- row
  }
  do.call(rbind, rows)
}

# Rebuild delta-zero data with the calibration generator; keep positive-signal data in raw order.
prepare_bank_nonnull_cp <- function(config = default_config_nonnull_cp()) {
  null_config <- null_config_nonnull_cp(config)
  attempt_limit <- raw_attempt_limit_nonnull_cp(config)
  datasets <- list()
  raw_attempts <- integer(length(config$delta_grid))
  for (delta_index in seq_along(config$delta_grid)) {
    if (config$delta_grid[delta_index] == 0) {
      null_bank <- prepare_bank_null_cp(null_config)
      datasets <- c(datasets, lapply(null_bank$datasets, convert_null_nonnull_cp, config = config))
      raw_attempts[delta_index] <- null_bank$raw_attempts
      next
    }
    retained <- 0L
    raw_id <- 0L
    while (retained < config$datasets_per_delta && raw_id < attempt_limit) {
      raw_id <- raw_id + 1L
      dataset <- generate_attempt_nonnull_cp(config$delta_grid[delta_index], raw_id, config)
      if (is.null(dataset)) next
      retained <- retained + 1L
      dataset$dataset_id <- retained
      dataset$source_batch <- 1L + as.integer(retained > config$rng$complementary_bank_boundary)
      datasets[[length(datasets) + 1L]] <- dataset
    }
    if (retained < config$datasets_per_delta) {
      stop("Raw-attempt limit reached for delta=", config$delta_grid[delta_index],
           ": retained ", retained, "/", config$datasets_per_delta,
           " eligible datasets after ", raw_id, " attempts.")
    }
    raw_attempts[delta_index] <- raw_id
  }
  bank <- list(config = config, datasets = datasets,
               rng_assignments = assign_rng_nonnull_cp(datasets, config), raw_attempts = raw_attempts)
  bank
}

result_row_nonnull_cp <- function(dataset, method_id, p_value, config,
                                         tail_count = NA_integer_,
                                         production_draws = NA_integer_,
                                         result_type = "analytic",
                                         selector_evaluations = 0L,
                                         statistic_observed = NULL) {
  if (is.null(statistic_observed)) {
    statistic_observed <- if (method_id == "data_fission") {
      fission_statistic_nonnull_cp(dataset, config)
    } else dataset$direction$observed
  }
  data.frame(
    delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
    method_id = method_id, target_id = "selected_changepoint_direction",
    target_dimension = dataset$selected_k,
    statistic_observed = statistic_observed,
    p_value = p_value, reject = as.integer(p_value <= config$test_level),
    tail_count = tail_count, production_draws = production_draws,
    result_type = result_type, status = "complete",
    selector_evaluations = selector_evaluations, stringsAsFactors = FALSE
  )
}

summarize_nonnull_cp <- function(data, config) {
  rows <- list()
  for (delta in config$delta_grid) for (method_id in config$methods) {
    selected <- data[data$delta == delta & data$method_id == method_id, , drop = FALSE]
    n <- nrow(selected)
    rejections <- sum(selected$reject)
    power <- rejections / n
    # Power error bars are +/- one standard error across datasets.
    se <- sqrt(power * (1 - power) / n)
    rows[[length(rows) + 1L]] <- data.frame(
      delta = delta, method_id = method_id, n = n, rejections = rejections,
      power = power, se = se, lower = power - se, upper = power + se,
      stringsAsFactors = FALSE
    )
  }
  list(power = do.call(rbind, rows))
}
