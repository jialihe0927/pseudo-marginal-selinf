default_config_nonnull_cp <- function() {
  list(
    application = "changepoint", delta_grid = c(0.4, 0.6, 0.8, 1),
    datasets_per_delta = 200L, n = 200L, randomization_sd = 0.1,
    signal_direction = replace(numeric(200L), 101:140, 1),
    minimum_segment = 10L, maximum_changepoints = 4L,
    penalty_coefficient = 1.25, tie_tolerance = 1e-12,
    test_level = 0.05,
    methods = c("full_selector", "full_path", "data_fission"),
    sampled_methods = c("full_selector", "full_path"),
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
      inference_bases = c(126000000L, 226000000L),
      complementary_bank_boundary = 100L, inference_stride = 7919L
    )
  )
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

build_target_nonnull_cp <- function(dataset, config) {
  segments <- segments_nonnull_cp(config$n, dataset$selected_changepoints)
  if (any(segments$end - segments$start + 1L < config$minimum_segment)) {
    stop("Observed segments violate the minimum length.")
  }
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

selection_event_nonnull_cp <- function(T, U, dataset, method_id, config) {
  if (!method_id %in% config$sampled_methods) stop("Unknown changepoint sampling method: ", method_id)
  target <- dataset$target
  Y <- drop(target$R_obs + target$V %*% (target$Sigma_T_inv %*% as.numeric(T)))
  selected <- select_nonnull_cp(Y + config$randomization_sd * U, config)
  if (method_id == "full_selector") {
    return(identical(selected$selected_changepoints, dataset$selected_changepoints))
  }
  # All four ordered split locations and signs are fixed, including those after K_obs.
  selected$selected_k == dataset$selected_k &&
    identical(selected$observed_path, dataset$observed_path) &&
    identical(selected$observed_signs, dataset$observed_signs)
}

generate_attempt_nonnull_cp <- function(delta, raw_id, config) {
  group <- match(delta, config$delta_grid)
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
  for (dataset in datasets) for (method_id in config$sampled_methods) {
    group <- match(dataset$delta, config$delta_grid)
    method <- match(method_id, config$sampled_methods)
    source_batch <- 1L + as.integer(dataset$dataset_id > config$rng$complementary_bank_boundary)
    seed <- as.integer(config$rng$inference_bases[source_batch] + group * 10000000L +
                         method * 1000000L + dataset$dataset_id * config$rng$inference_stride)
    row <- data.frame(
      delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
      method_id = method_id, source_batch = source_batch, seed = seed,
      stringsAsFactors = FALSE
    )
    row$rng_state <- I(list(with_rng_nonnull_cp(function() {
      as.integer(get(".Random.seed", envir = .GlobalEnv))
    }, seed = seed)))
    rows[[length(rows) + 1L]] <- row
  }
  do.call(rbind, rows)
}

prepare_bank_nonnull_cp <- function(config = default_config_nonnull_cp()) {
  validate_config_nonnull_cp(config)
  datasets <- list()
  raw_attempts <- integer(length(config$delta_grid))
  for (delta_index in seq_along(config$delta_grid)) {
    retained <- 0L
    raw_id <- 0L
    while (retained < config$datasets_per_delta) {
      raw_id <- raw_id + 1L
      dataset <- generate_attempt_nonnull_cp(config$delta_grid[delta_index], raw_id, config)
      if (is.null(dataset)) next
      retained <- retained + 1L
      dataset$dataset_id <- retained
      dataset$source_batch <- 1L + as.integer(retained > config$rng$complementary_bank_boundary)
      datasets[[length(datasets) + 1L]] <- dataset
    }
    raw_attempts[delta_index] <- raw_id
  }
  bank <- list(config = config, datasets = datasets,
               rng_assignments = assign_rng_nonnull_cp(datasets, config), raw_attempts = raw_attempts)
  validate_bank_nonnull_cp(bank)
  bank
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

fission_pvalue_nonnull_cp <- function(dataset, config) {
  if (dataset$direction$structural_zero) return(1)
  score <- fission_statistic_nonnull_cp(dataset, config)
  stats::pnorm(score / sqrt(101), lower.tail = FALSE)
}
validate_config_nonnull_cp <- function(config) {
  reference <- default_config_nonnull_cp()
  if (!identical(config$application, "changepoint") ||
      !identical(config$methods, reference$methods) ||
      !identical(config$sampled_methods, reference$sampled_methods)) {
    stop("The application or method definitions are incompatible.")
  }
  n <- config$datasets_per_delta
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop("datasets_per_delta must be a positive integer.")
  }
  if (!identical(config$delta_grid, reference$delta_grid) ||
      length(config$signal_direction) != config$n) {
    stop("The signal grid or direction is incompatible.")
  }
  pmmh <- config$pmmh
  lengths <- unlist(pmmh[c("warmup", "stabilization", "production", "total_transitions")])
  if (any(!is.finite(lengths)) || any(lengths != as.integer(lengths)) ||
      pmmh$warmup <= 0L || pmmh$production <= 0L || pmmh$stabilization != 0L ||
      pmmh$warmup + pmmh$production != pmmh$total_transitions ||
      pmmh$auxiliary_samples != 1L) {
    stop("Invalid PMMH lengths; stabilization must be zero and auxiliary_samples one.")
  }
  if (!is.finite(pmmh$rho) || pmmh$rho < 0 || pmmh$rho >= 1 ||
      pmmh$scale_min <= 0 || pmmh$scale_min > pmmh$scale_initial ||
      pmmh$scale_initial > pmmh$scale_max || pmmh$target_acceptance <= 0 ||
      pmmh$target_acceptance >= 1) {
    stop("Invalid auxiliary correlation, proposal scales, or adaptation target.")
  }
  if (config$randomization_sd != 0.1 || config$minimum_segment != 10L ||
      config$maximum_changepoints != 4L || config$penalty_coefficient != 1.25) {
    stop("This implementation supports the recorded bounded selector and fission scaling only.")
  }
  invisible(config)
}

validate_bank_nonnull_cp <- function(bank) {
  if (!all(c("config", "datasets", "rng_assignments") %in% names(bank))) {
    stop("The data bank is incomplete.")
  }
  config <- bank$config
  validate_config_nonnull_cp(config)
  datasets <- bank$datasets
  if (length(datasets) != length(config$delta_grid) * config$datasets_per_delta) {
    stop("The data bank has the wrong number of datasets.")
  }
  identity <- do.call(rbind, lapply(datasets, function(dataset) {
    required <- c("delta", "dataset_id", "raw_id", "Y", "W", "selected_k", "target", "direction")
    if (!all(required %in% names(dataset))) stop("A dataset is incomplete.")
    if (length(dataset$Y) != config$n || length(dataset$W) != config$n ||
        any(!is.finite(c(dataset$Y, dataset$W, dataset$target$T_obs)))) {
      stop("Dataset inputs or observed target values are invalid.")
    }
    data.frame(delta = dataset$delta, dataset_id = dataset$dataset_id,
               raw_id = dataset$raw_id, dimension = dataset$selected_k)
  }))
  expected <- expand.grid(delta = config$delta_grid, dataset_id = seq_len(config$datasets_per_delta))
  key <- paste(identity$delta, identity$dataset_id, sep = ":")
  if (anyNA(identity) || anyDuplicated(key) ||
      !setequal(key, paste(expected$delta, expected$dataset_id, sep = ":"))) {
    stop("Dataset identities do not cover the signal grid exactly.")
  }
  for (delta in config$delta_grid) {
    rows <- identity[identity$delta == delta, , drop = FALSE]
    if (!identical(as.integer(rows$dataset_id), seq_len(config$datasets_per_delta)) ||
        any(rows$raw_id < 1L | rows$raw_id != as.integer(rows$raw_id)) ||
        anyDuplicated(rows$raw_id) || any(diff(rows$raw_id) <= 0L)) {
      stop("Raw attempts must be unique and retained in increasing order within each delta.")
    }
    if (any(rows$dimension < 1L |
            rows$dimension > config$maximum_changepoints)) {
      stop("A retained dataset violates the eligibility rule.")
    }
  }
  for (dataset in datasets) {
    if (length(dataset$selected_changepoints) != dataset$selected_k ||
        length(dataset$observed_path) != config$maximum_changepoints ||
        length(dataset$observed_signs) != config$maximum_changepoints ||
        !identical(sort(dataset$observed_path[seq_len(dataset$selected_k)]),
                   dataset$selected_changepoints)) {
      stop("The observed path does not reconstruct the selected changepoints.")
    }
    covariance <- dataset$target$Sigma_T
    required_target <- c("V", "T_obs", "R_obs", "Sigma_T", "Sigma_T_inv")
    if (!all(required_target %in% names(dataset$target)) ||
        !identical(dim(covariance), rep(as.integer(dataset$selected_k), 2L)) ||
        any(!is.finite(covariance))) {
      stop("A conditional target has incompatible fields or covariance.")
    }
  }
  assignments <- bank$rng_assignments
  required <- c("delta", "dataset_id", "raw_id", "method_id", "rng_state")
  if (!all(required %in% names(assignments))) stop("Random-stream assignments are incomplete.")
  expected_jobs <- expand.grid(delta = config$delta_grid,
                               dataset_id = seq_len(config$datasets_per_delta),
                               method_id = config$sampled_methods, stringsAsFactors = FALSE)
  job_key <- paste(assignments$delta, assignments$dataset_id, assignments$method_id, sep = ":")
  expected_key <- paste(expected_jobs$delta, expected_jobs$dataset_id, expected_jobs$method_id, sep = ":")
  match_dataset <- match(paste(assignments$delta, assignments$dataset_id, sep = ":"), key)
  if (nrow(assignments) != nrow(expected_jobs) || anyNA(assignments[, required[-length(required)]]) ||
      anyDuplicated(job_key) || !setequal(job_key, expected_key) ||
      anyNA(match_dataset) || any(assignments$raw_id != identity$raw_id[match_dataset]) ||
      !all(vapply(assignments$rng_state, function(state) {
        is.integer(state) && length(state) == 7L && !anyNA(state)
      }, logical(1L)))) {
    stop("Random-stream assignments must cover every paired sampling task exactly once.")
  }
  invisible(bank)
}

summarize_nonnull_cp <- function(data, config) {
  rows <- list()
  for (delta in config$delta_grid) for (method_id in config$methods) {
    selected <- data[data$delta == delta & data$method_id == method_id, , drop = FALSE]
    n <- nrow(selected)
    rejections <- sum(selected$reject)
    power <- rejections / n
    se <- sqrt(power * (1 - power) / n)
    rows[[length(rows) + 1L]] <- data.frame(
      delta = delta, method_id = method_id, n = n, rejections = rejections,
      power = power, se = se, lower = power - se, upper = power + se,
      stringsAsFactors = FALSE
    )
  }
  list(power = do.call(rbind, rows))
}

validate_result_nonnull_cp <- function(result) {
  validate_bank_nonnull_cp(result)
  if (!all(c("dataset_results", "summary") %in% names(result))) stop("Inference results are missing.")
  config <- result$config
  data <- result$dataset_results
  required <- c("delta", "dataset_id", "raw_id", "method_id", "target_id",
                "target_dimension", "statistic_observed", "p_value", "reject",
                "tail_count", "production_draws", "result_type", "status")
  if (!all(required %in% names(data))) stop("Dataset-level result fields are incomplete.")
  expected <- expand.grid(delta = config$delta_grid, dataset_id = seq_len(config$datasets_per_delta),
                          method_id = config$methods, stringsAsFactors = FALSE)
  key <- paste(data$delta, data$dataset_id, data$method_id, sep = ":")
  expected_key <- paste(expected$delta, expected$dataset_id, expected$method_id, sep = ":")
  if (nrow(data) != nrow(expected) || anyDuplicated(key) ||
      !setequal(key, expected_key) || anyNA(data$status) || any(data$status != "complete")) {
    stop("Results must contain every paired dataset-method combination with complete status.")
  }
  identity <- vapply(result$datasets, function(dataset) paste(dataset$delta, dataset$dataset_id, sep = ":"), character(1L))
  positions <- match(paste(data$delta, data$dataset_id, sep = ":"), identity)
  dimension <- vapply(result$datasets, `[[`, integer(1L), "selected_k")
  raw_id <- vapply(result$datasets, `[[`, integer(1L), "raw_id")
  if (anyNA(positions) || anyNA(data$target_dimension) || anyNA(data$raw_id) ||
      any(data$target_dimension != dimension[positions]) || any(data$raw_id != raw_id[positions])) {
    stop("Result identities or target dimensions do not match the bank.")
  }
  if (any(!is.finite(data$p_value)) || any(data$p_value < 0 | data$p_value > 1) ||
      anyNA(data$reject) || any(data$reject != as.integer(data$p_value <= config$test_level))) {
    stop("Invalid p-values or rejection indicators.")
  }
  sampled <- data$method_id %in% config$sampled_methods
  zero <- vapply(result$datasets, function(dataset) dataset$direction$structural_zero, logical(1L))[positions]
  if (any(data$p_value[zero] != 1) || any(data$result_type[zero] != "structural_zero")) {
    stop("Structural-zero directions must retain their recorded p=1 convention.")
  }
  sampled <- sampled & !zero
  if (anyNA(data$result_type) || any(data$result_type[!sampled & !zero] != "analytic")) {
    stop("Nondegenerate Data Fission must be identified as an analytic result.")
  }
  if (any(data$result_type[sampled] != "simulation") ||
      any(!is.finite(data$statistic_observed[sampled])) ||
      any(!is.finite(data$tail_count[sampled])) ||
      any(data$tail_count[sampled] != as.integer(data$tail_count[sampled])) ||
      any(data$tail_count[sampled] < 0 | data$tail_count[sampled] > config$pmmh$production) ||
      anyNA(data$production_draws[sampled]) ||
      any(data$production_draws[sampled] != config$pmmh$production) ||
      any(abs(data$p_value[sampled] - data$tail_count[sampled] / data$production_draws[sampled]) > 1e-15)) {
    stop("Sampling p-values must match their valid inclusive-tail counts and production denominators.")
  }
  if (any(!is.na(data$tail_count[!sampled])) || any(!is.na(data$production_draws[!sampled]))) {
    stop("Analytic or degenerate results must not contain fabricated sampling counts.")
  }
  analytic <- data$method_id == "data_fission" & !zero
  expected_fission <- stats::pnorm(
    data$statistic_observed[analytic] / sqrt(101), lower.tail = FALSE
  )
  if (any(!is.finite(data$statistic_observed[analytic])) ||
      any(abs(expected_fission - data$p_value[analytic]) > 1e-15)) {
    stop("Data Fission p-values do not match their observed statistic and null variance.")
  }
  reconstructed <- summarize_nonnull_cp(data, config)$power
  if (!isTRUE(all.equal(reconstructed, result$summary$power, tolerance = 1e-15, check.attributes = FALSE))) {
    stop("The power summary does not match the paired dataset-level decisions.")
  }
  invisible(result)
}

result_row_nonnull_cp <- function(dataset, method_id, p_value, config,
                                         tail_count = NA_integer_,
                                         production_draws = NA_integer_,
                                         result_type = "analytic",
                                         selector_evaluations = 0L) {
  data.frame(
    delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
    method_id = method_id, target_id = "selected_changepoint_direction",
    target_dimension = dataset$selected_k,
    statistic_observed = if (method_id == "data_fission") {
      fission_statistic_nonnull_cp(dataset, config)
    } else {
      dataset$direction$observed
    },
    p_value = p_value, reject = as.integer(p_value <= config$test_level),
    tail_count = tail_count, production_draws = production_draws,
    result_type = result_type, status = "complete",
    selector_evaluations = selector_evaluations, stringsAsFactors = FALSE
  )
}
