# Configuration.
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

validate_config_null_cp <- function(config) {
  required <- c(
    "n", "randomization_sd", "minimum_segment", "maximum_changepoints",
    "penalty_coefficient", "dataset_count", "methods", "pmmh"
  )
  if (!all(required %in% names(config))) {
    stop("The changepoint configuration is incomplete.")
  }

  pmmh <- config$pmmh
  if (pmmh$warmup + pmmh$stabilization + pmmh$production != pmmh$total_transitions) {
    stop("Warm-up, stabilization, and production must sum to total_transitions.")
  }
  if (pmmh$stabilization != 0L) {
    stop("This implementation supports stabilization = 0 only.")
  }
  if (pmmh$production <= 0L || pmmh$auxiliary_samples != 1L) {
    stop("Production must be positive and auxiliary_samples must equal one.")
  }
  if (pmmh$proposal_scale_min <= 0 ||
      pmmh$proposal_scale_min > pmmh$proposal_scale_initial ||
      pmmh$proposal_scale_initial > pmmh$proposal_scale_max) {
    stop("Proposal scales must satisfy 0 < min <= initial <= max.")
  }
  if (!is.finite(pmmh$rho) || pmmh$rho < 0 || pmmh$rho >= 1) {
    stop("The auxiliary correlation must lie in [0, 1).")
  }
  invisible(config)
}

sampled_method_ids_null_cp <- function(config) {
  if (!is.null(config$sampled_method_ids)) return(as.character(config$sampled_method_ids))
  setdiff(as.character(config$methods$method_id), "naive_gaussian")
}

# Adaptive and fixed-K events retain the original split ordering and ties.
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
  rows <- vector("list", config$maximum_changepoints)
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
    rows[[step]] <- data.frame(step = step, winner, row.names = NULL)
    rss[step + 1L] <- piecewise_rss_null_cp(series, path)
  }

  list(path = path, path_table = do.call(rbind, rows), rss = rss)
}

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
  if (selected_k < 1L || selected_k > config$maximum_changepoints) {
    stop("Fixed K must lie between one and maximum_changepoints.")
  }

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

# Conditional Gaussian reconstruction and quadratic target.
reconstruct_response_null_cp <- function(target_value, target) {
  as.numeric(
    target$R_obs +
    target$V %*% (target$Sigma_T_inv %*% as.numeric(target_value))
  )
}

target_statistic_null_cp <- function(target_value, target) {
  drop(crossprod(target_value, target$Sigma_T_inv %*% target_value))
}

build_target_null_cp <- function(dataset, config) {
  selected <- sort(as.integer(dataset$selected_changepoints))
  selected_k <- as.integer(dataset$selected_k)
  if (selected_k < 1L || selected_k > config$maximum_changepoints || length(selected) != selected_k) stop("Target requires an eligible changepoint selection.")
  segments <- segments_null_cp(config$n, selected)
  if (any(segments$end - segments$start + 1L < config$minimum_segment)) stop("Observed segments violate the minimum length.")
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

# Global-null bank and preassigned streams.
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

assign_rng_null_cp <- function(datasets, config) {
  sampled <- sampled_method_ids_null_cp(config)
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

prepare_bank_null_cp <- function(config = default_config_null_cp()) {
  validate_config_null_cp(config)
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
  validate_bank_null_cp(bank)
  bank
}

# Paired summaries and frozen-input checks.
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

validate_bank_null_cp <- function(bank) {
  required <- c("config", "datasets", "rng_assignments")
  if (!all(required %in% names(bank))) stop("The changepoint data bank is incomplete.")
  config <- bank$config; validate_config_null_cp(config)
  datasets <- bank$datasets
  if (length(datasets) != config$dataset_count) stop("The data bank has the wrong number of datasets.")
  ids <- vapply(datasets, `[[`, integer(1L), "dataset_id")
  raw_ids <- vapply(datasets, `[[`, integer(1L), "raw_id")
  selected_k <- vapply(datasets, `[[`, integer(1L), "selected_k")
  if (!identical(ids, seq_len(config$dataset_count)) || anyDuplicated(raw_ids) || !all(diff(raw_ids) > 0L)) stop("Dataset or raw-attempt identities are invalid.")
  if (any(selected_k < 1L | selected_k > config$maximum_changepoints)) stop("A retained dataset violates the K selection rule.")
  assignments <- bank$rng_assignments
  expected <- expand.grid(dataset_id = ids, method_id = sampled_method_ids_null_cp(config), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- paste(assignments$dataset_id, assignments$method_id, sep = ":")
  expected_key <- paste(expected$dataset_id, expected$method_id, sep = ":")
  if (nrow(assignments) != nrow(expected) || anyDuplicated(key) || !setequal(key, expected_key) ||
      !all(assignments$raw_id == raw_ids[assignments$dataset_id]) || length(assignments$rng_state) != nrow(assignments)) {
    stop("The inference random-stream assignments are incomplete or inconsistent.")
  }
  invisible(bank)
}

validate_result_null_cp <- function(result) {
  required <- c("config", "datasets", "rng_assignments", "dataset_results", "summary")
  if (!all(required %in% names(result))) {
    stop("The changepoint result object is incomplete.")
  }
  validate_config_null_cp(result$config)

  validate_bank_null_cp(result[c("config", "datasets", "rng_assignments")])
  data <- result$dataset_results
  method_ids <- as.character(result$config$methods$method_id)
  expected <- expand.grid(dataset_id = seq_len(result$config$dataset_count), method_id = method_ids,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- paste(data$dataset_id, data$method_id, sep = ":")
  expected_key <- paste(expected$dataset_id, expected$method_id, sep = ":")
  if (nrow(data) != nrow(expected) || anyDuplicated(key) || !setequal(key, expected_key)) {
    stop("Results must contain exactly one paired record for every dataset-method combination.")
  }
  if (any(!is.finite(data$p_value)) || any(data$p_value < 0 | data$p_value > 1)) {
    stop("P-values must be finite and lie in [0, 1].")
  }

  dimensions <- vapply(result$datasets, `[[`, integer(1L), "selected_k")
  if (any(data$target_dimension != dimensions[data$dataset_id])) stop("A result target dimension does not match its dataset.")
  if (any(data$status != "complete")) stop("Incomplete results cannot enter the formal summary.")
  sampled <- data$method_id %in% sampled_method_ids_null_cp(result$config)
  if (any(data$result_type[sampled] != "simulation") || any(data$result_type[!sampled] != "analytic")) {
    stop("Result types do not match the configured methods.")
  }
  if (any(!is.na(data$tail_count[!sampled])) || any(!is.na(data$production_draws[!sampled]))) {
    stop("Analytic methods must not carry sampling tail counts.")
  }
  if (any(!is.finite(data$tail_count[sampled])) ||
      any(data$tail_count[sampled] != as.integer(data$tail_count[sampled])) ||
      any(data$tail_count[sampled] < 0L) ||
      any(data$tail_count[sampled] > data$production_draws[sampled])) {
    stop("Sampling tail counts must be valid integers within their denominators.")
  }
  if (any(data$production_draws[sampled] != result$config$pmmh$production)) {
    stop("Sampling denominators do not match the configured production length.")
  }
  reconstructed_p <- data$tail_count[sampled] / data$production_draws[sampled]
  if (any(abs(reconstructed_p - data$p_value[sampled]) > 1e-15)) {
    stop("A sampling p-value does not match its tail count and denominator.")
  }
  if (any(as.integer(data$reject) != as.integer(data$p_value <= result$config$test_level))) stop("Rejection indicators do not match the configured test level.")
  summary <- summarize_null_results_cp(data, result$config$test_level, method_ids)
  if (!isTRUE(all.equal(summary, result$summary, check.attributes = FALSE))) stop("The summary does not match the dataset-level results.")
  invisible(result)
}
