# Configuration and delta-zero compatibility

default_config_nonnull_lasso <- function() {
  list(
    application = "cv_lasso", delta_grid = c(0, 0.2, 0.3, 0.4, 0.5),
    datasets_per_delta = 200L, n = 200L, p = 100L,
    design_correlation = 0.5, response_sd = 1, randomization_sd = 0.15,
    signal_direction = c(1, -0.8, 0.6, -0.5, 0.4, rep(0, 95L)),
    folds = 10L, dimension_range = c(1L, 10L), raw_attempt_limit = 10000L,
    lambda_grid = c(
      0.230497466833294, 0.218176011289412, 0.20651321056202701,
      0.19547385564796799, 0.18502461967393699, 0.17513395728550901,
      0.16577200941439399, 0.15691051313644699, 0.14852271634830799,
      0.14058329700507699, 0.13306828667521001, 0.125954998181871,
      0.11922195711227, 0.112848836988239, 0.106816397902325,
      0.101106428434145, 0.095701690671655695, 0.090585868171367906,
      0.085743516700387698, 0.081160017611592999, 0.076821533711181597,
      0.072714967485359702, 0.068827921560062594, 0.065148661274337105,
      0.061666079254398098, 0.058369661881412099, 0.0552494575517768,
      0.0522960466340762, 0.049500513032014699, 0.046854417267481199,
      0.044349771002482198, 0.041979012923029503
    ),
    coefficient_tolerance = 1e-8, glmnet_threshold = 1e-7,
    rank_tolerance = 1e-10, test_level = 0.05,
    methods = c("model_only", "full_conditioning", "data_fission"),
    sampled_methods = c("full_conditioning", "model_only"),
    pmmh = list(
      warmup = 200L, stabilization = 0L, production = 5000L,
      total_transitions = 5200L, auxiliary_samples = 1L, rho = 0.8,
      scale_initial = 0.8, scale_min = 0.05, scale_max = 5,
      target_acceptance = 0.20
    ),
    rng = list(
      kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      generation_seeds = c(20461016L, 20561016L, 20661016L, 20761016L),
      inference_master = 20261022L,
      null_generation_master = 2026999901L,
      null_proposed_master = 2026999902L,
      null_extra_master = 20261028L,
      # The preceding signal group occupied streams 1:400 in the formal run.
      inference_stream_offset = 400L
    )
  )
}

validate_config_nonnull_lasso <- function(config) {
  if (!is.numeric(config$response_sd) || length(config$response_sd) != 1L ||
      !is.finite(config$response_sd) || config$response_sd != 1) {
    stop("CV-Lasso inference currently requires response_sd = 1 (unit response noise).")
  }
  reference <- default_config_nonnull_lasso()
  if (!identical(config$application, "cv_lasso") ||
      !identical(config$methods, reference$methods) ||
      !identical(config$sampled_methods, reference$sampled_methods)) {
    stop("The application or method definitions are incompatible.")
  }
  n <- config$datasets_per_delta
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop("datasets_per_delta must be a positive integer.")
  }
  if (!identical(config$delta_grid, reference$delta_grid) ||
      length(config$signal_direction) != config$p) {
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
  if (config$n %% config$folds != 0L || length(config$lambda_grid) != 32L ||
      any(!is.finite(config$lambda_grid)) || any(diff(config$lambda_grid) >= 0)) {
    stop("Invalid fixed folds or lambda grid.")
  }
  null_config_nonnull_lasso(config)
  invisible(config)
}

null_config_nonnull_lasso <- function(config) {
  # Delta zero reuses the null bank and kernel only for these shared settings.
  null <- default_config_null_lasso()
  expected_rng <- unname(unlist(null$rng[c("kind", "normal_kind", "sample_kind")]))
  if (!identical(config$rng$kind, expected_rng)) {
    stop("Delta-zero reuse requires the null RNG kinds: L'Ecuyer-CMRG, Inversion, Rejection.")
  }
  shared <- c(
    "n", "p", "design_correlation", "response_sd", "randomization_sd",
    "folds", "lambda_grid", "coefficient_tolerance", "glmnet_threshold",
    "dimension_range"
  )
  for (field in shared) {
    if (!is.numeric(config[[field]]) ||
        !identical(as.numeric(config[[field]]), as.numeric(null[[field]]))) {
      stop("Delta-zero reuse requires the null configuration's ", field, ".")
    }
  }

  null$dataset_count <- config$datasets_per_delta
  null$test_level <- config$test_level
  null$rng$raw_master_seed <- config$rng$null_generation_master
  null$rng$formal_master_seed <- config$rng$null_proposed_master
  # Different field names refer to the same scalar proposal and adaptation settings.
  pmmh <- config$pmmh
  null$pmmh$warmup <- pmmh$warmup
  null$pmmh$stabilization <- pmmh$stabilization
  null$pmmh$production <- pmmh$production
  null$pmmh$total_transitions <- pmmh$total_transitions
  null$pmmh$auxiliary_samples <- pmmh$auxiliary_samples
  null$pmmh$rho <- pmmh$rho
  null$pmmh$target_acceptance <- pmmh$target_acceptance
  null$pmmh$proposal_scale_initial <- pmmh$scale_initial
  null$pmmh$proposal_scale_min <- pmmh$scale_min
  null$pmmh$proposal_scale_max <- pmmh$scale_max
  validate_config_null_lasso(null)
  null
}

# Selector and conditioning events

select_nonnull_lasso <- function(Y, W, X, fold_id, config, context = "selector") {
  # Warnings remain visible, with the dataset/method context when available.
  fit_context <- context
  forward_warning <- function(w) {
    message(fit_context, ": ", conditionMessage(w))
  }
  fit_context <- paste(context, "CV fit")
  fit <- withCallingHandlers(
    glmnet::cv.glmnet(
      x = X, y = as.numeric(Y) + as.numeric(W), alpha = 1,
      lambda = config$lambda_grid, foldid = fold_id, nfolds = config$folds,
      type.measure = "deviance", intercept = FALSE, standardize = FALSE
    ),
    warning = forward_warning
  )
  index <- as.integer(which.min(fit$cvm))
  if (abs(fit$lambda[index] - config$lambda_grid[index]) > 1e-10) {
    stop("CV changed the fixed lambda grid: ", context)
  }
  fit_context <- paste(context, "full-data refit")
  refit <- withCallingHandlers(
    glmnet::glmnet(
      x = X, y = as.numeric(Y) + as.numeric(W), alpha = 1,
      lambda = config$lambda_grid[index], intercept = FALSE,
      standardize = FALSE, thresh = config$glmnet_threshold
    ),
    warning = forward_warning
  )
  coefficients <- as.numeric(stats::coef(refit))[-1L]
  if (length(coefficients) != config$p || any(!is.finite(coefficients))) {
    stop("Invalid Lasso coefficients: ", context)
  }
  selected <- as.integer(which(abs(coefficients) > config$coefficient_tolerance))
  list(
    selected_model = selected, selected_signs = as.integer(sign(coefficients[selected])),
    lambda_index = index, lambda_observed = as.numeric(config$lambda_grid[index])
  )
}

selection_event_nonnull_lasso <- function(T, U, dataset, method_id, config) {
  if (!method_id %in% config$sampled_methods) stop("Unknown Lasso sampling method: ", method_id)
  Y <- as.numeric(dataset$target$C_M %*% as.numeric(T) + as.numeric(dataset$target$Z_obs))
  selected <- select_nonnull_lasso(
    Y, config$randomization_sd * as.numeric(U), dataset$X, dataset$fold_id, config,
    paste0("delta=", dataset$delta, ", dataset_id=", dataset$dataset_id, ", method=", method_id)
  )
  same_model <- identical(selected$selected_model, dataset$selected_model)
  if (method_id == "model_only") return(same_model)
  same_model && identical(selected$selected_signs, dataset$selected_signs) &&
    identical(selected$lambda_index, dataset$lambda_index)
}

# Conditional target, direction and analytic data splitting

build_target_nonnull_lasso <- function(dataset, config) {
  X_M <- dataset$X[, dataset$selected_model, drop = FALSE]
  if (qr(X_M, tol = config$rank_tolerance)$rank != dataset$d) {
    stop("Selected design is rank deficient.")
  }
  gram <- crossprod(X_M)
  inverse <- chol2inv(chol(gram))
  H_M <- inverse %*% t(X_M)
  T_obs <- as.numeric(H_M %*% dataset$Y)
  # Preserve the formal reconstruction path rather than replacing C_M by X_M.
  C_M <- t(H_M) %*% solve(H_M %*% t(H_M))
  Z_obs <- as.numeric(dataset$Y - C_M %*% T_obs)
  list(
    X_M = X_M, H_M = H_M, T_obs = T_obs,
    Sigma_M = config$response_sd^2 * inverse, C_M = C_M, Z_obs = Z_obs
  )
}

build_direction_nonnull_lasso <- function(dataset, config) {
  X_M <- dataset$target$X_M
  gram <- crossprod(X_M)
  rhs <- as.numeric(crossprod(X_M, dataset$X %*% config$signal_direction))
  factor <- chol(gram)
  a <- as.numeric(backsolve(factor, forwardsolve(t(factor), rhs)))
  precision_a <- as.numeric(gram %*% a)
  denominator <- sqrt(sum(a * precision_a))
  if (!is.finite(denominator) || denominator <= 1e-12) stop("Degenerate projected direction.")
  list(a = a, precision_a = precision_a, denominator = denominator)
}

statistic_nonnull_lasso <- function(T, direction) {
  if (is.null(dim(T))) return(sum(direction$precision_a * as.numeric(T)) / direction$denominator)
  as.numeric(as.matrix(T) %*% direction$precision_a) / direction$denominator
}

fission_statistic_nonnull_lasso <- function(dataset, config) {
  Y_fission <- as.numeric(dataset$Y - config$randomization_sd^(-2) * dataset$W)
  # The formal matched baseline solves the normal equations directly.
  gram <- crossprod(dataset$target$X_M)
  T_fission <- as.numeric(solve(gram, crossprod(dataset$target$X_M, Y_fission)))
  statistic_nonnull_lasso(T_fission, dataset$direction)
}

fission_pvalue_nonnull_lasso <- function(
    dataset, config, score = fission_statistic_nonnull_lasso(dataset, config)) {
  stats::pnorm(score / sqrt(1 + config$randomization_sd^(-2)), lower.tail = FALSE)
}

convert_null_nonnull_lasso <- function(dataset, config, selected_signs) {
  converted <- list(
    delta = 0, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
    X = dataset$design, Y = dataset$response, W = dataset$randomization,
    fold_id = dataset$fold_id, selected_model = dataset$selected_model,
    selected_signs = as.integer(selected_signs), lambda_index = dataset$lambda_index,
    lambda_observed = dataset$lambda_observed, d = dataset$target_dimension,
    generation_rng_state = dataset$generation_rng_state
  )
  converted$target <- build_target_nonnull_lasso(converted, config)
  converted$direction <- build_direction_nonnull_lasso(converted, config)
  converted
}

# Raw attempts, eligible bank and task streams

with_rng_nonnull_lasso <- function(fun, seed = NULL, state = NULL,
                                   kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection")) {
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
  do.call(RNGkind, as.list(kind))
  if (is.null(state)) set.seed(as.integer(seed))
  else assign(".Random.seed", as.integer(state), envir = .GlobalEnv)
  fun()
}

rng_streams_nonnull_lasso <- function(seed, count) {
  with_rng_nonnull_lasso(function() {
    streams <- vector("list", count)
    state <- get(".Random.seed", envir = .GlobalEnv)
    for (index in seq_len(count)) {
      streams[[index]] <- as.integer(state)
      state <- parallel::nextRNGStream(state)
    }
    streams
  }, seed = seed)
}

generate_attempt_nonnull_lasso <- function(delta, raw_id, state, config) {
  seeds <- with_rng_nonnull_lasso(function() sample.int(2^30, 4L), state = state)
  component <- function(index, fun) {
    with_rng_nonnull_lasso(
      fun, seed = seeds[index], kind = c("Mersenne-Twister", "Inversion", "Rejection")
    )
  }
  X <- component(1L, function() {
    innovations <- matrix(rnorm(config$n * config$p), config$n, config$p)
    design <- matrix(0, config$n, config$p)
    design[, 1L] <- innovations[, 1L]
    for (column in 2:config$p) {
      design[, column] <- config$design_correlation * design[, column - 1L] +
        sqrt(1 - config$design_correlation^2) * innovations[, column]
    }
    centered <- sweep(design, 2L, colMeans(design), "-")
    sweep(centered, 2L, sqrt(colSums(centered^2)) / sqrt(config$n), "/")
  })
  epsilon <- component(2L, function() rnorm(config$n, sd = config$response_sd))
  Y <- as.numeric(X %*% as.numeric(delta * config$signal_direction) + epsilon)
  W <- component(3L, function() rnorm(config$n, sd = config$randomization_sd))
  folds <- component(4L, function() sample(rep(seq_len(config$folds), length.out = config$n)))
  observed <- select_nonnull_lasso(Y, W, X, folds, config, paste0("raw_id=", raw_id))
  dataset <- c(list(
    delta = delta, raw_id = as.integer(raw_id), X = X, Y = Y, W = W,
    fold_id = as.integer(folds), generation_rng_state = as.integer(state), component_seeds = seeds
  ), observed)
  dataset$d <- as.integer(length(dataset$selected_model))
  if (!dataset$d %in% seq.int(config$dimension_range[1L], config$dimension_range[2L])) return(NULL)
  dataset$target <- build_target_nonnull_lasso(dataset, config)
  replay <- select_nonnull_lasso(Y, W, X, folds, config, paste0("raw_id=", raw_id, ", replay"))
  if (!identical(observed, replay)) stop("Observed selector replay failed.")
  dataset$direction <- build_direction_nonnull_lasso(dataset, config)
  dataset
}

assign_rng_nonnull_lasso <- function(datasets, config) {
  count <- sum(vapply(datasets, function(dataset) dataset$delta > 0, logical(1))) * length(config$sampled_methods)
  streams <- rng_streams_nonnull_lasso(
    config$rng$inference_master, config$rng$inference_stream_offset + count
  )
  rows <- list()
  positive_index <- 0L
  null_proposed <- rng_streams_nonnull_lasso(config$rng$null_proposed_master, 2L * config$datasets_per_delta)
  null_extra <- rng_streams_nonnull_lasso(config$rng$null_extra_master, config$datasets_per_delta)
  for (dataset in datasets) for (method_id in config$sampled_methods) {
    index <- length(rows) + 1L
    if (dataset$delta == 0) {
      stream_id <- if (method_id == "model_only") 2L * dataset$dataset_id - 1L else dataset$dataset_id
      state <- if (method_id == "model_only") null_proposed[[stream_id]] else null_extra[[stream_id]]
    } else {
      positive_index <- positive_index + 1L
      stream_id <- config$rng$inference_stream_offset + positive_index
      state <- streams[[stream_id]]
    }
    row <- data.frame(
      delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
      method_id = method_id, stream_id = stream_id, stringsAsFactors = FALSE
    )
    row$rng_state <- I(list(state))
    rows[[index]] <- row
  }
  do.call(rbind, rows)
}

prepare_bank_nonnull_lasso <- function(config = default_config_nonnull_lasso()) {
  validate_config_nonnull_lasso(config)
  datasets <- list()
  raw_attempts <- integer(length(config$delta_grid))
  for (group in seq_along(config$delta_grid)) {
    if (config$delta_grid[group] == 0) {
      null_config <- null_config_nonnull_lasso(config)
      null_bank <- prepare_bank_null_lasso(null_config)
      for (dataset in null_bank$datasets) {
        observed <- select_nonnull_lasso(dataset$response, dataset$randomization,
                                         dataset$design, dataset$fold_id, config)
        datasets[[length(datasets) + 1L]] <- convert_null_nonnull_lasso(dataset, config, observed$selected_signs)
      }
      raw_attempts[group] <- null_bank$raw_attempts
      next
    }
    positive_group <- match(config$delta_grid[group], config$delta_grid[config$delta_grid > 0])
    streams <- rng_streams_nonnull_lasso(config$rng$generation_seeds[positive_group], config$raw_attempt_limit)
    retained <- 0L
    for (raw_id in seq_len(config$raw_attempt_limit)) {
      dataset <- generate_attempt_nonnull_lasso(config$delta_grid[group], raw_id, streams[[raw_id]], config)
      if (is.null(dataset)) next
      retained <- retained + 1L
      dataset$dataset_id <- retained
      datasets[[length(datasets) + 1L]] <- dataset
      if (retained == config$datasets_per_delta) break
    }
    if (retained != config$datasets_per_delta) stop("Insufficient eligible datasets within the raw-attempt limit.")
    raw_attempts[group] <- raw_id
  }
  bank <- list(config = config, datasets = datasets,
               rng_assignments = assign_rng_nonnull_lasso(datasets, config), raw_attempts = raw_attempts)
  validate_bank_nonnull_lasso(bank)
  bank
}

# Summaries and paired-result checks

result_row_nonnull_lasso <- function(dataset, method_id, p_value, config,
                                         tail_count = NA_integer_,
                                         production_draws = NA_integer_,
                                         result_type = "analytic",
                                         selector_evaluations = 0L,
                                         statistic_observed = NULL) {
  if (is.null(statistic_observed)) {
    statistic_observed <- if (method_id == "data_fission") {
      fission_statistic_nonnull_lasso(dataset, config)
    } else {
      statistic_nonnull_lasso(dataset$target$T_obs, dataset$direction)
    }
  }
  data.frame(
    delta = dataset$delta, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
    method_id = method_id, target_id = "selected_model_projection",
    target_dimension = dataset$d,
    statistic_observed = statistic_observed,
    p_value = p_value, reject = as.integer(p_value <= config$test_level),
    tail_count = tail_count, production_draws = production_draws,
    result_type = result_type, status = "complete",
    selector_evaluations = selector_evaluations, stringsAsFactors = FALSE
  )
}

summarize_nonnull_lasso <- function(data, config) {
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

validate_bank_nonnull_lasso <- function(bank) {
  if (!all(c("config", "datasets", "rng_assignments") %in% names(bank))) {
    stop("The data bank is incomplete.")
  }
  config <- bank$config
  validate_config_nonnull_lasso(config)
  datasets <- bank$datasets
  if (length(datasets) != length(config$delta_grid) * config$datasets_per_delta) {
    stop("The data bank has the wrong number of datasets.")
  }
  identity <- do.call(rbind, lapply(datasets, function(dataset) {
    required <- c("delta", "dataset_id", "raw_id", "Y", "W", "d", "target", "direction")
    if (!all(required %in% names(dataset))) stop("A dataset is incomplete.")
    if (length(dataset$Y) != config$n || length(dataset$W) != config$n ||
        any(!is.finite(c(dataset$Y, dataset$W, dataset$target$T_obs)))) {
      stop("Dataset inputs or observed target values are invalid.")
    }
    data.frame(delta = dataset$delta, dataset_id = dataset$dataset_id,
               raw_id = dataset$raw_id, dimension = dataset$d)
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
    if (any(rows$dimension < config$dimension_range[1L] |
            rows$dimension > config$dimension_range[2L])) {
      stop("A retained dataset violates the eligibility rule.")
    }
  }
  for (dataset in datasets) {
    if (length(dataset$selected_model) != dataset$d ||
        length(dataset$selected_signs) != dataset$d ||
        !identical(dim(dataset$X), c(config$n, config$p)) ||
        !identical(as.integer(tabulate(dataset$fold_id, config$folds)),
                   rep(as.integer(config$n / config$folds), config$folds))) {
      stop("Selected model or fixed-fold dimensions are invalid.")
    }
    covariance <- dataset$target$Sigma_M
    required_target <- c("H_M", "T_obs", "C_M", "Z_obs", "Sigma_M")
    if (!all(required_target %in% names(dataset$target)) ||
        !identical(dim(covariance), rep(as.integer(dataset$d), 2L)) ||
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

validate_result_nonnull_lasso <- function(result) {
  validate_bank_nonnull_lasso(result)
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
  dimension <- vapply(result$datasets, `[[`, integer(1L), "d")
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
  if (anyNA(data$result_type) || any(data$result_type[!sampled] != "analytic")) {
    stop("Data Fission must be identified as an analytic result.")
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
  analytic <- data$method_id == "data_fission"
  expected_fission <- stats::pnorm(
    data$statistic_observed[analytic] / sqrt(1 + config$randomization_sd^(-2)), lower.tail = FALSE
  )
  if (any(!is.finite(data$statistic_observed[analytic])) ||
      any(abs(expected_fission - data$p_value[analytic]) > 1e-15)) {
    stop("Data Fission p-values do not match their observed statistic and null variance.")
  }
  reconstructed <- summarize_nonnull_lasso(data, config)$power
  if (!isTRUE(all.equal(reconstructed, result$summary$power, tolerance = 1e-15, check.attributes = FALSE))) {
    stop("The power summary does not match the paired dataset-level decisions.")
  }
  invisible(result)
}
