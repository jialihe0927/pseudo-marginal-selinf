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
      # Reserve this stream prefix before assigning positive-signal tasks.
      inference_stream_offset = 400L
    )
  )
}

null_config_nonnull_lasso <- function(config) {
  # Construct the calibration configuration used to regenerate delta-zero data and chains.
  null <- default_config_null_lasso()
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
  null
}

# Run fixed-fold CV and the full-data refit; retain penalty index, model and active signs.
select_nonnull_lasso <- function(Y, W, X, fold_id, config) {
  fit <- glmnet::cv.glmnet(
    x = X, y = as.numeric(Y) + as.numeric(W), alpha = 1,
    lambda = config$lambda_grid, foldid = fold_id, nfolds = config$folds,
    type.measure = "deviance", intercept = FALSE, standardize = FALSE
  )
  index <- as.integer(which.min(fit$cvm))
  refit <- glmnet::glmnet(
    x = X, y = as.numeric(Y) + as.numeric(W), alpha = 1,
    lambda = config$lambda_grid[index], intercept = FALSE,
    standardize = FALSE, thresh = config$glmnet_threshold
  )
  coefficients <- as.numeric(stats::coef(refit))[-1L]
  selected <- as.integer(which(abs(coefficients) > config$coefficient_tolerance))
  list(
    selected_model = selected, selected_signs = as.integer(sign(coefficients[selected])),
    lambda_index = index, lambda_observed = as.numeric(config$lambda_grid[index])
  )
}

# Proposed fixes the model only; Extra conditioning also fixes the penalty index and signs.
selection_event_nonnull_lasso <- function(T, U, dataset, method_id, config) {
  Y <- as.numeric(dataset$target$C_M %*% as.numeric(T) + as.numeric(dataset$target$Z_obs))
  selected <- select_nonnull_lasso(
    Y, config$randomization_sd * as.numeric(U), dataset$X, dataset$fold_id, config
  )
  same_model <- identical(selected$selected_model, dataset$selected_model)
  if (method_id == "model_only") return(same_model)
  same_model && identical(selected$selected_signs, dataset$selected_signs) &&
    identical(selected$lambda_index, dataset$lambda_index)
}

build_target_nonnull_lasso <- function(dataset, config) {
  X_M <- dataset$X[, dataset$selected_model, drop = FALSE]
  gram <- crossprod(X_M)
  inverse <- chol2inv(chol(gram))
  H_M <- inverse %*% t(X_M)
  T_obs <- as.numeric(H_M %*% dataset$Y)
  # C_M maps target coefficients back to Y; Z_obs is held fixed during conditional sampling.
  C_M <- t(H_M) %*% solve(H_M %*% t(H_M))
  Z_obs <- as.numeric(dataset$Y - C_M %*% T_obs)
  list(
    X_M = X_M, H_M = H_M, T_obs = T_obs,
    Sigma_M = config$response_sd^2 * inverse, C_M = C_M, Z_obs = Z_obs
  )
}

# Project the fixed signal template onto the selected model; delta does not enter the direction.
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

# Standardize the fixed directional score for one target vector or rows of target draws.
statistic_nonnull_lasso <- function(T, direction) {
  if (is.null(dim(T))) return(sum(direction$precision_a * as.numeric(T)) / direction$denominator)
  as.numeric(as.matrix(T) %*% direction$precision_a) / direction$denominator
}

fission_statistic_nonnull_lasso <- function(dataset, config) {
  # This inference view is independent of the selection view Y + W under unit Gaussian noise.
  Y_fission <- as.numeric(dataset$Y - config$randomization_sd^(-2) * dataset$W)
  # Compute the selected-model coefficients of the independent inference view.
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
  # Separate seeds drive the design, response noise, randomization and fold assignment.
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
  observed <- select_nonnull_lasso(Y, W, X, folds, config)
  dataset <- c(list(
    delta = delta, raw_id = as.integer(raw_id), X = X, Y = Y, W = W,
    fold_id = as.integer(folds), generation_rng_state = as.integer(state), component_seeds = seeds
  ), observed)
  dataset$d <- as.integer(length(dataset$selected_model))
  if (!dataset$d %in% seq.int(config$dimension_range[1L], config$dimension_range[2L])) return(NULL)
  dataset$target <- build_target_nonnull_lasso(dataset, config)
  dataset$direction <- build_direction_nonnull_lasso(dataset, config)
  dataset
}

# Pair streams with delta/dataset/method identities in the fixed bank order.
# Delta-zero Proposed keeps the calibration stream; Extra conditioning has its own stream.
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

# Regenerate delta zero through the calibration generator; retain other signals in raw order.
prepare_bank_nonnull_lasso <- function(config = default_config_nonnull_lasso()) {
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
  bank
}

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
