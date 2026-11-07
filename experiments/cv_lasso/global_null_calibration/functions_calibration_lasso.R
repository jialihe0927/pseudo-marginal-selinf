default_config_null_lasso <- function() {
  list(
    experiment = "CV-Lasso global-null calibration",
    n = 200L,
    p = 100L,
    design_correlation = 0.5,
    response_sd = 1,
    randomization_sd = 0.15,
    folds = 10L,
    lambda_grid = c(
      0.230497466833294, 0.218176011289412, 0.206513210562027,
      0.195473855647968, 0.185024619673937, 0.175133957285509,
      0.165772009414394, 0.156910513136447, 0.148522716348308,
      0.140583297005077, 0.133068286675210, 0.125954998181871,
      0.119221957112270, 0.112848836988239, 0.106816397902325,
      0.101106428434145, 0.0957016906716557, 0.0905858681713679,
      0.0857435167003877, 0.0811600176115930, 0.0768215337111816,
      0.0727149674853597, 0.0688279215600626, 0.0651486612743371,
      0.0616660792543981, 0.0583696618814121, 0.0552494575517768,
      0.0522960466340762, 0.0495005130320147, 0.0468544172674812,
      0.0443497710024822, 0.0419790129230295
    ),
    glmnet_threshold = 1e-7,
    coefficient_tolerance = 1e-8,
    lambda_tolerance = 1e-10,
    dimension_range = c(1L, 10L),
    dataset_count = 200L,
    raw_attempt_limit = 5000L,
    retention_rule = "first eligible datasets in raw-attempt order",
    methods = data.frame(
      method_id = c("model_only", "fixed_lambda", "naive_gaussian"),
      display_name = c("Model-only", "Fixed-lambda", "Naive"),
      event = c(
        "rerun fixed-fold CV and full-data refit; selected model equals observed model",
        "fix observed lambda and rerun full-data refit; selected model equals observed model",
        "no selection adjustment; analytic chi-square reference"
      ),
      stringsAsFactors = FALSE
    ),
    sampled_method_ids = c("model_only", "fixed_lambda"),
    statistic = "Q = T' Sigma_T^{-1} T",
    p_value = "inclusive upper-tail count / 5000; no +1 correction",
    test_level = 0.05,
    rejection_rule = "p <= 0.05",
    pmmh = list(
      warmup = 200L, stabilization = 0L, production = 5000L,
      total_transitions = 5200L, auxiliary_samples = 1L, rho = 0.8,
      proposal_scale_initial = 0.8, proposal_scale_min = 0.05,
      proposal_scale_max = 5, target_acceptance = 0.20,
      adaptation = "h <- exp(log(h) + 0.5/(iteration+10)^0.6 * (acceptance_probability - 0.20))",
      production_adaptation = FALSE,
      proposal = "T' = T + h/sqrt(d) G xi, with G G' = Sigma_T"
    ),
    rng = list(
      kind = "L'Ecuyer-CMRG", normal_kind = "Inversion", sample_kind = "Rejection",
      raw_master_seed = 2026999901L, formal_master_seed = 2026999902L,
      assignment = "one preallocated L'Ecuyer-CMRG stream per dataset-method job"
    )
  )
}

selector_setup_null_lasso <- function(dataset) {
  list(
    design = dataset$design,
    fold_id = dataset$fold_id,
    lambda_grid = as.numeric(dataset$lambda_grid)
  )
}

fit_fixed_lambda_null_lasso <- function(response, lambda, setup, config) {
  fit <- glmnet::glmnet(
    x = setup$design,
    y = as.numeric(response),
    alpha = 1,
    lambda = as.numeric(lambda),
    intercept = FALSE,
    standardize = FALSE,
    thresh = config$glmnet_threshold
  )
  coefficients <- as.numeric(stats::coef(fit))[-1L]
  selected <- which(abs(coefficients) > config$coefficient_tolerance)
  list(selected = as.integer(selected), coefficients = coefficients)
}

# Select the first CV minimum with fixed folds, then refit at that penalty on all data.
select_cv_lasso_null <- function(response, randomization, setup, config) {
  selection_response <- as.numeric(response) + as.numeric(randomization)
  fit <- glmnet::cv.glmnet(
    x = setup$design,
    y = selection_response,
    alpha = 1,
    lambda = setup$lambda_grid,
    foldid = setup$fold_id,
    nfolds = config$folds,
    type.measure = "deviance",
    intercept = FALSE,
    standardize = FALSE
  )
  lambda_index <- as.integer(which.min(as.numeric(fit$cvm)))
  refit <- fit_fixed_lambda_null_lasso(
    selection_response, setup$lambda_grid[lambda_index], setup, config
  )
  list(selected = refit$selected, lambda_index = lambda_index)
}

select_fixed_lambda_null_lasso <- function(response, randomization, lambda, setup, config) {
  fit_fixed_lambda_null_lasso(
    as.numeric(response) + as.numeric(randomization), lambda, setup, config
  )
}

same_model_null_lasso <- function(left, right) {
  identical(sort(as.integer(left)), sort(as.integer(right)))
}

# Proposed reruns CV and matches only the model; Fixed-lambda holds the penalty fixed.
# Neither event additionally fixes coefficient signs.
selection_event_null_lasso <- function(target_value, auxiliary, dataset, method_id, config) {
  response <- reconstruct_response_null_lasso(target_value, dataset$target)
  randomization <- config$randomization_sd * as.numeric(auxiliary)
  setup <- selector_setup_null_lasso(dataset)

  selected <- if (identical(method_id, "model_only")) {
    select_cv_lasso_null(response, randomization, setup, config)$selected
  } else {
    select_fixed_lambda_null_lasso(
      response,
      randomization,
      dataset$lambda_observed,
      setup,
      config
    )$selected
  }

  same_model_null_lasso(selected, dataset$selected_model)
}

# Selected-model coefficients T = HY, their covariance, and the residual Y - X_M T.
build_target_null_lasso <- function(dataset, config) {
  selected <- sort(as.integer(dataset$selected_model))
  dimension <- length(selected)
  X_M <- dataset$design[, selected, drop = FALSE]
  G_M <- crossprod(X_M)
  Sigma_M <- chol2inv(chol(G_M))
  H <- Sigma_M %*% t(X_M)
  T_obs <- as.numeric(H %*% dataset$response)
  Z_obs <- as.numeric(dataset$response - X_M %*% T_obs)
  Q_obs <- as.numeric(crossprod(T_obs, G_M %*% T_obs))
  list(
    dataset_id = dataset$dataset_id, raw_id = dataset$raw_id, M_obs = selected,
    d = dimension, X_M = X_M, G_M = G_M, Sigma_M = Sigma_M, H = H,
    C = X_M, T_obs = T_obs, Z_obs = Z_obs, Q_obs = Q_obs,
    lambda_obs = dataset$lambda_observed, lambda_index = dataset$lambda_index,
    foldid = dataset$fold_id, lambda_grid = dataset$lambda_grid
  )
}

# Vary T while holding the observed nuisance residual Z_obs fixed.
reconstruct_response_null_lasso <- function(target_value, target) {
  as.numeric(target$X_M %*% as.numeric(target_value) + target$Z_obs)
}

target_statistic_null_lasso <- function(target_value, target) {
  # Q = T^T Sigma_M^(-1) T, evaluated as the squared norm of whitened T.
  factor <- chol(target$Sigma_M)
  standardized <- backsolve(factor, as.numeric(target_value), transpose = TRUE)
  sum(standardized^2)
}

standardize_design_null_lasso <- function(design) {
  centered <- sweep(design, 2L, colMeans(design), "-")
  scale <- sqrt(colSums(centered^2)) / sqrt(nrow(centered))
  sweep(centered, 2L, scale, "/")
}

# Generate an AR(1) design with centered/scaled columns and unit-noise null responses.
generate_global_null_lasso <- function(config) {
  innovations <- matrix(rnorm(config$n * config$p), config$n, config$p)
  design_raw <- matrix(0, config$n, config$p)
  design_raw[, 1L] <- innovations[, 1L]

  if (config$p > 1L) {
    for (column in 2:config$p) {
      design_raw[, column] <-
        config$design_correlation * design_raw[, column - 1L] +
        sqrt(1 - config$design_correlation^2) * innovations[, column]
    }
  }

  response <- rnorm(config$n)
  randomization <- rnorm(config$n, sd = config$randomization_sd)
  fold_id <- sample(rep(seq_len(config$folds), length.out = config$n))

  list(
    design = standardize_design_null_lasso(design_raw),
    response = response,
    randomization = randomization,
    fold_id = as.integer(fold_id)
  )
}

rng_streams_null_lasso <- function(master_seed, n_streams) {
  if (n_streams < 1L) return(list())
  old_kind <- RNGkind()
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(as.integer(master_seed))
  state <- get(".Random.seed", envir = .GlobalEnv)
  streams <- vector("list", n_streams)
  for (index in seq_len(n_streams)) {
    streams[[index]] <- as.integer(state)
    state <- parallel::nextRNGStream(state)
  }
  streams
}

with_rng_state_null_lasso <- function(state, fun) {
  old_kind <- RNGkind()
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", as.integer(state), envir = .GlobalEnv)
  value <- fun()
  list(value = value, final_state = as.integer(get(".Random.seed", envir = .GlobalEnv)))
}

generate_raw_attempt_null_lasso <- function(raw_id, rng_state, config) {
  generated <- with_rng_state_null_lasso(rng_state, function() generate_global_null_lasso(config))
  dataset <- c(
    list(raw_id = as.integer(raw_id), generation_stream_id = as.integer(raw_id),
         generation_rng_state = as.integer(rng_state)),
    generated$value,
    list(lambda_grid = config$lambda_grid)
  )
  setup <- selector_setup_null_lasso(dataset)
  observed <- select_cv_lasso_null(
    dataset$response, dataset$randomization, setup, config
  )
  c(dataset, list(
    lambda_observed = dataset$lambda_grid[observed$lambda_index],
    lambda_index = observed$lambda_index,
    selected_model = observed$selected,
    target_dimension = length(observed$selected)
  ))
}

# Preallocate streams in retained dataset/method order before any tasks are run.
assign_rng_null_lasso <- function(datasets, config) {
  sampled <- config$sampled_method_ids
  streams <- rng_streams_null_lasso(config$rng$formal_master_seed, length(datasets) * length(sampled))
  records <- vector("list", length(streams))
  job <- 0L
  for (dataset in datasets) for (method_id in sampled) {
    job <- job + 1L
    records[[job]] <- data.frame(
      job_id = job, dataset_id = dataset$dataset_id, raw_id = dataset$raw_id,
      method_id = method_id, rng_stream_id = job
    )
    records[[job]]$rng_state <- I(list(streams[[job]]))
  }
  do.call(rbind, records)
}

# Keep the first raw attempts whose selected model satisfies the dimension rule.
prepare_bank_null_lasso <- function(config = default_config_null_lasso()) {
  streams <- rng_streams_null_lasso(config$rng$raw_master_seed, config$raw_attempt_limit)
  datasets <- list()
  raw_id <- 0L
  while (length(datasets) < config$dataset_count && raw_id < config$raw_attempt_limit) {
    raw_id <- raw_id + 1L
    candidate <- generate_raw_attempt_null_lasso(raw_id, streams[[raw_id]], config)
    if (candidate$target_dimension < config$dimension_range[1L] ||
        candidate$target_dimension > config$dimension_range[2L]) next
    candidate$dataset_id <- length(datasets) + 1L
    candidate$target <- build_target_null_lasso(candidate, config)
    datasets[[candidate$dataset_id]] <- candidate
  }
  if (length(datasets) != config$dataset_count) stop("The raw-attempt limit was reached before the bank was complete.")
  config$raw_attempts <- raw_id
  bank <- list(config = config, datasets = datasets, rng_assignments = assign_rng_null_lasso(datasets, config), raw_attempts = raw_id)
  bank
}

summarize_null_results_lasso <- function(dataset_results, alpha = 0.05, methods = NULL) {
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
