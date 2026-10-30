same_model_lasso <- function(a, b) {
  identical(sort(as.integer(a)), sort(as.integer(b)))
}

validate_efficiency_config_lasso <- function(config) {
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

mean_ess_lasso <- function(x, method, coordinate) {
  label <- paste(method, coordinate, sep = ", ")
  estimate <- mean_diagnostics(x, label)
  if (nzchar(estimate$undefined_reason)) {
    message(label, ": ", estimate$undefined_reason)
  }
  estimate$ess_mean
}

canonical_selection_lasso <- function(lambda_index, lambda_value, beta, tolerance) {
  active <- sort(as.integer(which(abs(beta) > tolerance)))
  list(
    lambda_index = as.integer(lambda_index),
    lambda_value = as.numeric(lambda_value),
    active = active,
    signs = as.integer(sign(beta[active])),
    beta = as.numeric(beta)
  )
}

fixed_lasso_fit <- function(response, lambda, selector) {
  cfg <- selector$config
  fit <- glmnet::glmnet(
    x = selector$X,
    y = as.numeric(response),
    alpha = cfg$alpha,
    lambda = as.numeric(lambda),
    intercept = cfg$intercept,
    standardize = cfg$standardize,
    thresh = cfg$glmnet_thresh
  )
  canonical_selection_lasso(
    NA_integer_,
    lambda,
    as.numeric(stats::coef(fit))[-1L],
    cfg$coef_tol
  )
}

cv_lasso_selector <- function(response, randomization, selector) {
  cfg <- selector$config
  randomized_response <- as.numeric(response) + as.numeric(randomization)
  if (length(randomized_response) != nrow(selector$X) ||
      any(!is.finite(randomized_response))) {
    stop("Invalid candidate response or randomization.")
  }

  fit <- glmnet::cv.glmnet(
    x = selector$X,
    y = randomized_response,
    alpha = cfg$alpha,
    lambda = cfg$lambda_grid,
    foldid = selector$foldid,
    nfolds = cfg$K,
    type.measure = cfg$type_measure,
    intercept = cfg$intercept,
    standardize = cfg$standardize
  )
  index <- as.integer(which.min(fit$cvm))
  lambda <- as.numeric(cfg$lambda_grid[index])
  tolerance <- cfg$lambda_grid_internal_tol * max(1, abs(lambda))
  if (abs(as.numeric(fit$lambda[index]) - lambda) > tolerance) {
    stop("cv.glmnet returned a lambda inconsistent with the frozen grid.")
  }

  selection <- fixed_lasso_fit(randomized_response, lambda, selector)
  selection$lambda_index <- index
  selection$cvm <- as.numeric(fit$cvm)
  selection
}

reconstruct_response_lasso <- function(target, benchmark) {
  as.numeric(benchmark$C %*% as.numeric(target) + benchmark$Z_obs)
}

selection_event_lasso <- function(target, auxiliary, benchmark) {
  response <- reconstruct_response_lasso(target, benchmark)
  selection <- cv_lasso_selector(
    response,
    benchmark$gamma * as.numeric(auxiliary),
    benchmark$selector
  )
  same_model_lasso(selection$active, benchmark$selected_model)
}

log_gaussian_lasso <- function(target, center, covariance_inverse, log_determinant) {
  delta <- as.numeric(target) - as.numeric(center)
  -0.5 * (
    length(delta) * log(2 * pi) +
      log_determinant +
      drop(crossprod(delta, covariance_inverse %*% delta))
  )
}

