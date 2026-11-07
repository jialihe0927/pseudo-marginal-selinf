same_model_lasso <- function(a, b) {
  identical(sort(as.integer(a)), sort(as.integer(b)))
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

# Choose the first CV minimum on the fixed folds/grid, then refit on all data.
cv_lasso_selector <- function(response, randomization, selector) {
  cfg <- selector$config
  randomized_response <- as.numeric(response) + as.numeric(randomization)

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

  selection <- fixed_lasso_fit(randomized_response, lambda, selector)
  selection$lambda_index <- index
  selection$cvm <- as.numeric(fit$cvm)
  selection
}

# Reconstruct the response with the observed nuisance residual Z_obs fixed.
reconstruct_response_lasso <- function(target, benchmark) {
  as.numeric(benchmark$C %*% as.numeric(target) + benchmark$Z_obs)
}

# Condition only on the active set, not the CV penalty or coefficient signs.
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

