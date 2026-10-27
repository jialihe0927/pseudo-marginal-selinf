validation_methods_lasso <- function() {
  c("Rejection", "Crude PMMH", "Correlated PMMH")
}

validate_samples_lasso <- function(object) {
  required <- c("config", "samples", "sample_info")
  if (!all(required %in% names(object))) stop("The validation object is incomplete.")
  config <- object$config
  if (!identical(config$application, "cv_lasso") ||
      !identical(config$coordinates, paste0("T", seq_len(6L)))) {
    stop("The application or coordinate order is incompatible.")
  }
  methods <- validation_methods_lasso()
  if (!identical(names(object$samples), methods) ||
      !all(methods %in% names(object$sample_info))) {
    stop("Rejection, PM and CPM samples and metadata are all required.")
  }
  for (method in methods) {
    samples <- object$samples[[method]]
    info <- object$sample_info[[method]]
    if (!is.matrix(samples) || !is.numeric(samples) ||
        !identical(colnames(samples), config$coordinates) ||
        ncol(samples) != length(config$coordinates) || nrow(samples) < 2L ||
        any(!is.finite(samples)) || nrow(samples) != info$sample_count ||
        !identical(info$target_id, config$target_id)) {
      stop("Invalid sample dimensions, coordinate order, target or values for ", method)
    }
    if (method == "Rejection") {
      if (info$selector_evaluations != config$sampler$rejection_attempts ||
          info$selector_evaluations != 3000000L) {
        stop("The reference budget must be 3,000,000 proposals.")
      }
      if (!is.null(info$accepted_identity)) {
        ids <- info$accepted_identity$attempt_id
        if (length(ids) != nrow(samples) || anyNA(ids) || anyDuplicated(ids) ||
            any(diff(ids) <= 0L) || any(ids < 1L | ids > info$selector_evaluations)) {
          stop("Invalid ordered rejection attempt identities.")
        }
      }
    } else {
      expected <- seq.int(config$sampler$adaptive_transitions + 1L,
                          config$sampler$total_transitions)
      if (nrow(samples) != config$sampler$production_transitions ||
          !identical(as.integer(info$production_iterations), expected) ||
          config$sampler$stabilization_transitions != 0L ||
          length(info$accepted) != nrow(samples) ||
          length(info$proposed_zero) != nrow(samples) ||
          anyNA(info$accepted) || anyNA(info$proposed_zero) ||
          info$selector_evaluations != config$sampler$total_transitions) {
        stop("Production interval or transition metadata is inconsistent for ", method)
      }
      repeats <- sum(rowSums(samples[-1L, , drop = FALSE] !=
                               samples[-nrow(samples), , drop = FALSE]) == 0L)
      if (repeats != info$adjacent_repeats) stop("Repeated production states changed for ", method)
      # A rejected proposal after the first retained state must preserve its predecessor.
      rejected <- which(!as.logical(info$accepted[-1L])) + 1L
      if (length(rejected) &&
          any(samples[rejected, , drop = FALSE] != samples[rejected - 1L, , drop = FALSE])) {
        stop("Rejected production transitions did not retain their states for ", method)
      }
    }
  }
  covariance <- config$target$covariance
  if (!identical(dim(covariance), rep(6L, 2L)) ||
      any(!is.finite(covariance)) || any(eigen(covariance, symmetric = TRUE)$values <= 0)) {
    stop("The frozen Gaussian covariance is invalid.")
  }
  if (max(abs(covariance - config$benchmark$Sigma_T)) > 1e-9 ||
      max(abs(config$target$gaussian_center - config$benchmark$gaussian_center)) > 1e-9) {
    stop("The diagnostic target and sampler target differ.")
  }
  invisible(object)
}

statistics_validation_lasso <- function(samples, config) {
  covariance <- config$target$covariance
  quadratic <- rowSums(samples * t(solve(covariance, t(samples))))
  directional <- as.numeric(solve(covariance, config$target$direction) %*% t(samples)) /
    sqrt(config$target$direction_information)
  data.frame(samples, Q = quadratic, S_proj = directional, check.names = FALSE)
}

ess_validation_lasso <- function(x, label = "statistic") {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for mean ESS and MCSE.")
  }
  x <- as.numeric(x)
  n <- length(x)
  undefined <- function(reason, package_warnings = "") {
    list(ess_mean = NA_real_, mcse_mean = NA_real_, tau = NA_real_,
         undefined_reason = reason, package_warnings = package_warnings)
  }
  if (any(!is.finite(x))) stop("Nonfinite ESS/MCSE input for ", label)
  if (n < 4L) return(undefined("insufficient observations"))
  if (length(unique(x)) < 2L) return(undefined("constant sequence"))

  # One variable, one chain: rows are iterations, not coordinates or methods.
  draws <- matrix(as.numeric(x), ncol = 1L)
  warning_messages <- character()
  estimates <- withCallingHandlers({
    list(
      ess_mean = as.numeric(posterior::ess_basic(draws, split = TRUE)),
      mcse_mean = as.numeric(posterior::mcse_mean(draws))
    )
  }, warning = function(condition) {
    warning_messages <<- c(warning_messages, conditionMessage(condition))
    message("posterior warning for ", label, ": ", conditionMessage(condition))
    # Do not muffle the package warning or catch unexpected errors.
  })
  package_warnings <- paste(unique(warning_messages), collapse = " | ")
  if (length(estimates$ess_mean) != 1L || length(estimates$mcse_mean) != 1L) {
    stop("Unexpected posterior result shape for ", label)
  }
  if (!is.finite(estimates$ess_mean) || estimates$ess_mean <= 0 ||
      !is.finite(estimates$mcse_mean) || estimates$mcse_mean <= 0) {
    return(undefined("posterior returned a nonfinite or nonpositive estimate", package_warnings))
  }
  expected_mcse <- stats::sd(x) / sqrt(estimates$ess_mean)
  if (!isTRUE(all.equal(estimates$mcse_mean, expected_mcse, tolerance = 1e-10))) {
    stop("posterior ESS/MCSE identity failed for ", label,
         "; check the installed posterior interface.")
  }
  list(
    ess_mean = estimates$ess_mean, mcse_mean = estimates$mcse_mean,
    tau = n / estimates$ess_mean,
    undefined_reason = "", package_warnings = package_warnings
  )
}

analyze_validation_lasso <- function(object) {
  validate_samples_lasso(object)
  config <- object$config
  config$diagnostic$ess <- "posterior::ess_basic(matrix(x, ncol = 1L), split = TRUE)"
  config$diagnostic$mean_mcse <- "posterior::mcse_mean(matrix(x, ncol = 1L))"
  config$diagnostic$indicator_mcse <- "split basic ESS and mean MCSE on each specified upper-tail indicator"
  config$diagnostic$holding_run <- paste(
    "Maximum number of consecutive exactly identical full state vectors within retained production;",
    "direct numeric comparison, including the first state of each run."
  )
  config$diagnostic$reference_independence <- paste(
    "Conditional on the frozen target, rejection and the two kernels use separately assigned",
    "RNG streams with no saved starting-state reuse; combined variances assume their independence."
  )
  object$config <- config
  object$config$software$diagnostics <- c(
    R = as.character(getRversion()), posterior = as.character(utils::packageVersion("posterior"))
  )
  methods <- validation_methods_lasso()
  values <- lapply(object$samples, statistics_validation_lasso, config = config)
  statistics <- colnames(values[[1L]])
  distribution <- list()
  ess_rows <- list()
  agreement <- list()
  for (method in methods) for (statistic in statistics) {
    x <- values[[method]][, statistic]
    quantiles <- stats::quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9), names = FALSE)
    distribution[[length(distribution) + 1L]] <- data.frame(
      method = method, statistic = statistic, n = length(x), mean = mean(x),
      sd = stats::sd(x), q10 = quantiles[1L], q25 = quantiles[2L],
      median = quantiles[3L], q75 = quantiles[4L], q90 = quantiles[5L]
    )
    if (method != "Rejection") {
      estimate <- ess_validation_lasso(x, paste(method, statistic, sep = ", "))
      half <- floor(length(x) / 2L)
      ess_rows[[length(ess_rows) + 1L]] <- data.frame(
        method = method, statistic = statistic, n = length(x),
        ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean, tau = estimate$tau,
        undefined_reason = estimate$undefined_reason, package_warnings = estimate$package_warnings,
        first_half_mean = mean(x[seq_len(half)]),
        second_half_mean = mean(x[seq.int(half + 1L, length(x))])
      )
    }
  }
  ess <- do.call(rbind, ess_rows)
  for (statistic in c(config$coordinates, "Q", "S_proj")) {
    reference <- values$Rejection[, statistic]
    reference_se <- stats::sd(reference) / sqrt(length(reference))
    for (method in methods[-1L]) {
      x <- values[[method]][, statistic]
      mcse <- ess$mcse_mean[ess$method == method & ess$statistic == statistic]
      combined <- sqrt(mcse^2 + reference_se^2)
      grid <- sort(unique(c(reference, x)))
      agreement[[length(agreement) + 1L]] <- data.frame(
        statistic = statistic, method = method, reference_mean = mean(reference),
        mean = mean(x), mean_difference = mean(x) - mean(reference),
        reference_se = reference_se, mcse_mean = mcse, combined_mcse = combined,
        z_difference = (mean(x) - mean(reference)) / combined,
        sd_ratio = stats::sd(x) / stats::sd(reference),
        median_difference = stats::median(x) - stats::median(reference),
        q10_difference = unname(stats::quantile(x, 0.1) - stats::quantile(reference, 0.1)),
        q90_difference = unname(stats::quantile(x, 0.9) - stats::quantile(reference, 0.9)),
        ecdf_difference = max(abs(stats::ecdf(reference)(grid) - stats::ecdf(x)(grid)))
      )
    }
  }
  agreement <- do.call(rbind, agreement)
  tails <- tails_validation_lasso(values, config)
  covariance <- covariance_validation_lasso(object$samples, config)
  diagnostics <- list(ess_mcse = ess, mixing = mixing_validation_lasso(object, ess))
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Package 'posterior' is required.")
  bulk <- list()
  for (method in methods[-1L]) for (coordinate in config$coordinates) {
    x <- object$samples[[method]][, coordinate]
    if (length(unique(x)) < 2L) {
      message("Bulk ESS undefined for ", method, ", ", coordinate, ": constant sequence")
      estimate <- NA_real_
    } else estimate <- as.numeric(posterior::ess_bulk(matrix(x, ncol = 1L)))[1L]
    bulk[[length(bulk) + 1L]] <- data.frame(method = method, coordinate = coordinate, ess_bulk = estimate)
  }
  diagnostics$coordinate_bulk_ess <- do.call(rbind, bulk)
  diagnostics$efficiency_estimators <- do.call(rbind, lapply(methods[-1L], function(method) {
    bulk_minimum <- min(diagnostics$coordinate_bulk_ess$ess_bulk[
      diagnostics$coordinate_bulk_ess$method == method
    ])
    mean_minimum <- min(ess$ess_mean[ess$method == method & ess$statistic %in% config$coordinates])
    evaluations <- object$sample_info[[method]]$selector_evaluations
    data.frame(
      method = method, posterior_bulk_minimum = bulk_minimum,
      posterior_bulk_per_1000 = 1000 * bulk_minimum / evaluations,
      ess_mean_minimum = mean_minimum, ess_mean_per_1000 = 1000 * mean_minimum / evaluations,
      reference_per_1000 = 1000 * nrow(object$samples$Rejection) /
        object$sample_info$Rejection$selector_evaluations
    )
  }))
  diagnostics$distribution <- do.call(rbind, distribution)
  diagnostics$agreement <- agreement
  diagnostics$tails <- tails
  diagnostics$covariance <- covariance
  object$diagnostics <- diagnostics
  object$summary <- list(
    undefined_mean_diagnostics = sum(!is.finite(ess$ess_mean)),
    undefined_tail_indicators = sum(tails$method != "Rejection" & !is.finite(tails$ess_indicator)),
    maximum_absolute_mean_z = if (all(is.finite(agreement$z_difference))) {
      max(abs(agreement$z_difference))
    } else NA_real_,
    coordinate_minimum_ess_mean = setNames(vapply(methods[-1L], function(method) {
      min(ess$ess_mean[ess$method == method & ess$statistic %in% config$coordinates])
    }, numeric(1L)), methods[-1L])
  )
  object
}

maximum_holding_run_lasso <- function(draws) {
  if (!is.matrix(draws) || !is.numeric(draws) ||
      nrow(draws) < 1L || ncol(draws) < 1L ||
      any(!is.finite(draws))) {
    stop("Holding-run calculation requires a finite numeric draw matrix.")
  }
  n <- nrow(draws)
  if (n == 1L) return(1L)

  # Count saved production states in each exactly identical full-vector run.
  same_as_previous <- rowSums(
    draws[-1L, , drop = FALSE] != draws[-n, , drop = FALSE]
  ) == 0L
  run_starts <- c(1L, which(!same_as_previous) + 1L)
  as.integer(max(diff(c(run_starts, n + 1L))))
}

mixing_validation_lasso <- function(object, ess) {
  rows <- list()
  for (method in validation_methods_lasso()[-1L]) {
    info <- object$sample_info[[method]]
    estimates <- ess[ess$method == method, , drop = FALSE]
    holding <- maximum_holding_run_lasso(object$samples[[method]])
    coordinates <- estimates[estimates$statistic %in% object$config$coordinates, ]
    minimum <- min(coordinates$ess_mean)
    minimum_coordinate <- if (is.finite(minimum)) {
      coordinates$statistic[which.min(coordinates$ess_mean)]
    } else NA_character_
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, n = nrow(object$samples[[method]]),
      selector_evaluations = info$selector_evaluations,
      acceptance = mean(info$accepted), zero_estimate_rate = mean(info$proposed_zero),
      final_scale = info$final_scale, maximum_holding_run = holding,
      minimum_coordinate_ess_mean = minimum, minimum_coordinate = minimum_coordinate,
      ess_q_mean = estimates$ess_mean[estimates$statistic == "Q"],
      ess_direction_mean = estimates$ess_mean[estimates$statistic == "S_proj"]
    )
  }
  do.call(rbind, rows)
}

covariance_validation_lasso <- function(samples, config) {
  covariances <- lapply(samples, stats::cov)
  correlations <- lapply(samples, stats::cor)
  decomposition <- eigen(config$target$covariance, symmetric = TRUE)
  whitening <- decomposition$vectors %*% diag(1 / sqrt(decomposition$values)) %*%
    t(decomposition$vectors)
  rows <- list()
  for (method in validation_methods_lasso()[-1L]) {
    difference <- covariances[[method]] - covariances$Rejection
    correlation_difference <- correlations[[method]] - correlations$Rejection
    whitened_difference <- whitening %*% difference %*% whitening
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, covariance_frobenius = sqrt(sum(difference^2)),
      relative_frobenius = sqrt(sum(difference^2)) / sqrt(sum(covariances$Rejection^2)),
      correlation_maximum = max(abs(correlation_difference[upper.tri(correlation_difference)])),
      whitened_frobenius = sqrt(sum(whitened_difference^2))
    )
  }
  list(comparison = do.call(rbind, rows), covariance = covariances,
       correlation = correlations, eigenvalues = lapply(covariances, function(x) {
         eigen(x, symmetric = TRUE, only.values = TRUE)$values
       }))
}

tails_validation_lasso <- function(values, config) {
  rows <- list()
  for (statistic in c("Q", "S_proj")) for (level in config$diagnostic$tail_levels) {
    threshold <- stats::quantile(values$Rejection[, statistic], level)
    reference <- mean(values$Rejection[, statistic] >= threshold)
    for (method in names(values)) {
      indicator <- as.numeric(values[[method]][, statistic] >= threshold)
      # Chain MCSE is conditional on this empirical threshold; it excludes its uncertainty.
      estimate <- if (method == "Rejection") {
        list(ess_mean = NA_real_, mcse_mean = NA_real_,
             undefined_reason = "threshold estimated from reference samples; reference MCSE not computed",
             package_warnings = "")
      } else ess_validation_lasso(indicator, paste(method, statistic, level, "upper tail"))
      rows[[length(rows) + 1L]] <- data.frame(
        statistic = statistic, level = level, threshold = as.numeric(threshold), method = method,
        probability = mean(indicator), reference_probability = reference,
        difference = mean(indicator) - reference, ess_indicator = estimate$ess_mean,
        mcse_indicator = estimate$mcse_mean, undefined_reason = estimate$undefined_reason,
        package_warnings = estimate$package_warnings, n = length(indicator),
        successes = sum(indicator), failures = sum(1 - indicator),
        transitions = sum(diff(indicator) != 0)
      )
    }
  }
  do.call(rbind, rows)
}

regenerate_validation_lasso <- function(object, efficiency_directory) {
  # Explicitly requested sampling reuses the released kernels; offline analysis never calls this.
  config <- object$config
  sampler_config <- config$sampler
  sampler_environment <- new.env(parent = environment())
  sys.source(file.path(efficiency_directory, "functions_lasso.R"), envir = sampler_environment)
  sys.source(file.path(efficiency_directory, "samplers_lasso.R"), envir = sampler_environment)
  sampler_environment$validate_efficiency_config_lasso(sampler_config)
  if (!requireNamespace("glmnet", quietly = TRUE) || !requireNamespace("posterior", quietly = TRUE)) {
    stop("Packages glmnet and posterior are required for explicit sample generation.")
  }
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  samples <- list()
  info <- list()
  rejection <- sampler_environment$run_rejection_lasso(config$benchmark, sampler_config)
  samples$Rejection <- rejection$accepted
  colnames(samples$Rejection) <- config$coordinates
  info$Rejection <- list(
    sample_count = nrow(rejection$accepted), selector_evaluations = rejection$attempts,
    target_id = config$target_id, rng_kind = RNGkind(),
    source = "explicit regeneration with the released rejection kernel",
    accepted_identity = NULL
  )
  for (method in validation_methods_lasso()[-1L]) {
    run <- sampler_environment$run_pmmh_lasso(config$benchmark, sampler_config, method)
    samples[[method]] <- run$draws
    colnames(samples[[method]]) <- config$coordinates
    repeats <- sum(rowSums(run$draws[-1L, , drop = FALSE] !=
                             run$draws[-nrow(run$draws), , drop = FALSE]) == 0L)
    info[[method]] <- list(
      sample_count = nrow(run$draws), selector_evaluations = run$selector_evaluations,
      target_id = config$target_id,
      production_iterations = seq.int(sampler_config$adaptive_transitions + 1L,
                                      sampler_config$total_transitions),
      accepted = run$accepted, proposed_zero = !run$proposed_selected, held = !run$accepted,
      adjacent_repeats = repeats, final_scale = run$final_proposal_scale,
      initial_target = config$benchmark$T_obs,
      initial_auxiliary = config$benchmark$W_obs / config$benchmark$gamma,
      source = "explicit regeneration with the released PMMH kernel"
    )
  }
  object$samples <- samples
  object$sample_info <- info
  object$diagnostics <- NULL
  object$summary <- NULL
  validate_samples_lasso(object)
  analyze_validation_lasso(object)
}
