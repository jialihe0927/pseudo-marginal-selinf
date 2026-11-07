# Keep raw coordinates and uncentered Q; S_shift and Q_centered use the frozen Gaussian center.
statistics_validation_cp <- function(samples, config) {
  inverse <- config$target$precision
  quadratic <- function(x) rowSums((x %*% inverse) * x)
  # Subtract each Gaussian-center coordinate down its sample column.
  centered <- sweep(samples, 2L, config$target$gaussian_center, "-")
  directional <- drop(centered %*% inverse %*% config$target$direction /
                        sqrt(config$target$direction_information))
  cbind(samples, Q = quadratic(samples), S_shift = directional,
        Q_centered = quadratic(centered))
}

# Generate the rejection reference and production chains with their saved sampling configuration.
regenerate_validation_cp <- function(object, efficiency_directory) {
  config <- object$config
  sampler_config <- config$sampler
  sampler_environment <- new.env(parent = environment())
  sys.source(file.path(efficiency_directory, "functions_cp.R"), envir = sampler_environment)
  sys.source(file.path(efficiency_directory, "samplers_cp.R"), envir = sampler_environment)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  samples <- list()
  info <- list()
  rejection <- sampler_environment$run_rejection_cp(config$benchmark, sampler_config)
  samples$Rejection <- rejection$accepted
  colnames(samples$Rejection) <- config$coordinates
  info$Rejection <- list(
    sample_count = nrow(rejection$accepted), selector_evaluations = rejection$attempts,
    target_id = config$target_id, rng_kind = RNGkind(),
    source = "explicit regeneration with the released rejection kernel",
    accepted_identity = NULL
  )
  for (method in c("Crude PMMH", "Correlated PMMH")) {
    run <- sampler_environment$run_pmmh_cp(config$benchmark, sampler_config, method)
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
      initial_auxiliary = config$benchmark$U_obs,
      source = "explicit regeneration with the released PMMH kernel"
    )
  }
  object$samples <- samples
  object$sample_info <- info
  object$diagnostics <- NULL
  object$summary <- NULL
  analyze_validation_cp(object)
}

autocovariance_validation_cp <- function(x) {
  n <- length(x)
  centered <- x - mean(x)
  padding <- 2^ceiling(log2(2 * n))
  transformed <- fft(c(centered, rep(0, padding - n)))
  # Descriptive ACF only: ESS/MCSE are computed separately by posterior.
  Re(fft(transformed * Conj(transformed), inverse = TRUE))[seq_len(n)] /
    padding / rev(seq_len(n))
}

mixing_validation_cp <- function(object, ess) {
  rows <- list()
  for (method in c("Crude PMMH", "Correlated PMMH")) {
    info <- object$sample_info[[method]]
    estimates <- ess[ess$method == method, , drop = FALSE]
    holding <- maximum_holding_run(object$samples[[method]])
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
      ess_direction_mean = estimates$ess_mean[estimates$statistic == "S_shift"]
    )
  }
  do.call(rbind, rows)
}

# Gaussian target whitening supplies a common scale, not the selective covariance itself.
covariance_validation_cp <- function(samples, config) {
  covariances <- lapply(samples, stats::cov)
  correlations <- lapply(samples, stats::cor)
  decomposition <- eigen(config$target$covariance, symmetric = TRUE)
  whitening <- decomposition$vectors %*% diag(1 / sqrt(decomposition$values)) %*%
    t(decomposition$vectors)
  whitened <- lapply(covariances, function(covariance) whitening %*% covariance %*% whitening)
  rows <- list()
  for (method in c("Crude PMMH", "Correlated PMMH")) {
    difference <- covariances[[method]] - covariances$Rejection
    correlation_difference <- correlations[[method]] - correlations$Rejection
    whitened_difference <- whitened[[method]] - whitened$Rejection
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, covariance_frobenius = sqrt(sum(difference^2)),
      relative_frobenius = norm(whitened_difference, "F") / max(norm(whitened$Rejection, "F"), .Machine$double.eps),
      correlation_maximum = max(abs(correlation_difference)),
      whitened_frobenius = norm(whitened_difference, "F")
    )
  }
  list(comparison = do.call(rbind, rows), covariance = covariances,
       correlation = correlations, eigenvalues = lapply(covariances, function(x) {
         eigen(x, symmetric = TRUE, only.values = TRUE)$values
       }))
}

# Inclusive upper-tail indicators use fixed pre-selection Q and S_shift quantiles.
# At IID boundary counts retain Wilson intervals rather than a zero uncertainty claim.
tails_validation_cp <- function(values, config) {
  target <- config$target
  noncentrality <- drop(t(target$gaussian_center) %*% target$precision %*% target$gaussian_center)
  rows <- list()
  for (statistic in c("Q", "S_shift")) for (level in config$diagnostic$tail_levels) {
    threshold <- if (statistic == "Q") {
      stats::qchisq(level, df = length(config$coordinates), ncp = noncentrality)
    } else stats::qnorm(level)
    reference_indicator <- as.numeric(values$Rejection[, statistic] >= threshold)
    reference_probability <- mean(reference_indicator)
    n_reference <- length(reference_indicator)
    reference_se <- sqrt(reference_probability * (1 - reference_probability) / n_reference)
    boundary <- sum(reference_indicator) %in% c(0L, n_reference)
    z <- stats::qnorm(0.975)
    denominator <- 1 + z^2 / n_reference
    center <- (reference_probability + z^2 / (2 * n_reference)) / denominator
    half <- z * sqrt(reference_probability * (1 - reference_probability) / n_reference +
                       z^2 / (4 * n_reference^2)) / denominator
    for (method in names(values)) {
      indicator <- as.numeric(values[[method]][, statistic] >= threshold)
      probability <- mean(indicator)
      estimate <- if (method == "Rejection") {
        list(ess_mean = NA_real_, mcse_mean = if (boundary) NA_real_ else reference_se,
             undefined_reason = if (boundary) "IID reference boundary; Wilson interval retained" else
               "MCMC ESS not applicable to IID reference")
      } else mean_diagnostics(indicator, paste(method, statistic, level, "upper tail"))
      combined <- if (boundary || method == "Rejection") NA_real_ else sqrt(reference_se^2 + estimate$mcse_mean^2)
      standardized <- if (boundary || method == "Rejection") NA_real_ else
        (probability - reference_probability) / combined
      rows[[length(rows) + 1L]] <- data.frame(
        statistic = statistic, level = level, threshold = threshold, method = method,
        n = length(indicator), exceedances = sum(indicator), probability = probability,
        reference_probability = reference_probability, reference_se = reference_se,
        ess_indicator = estimate$ess_mean, mcse_indicator = estimate$mcse_mean, combined_mcse = combined,
        z_difference = standardized, reference_boundary = boundary,
        reference_wilson_lower = max(0, center - half),
        reference_wilson_upper = min(1, center + half),
        undefined_reason = estimate$undefined_reason,
        successes = sum(indicator), failures = sum(1 - indicator),
        transitions = sum(diff(indicator) != 0)
      )
    }
  }
  do.call(rbind, rows)
}

autocorrelation_validation_cp <- function(values) {
  rows <- list()
  for (method in c("Crude PMMH", "Correlated PMMH")) for (statistic in colnames(values[[method]])) {
    covariance <- autocovariance_validation_cp(values[[method]][, statistic])
    acf <- (covariance / covariance[1L])[-1L]
    lags <- seq_len(min(1000L, length(acf)))
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, statistic = statistic, lag = lags, acf = as.numeric(acf[lags])
    )
  }
  do.call(rbind, rows)
}

rolling_validation_cp <- function(values, config) {
  window <- config$diagnostic$rolling_window
  step <- config$diagnostic$rolling_step
  rows <- list()
  for (method in c("Crude PMMH", "Correlated PMMH")) for (statistic in c("T1", "Q", "S_shift")) {
    x <- values[[method]][, statistic]
    starts <- seq.int(1L, length(x) - window + 1L, by = step)
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, statistic = statistic, window = window, start = starts,
      end = starts + window - 1L,
      mean = vapply(starts, function(start) mean(x[seq.int(start, start + window - 1L)]), numeric(1L))
    )
  }
  do.call(rbind, rows)
}

# Analyze saved rejection vectors and production chains; warm-up is already removed.
# All holding states contribute to moments, tail indicators and marginal comparisons.
analyze_validation_cp <- function(object) {
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
  methods <- c("Rejection", "Crude PMMH", "Correlated PMMH")
  values <- lapply(object$samples, statistics_validation_cp, config = config)
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
      estimate <- mean_diagnostics(x, paste(method, statistic, sep = ", "))
      half <- floor(length(x) / 2L)
      ess_rows[[length(ess_rows) + 1L]] <- data.frame(
        method = method, statistic = statistic, n = length(x),
        ess_mean = estimate$ess_mean, mcse_mean = estimate$mcse_mean, tau = estimate$tau,
        undefined_reason = estimate$undefined_reason,
        first_half_mean = mean(x[seq_len(half)]),
        second_half_mean = mean(x[seq.int(half + 1L, length(x))])
      )
    }
  }
  ess <- do.call(rbind, ess_rows)
  for (statistic in c(config$coordinates, "Q", "S_shift")) {
    reference <- values$Rejection[, statistic]
    reference_se <- stats::sd(reference) / sqrt(length(reference))
    for (method in methods[-1L]) {
      x <- values[[method]][, statistic]
      mcse <- ess$mcse_mean[ess$method == method & ess$statistic == statistic]
      # Reference and chain variances are added under their separate-stream assumption.
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
  tails <- tails_validation_cp(values, config)
  covariance <- covariance_validation_cp(object$samples, config)
  diagnostics <- list(ess_mcse = ess, mixing = mixing_validation_cp(object, ess))
  diagnostics$autocorrelation <- autocorrelation_validation_cp(values)
  diagnostics$rolling_means <- rolling_validation_cp(values, config)
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
  object$summary$coordinate_maximum_mcse_mean <- setNames(vapply(methods[-1L], function(method) {
    max(ess$mcse_mean[ess$method == method & ess$statistic %in% config$coordinates])
  }, numeric(1L)), methods[-1L])
  object$config$sampler$ess <- "minimum coordinate posterior::ess_basic(split = TRUE)"
  object
}
