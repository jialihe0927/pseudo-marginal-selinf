#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source("samplers_fd.R", local = TRUE)
source("analysis_fd.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args)) args[[1L]] else "rerun_fd.rds"
if (length(args) > 1L) stop("Usage: Rscript run_fd.R [output.rds]")
if (file.exists(output_path) || dir.exists(output_path)) {
  stop("Output already exists: ", output_path)
}

b_values <- 2:8
fixed_cpm_caps <- c(
  `2` = 3000L,
  `3` = 3135L,
  `4` = 3360L,
  `5` = 4632L,
  `6` = 8084L,
  `7` = 10111L,
  `8` = 11944L
)

config <- list(
  b = b_values,
  regimes = c("fixed_a", "adaptive_a"),
  fixed_a = 0.9,
  adaptive_a = "1 - b^-2",
  lambda = 2,
  epsilon = 0.10,
  n_chains = 1000L,
  checkpoint_log_step = log(1.05),
  stable_checkpoints = 8L,
  warmup = 0L,
  observable = "I(X <= selected-distribution median)",
  variance_ratio = "4 * across-chain variance of chain averages",
  point_estimate = "earliest start of eight consecutive passing checkpoints",
  selector_cost = "one selector evaluation per transition",
  seed_scheme = "L'Ecuyer-CMRG; seed 9100000 + configuration index",
  simulation_plan = data.frame(
    regime = c(rep("fixed_a", 14L), rep("adaptive_a", 7L)),
    b = rep(b_values, 3L),
    method = c(
      rep("Independent PM", 7L),
      rep("CPM", 7L),
      rep("Independent PM", 7L)
    ),
    maximum_observations = c(
      rep(2000000L, 7L),
      unname(fixed_cpm_caps[as.character(b_values)]),
      rep(500000L, 7L)
    ),
    post_certification_checkpoints = c(
      rep(0L, 7L),
      rep(8L, 7L),
      rep(8L, 7L)
    ),
    stringsAsFactors = FALSE
  )
)

RNGkind("L'Ecuyer-CMRG")
simulation_rows <- vector("list", nrow(config$simulation_plan))

for (configuration_id in seq_len(nrow(config$simulation_plan))) {
  plan <- config$simulation_plan[configuration_id, ]
  a <- if (plan$regime == "fixed_a") config$fixed_a else 1 - plan$b^-2
  rho <- if (plan$method == "CPM") 1 - plan$b^-2 else 0
  proposal_scale <- config$lambda / (plan$b * sqrt(a))

  message(
    sprintf(
      "Running %s, %s, b=%d",
      plan$regime,
      plan$method,
      plan$b
    )
  )
  estimate <- estimate_fixed_accuracy(
    method = plan$method,
    a = a,
    b = plan$b,
    rho = rho,
    proposal_scale = proposal_scale,
    maximum_observations = plan$maximum_observations,
    n_chains = config$n_chains,
    epsilon = config$epsilon,
    stable_checkpoints = config$stable_checkpoints,
    post_certification_checkpoints = plan$post_certification_checkpoints,
    seed = 9100000L + configuration_id,
    checkpoint_log_step = config$checkpoint_log_step
  )
  estimate$regime <- plan$regime
  simulation_rows[[configuration_id]] <- estimate
}

simulation_estimates <- do.call(rbind, simulation_rows)
rejection_time <- vapply(
  b_values,
  rejection_fixed_accuracy_time,
  numeric(1L),
  epsilon = config$epsilon
)
rejection_estimates <- do.call(rbind, lapply(config$regimes, function(regime) {
  data.frame(
    b = b_values,
    method = "Rejection",
    a = if (regime == "fixed_a") config$fixed_a else 1 - b_values^-2,
    rho = NA_real_,
    proposal_scale = NA_real_,
    complexity_n = rejection_time,
    observations_simulated = NA_real_,
    transitions_per_chain = NA_real_,
    selector_evaluations = NA_real_,
    resolved = TRUE,
    estimate_type = "exact analytic crossing",
    result_type = "analytic",
    regime = regime,
    stringsAsFactors = FALSE
  )
}))

estimates <- rbind(simulation_estimates, rejection_estimates)
method_order <- c("Rejection", "Independent PM", "CPM")
estimates <- estimates[
  order(
    match(estimates$regime, config$regimes),
    estimates$b,
    match(estimates$method, method_order)
  ),
]
row.names(estimates) <- NULL

result <- list(
  config = config,
  estimates = estimates
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_path, version = 3)
message("Saved results to ", normalizePath(output_path, winslash = "/"))
