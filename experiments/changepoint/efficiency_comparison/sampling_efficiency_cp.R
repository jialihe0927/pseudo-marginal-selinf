#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (!requireNamespace("posterior", quietly = TRUE)) stop("Package 'posterior' is required.")
source("functions_cp.R", local = TRUE)
source("samplers_cp.R", local = TRUE)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[[1L]] %in% c("pmmh", "rejection", "all")) {
  stop("Usage: Rscript sampling_efficiency_cp.R {pmmh|rejection|all} output.rds")
}
mode <- arguments[[1L]]
output_path <- arguments[[2L]]
if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

paper_results <- readRDS("efficiency_cp.rds")
benchmark <- paper_results$benchmark
config <- paper_results$config
validate_efficiency_config_cp(config)

summarize_pmmh <- function(result) {
  coordinate_names <- paste0("T", seq_len(ncol(result$draws)))
  coordinate_ess <- vapply(
    seq_len(ncol(result$draws)),
    function(index) {
      bulk_ess_cp(
        result$draws[, index],
        method = result$method,
        coordinate = coordinate_names[index]
      )
    },
    numeric(1L)
  )
  minimum_ess <- if (any(!is.finite(coordinate_ess))) {
    NA_real_
  } else {
    min(coordinate_ess)
  }

  list(
    estimate = data.frame(
      method = result$method,
      selector_evaluations = result$selector_evaluations,
      effective_samples = minimum_ess,
      efficiency_per_1000 = 1000 * minimum_ess / result$selector_evaluations,
      production_acceptance = mean(result$accepted),
      zero_estimate_rate = mean(!result$proposed_selected),
      final_proposal_scale = result$final_proposal_scale,
      stringsAsFactors = FALSE
    ),
    coordinate_ess = data.frame(
      method = result$method,
      coordinate = coordinate_names,
      bulk_ess = coordinate_ess,
      stringsAsFactors = FALSE
    ),
    production_draws = result$draws
  )
}

parts <- list()
if (mode %in% c("pmmh", "all")) {
  methods <- c("Crude PMMH", "Correlated PMMH")
  runs <- lapply(
    methods,
    function(method) run_pmmh_cp(benchmark, config, method)
  )
  summaries <- lapply(runs, summarize_pmmh)
  parts$pmmh <- list(
    estimates = do.call(
      rbind,
      lapply(summaries, `[[`, "estimate")
    ),
    coordinate_ess = do.call(
      rbind,
      lapply(summaries, `[[`, "coordinate_ess")
    ),
    production_draws = setNames(
      lapply(summaries, `[[`, "production_draws"),
      methods
    )
  )
}
if (mode %in% c("rejection", "all")) {
  rejection <- run_rejection_cp(benchmark, config)
  parts$rejection <- list(
    estimates = data.frame(
      method = "Rejection",
      selector_evaluations = rejection$attempts,
      effective_samples = rejection$hits,
      efficiency_per_1000 = 1000 * rejection$hits / rejection$attempts,
      production_acceptance = NA_real_,
      zero_estimate_rate = NA_real_,
      final_proposal_scale = NA_real_,
      stringsAsFactors = FALSE
    ),
    accepted_draws = rejection$accepted
  )
}

result <- list(
  config = config,
  benchmark = benchmark,
  estimates = do.call(rbind, lapply(parts, `[[`, "estimates")),
  coordinate_ess = if (!is.null(parts$pmmh)) parts$pmmh$coordinate_ess else NULL,
  simulation_output = parts
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_path, version = 3)

