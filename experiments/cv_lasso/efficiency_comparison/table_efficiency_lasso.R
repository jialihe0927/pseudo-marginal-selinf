#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 1L) stop("Usage: Rscript table_efficiency_lasso.R [efficiency_lasso.rds]")
input_path <- if (length(arguments)) arguments[[1L]] else "efficiency_lasso.rds"
if (!file.exists(input_path)) stop("Results file not found: ", input_path)

result <- readRDS(input_path)
required <- c("method", "final_proposal_scale", "efficiency_per_1000",
              "production_acceptance", "zero_estimate_rate")
if (!all(required %in% names(result$estimates))) stop("Result object has incompatible estimates.")
order <- c("Rejection", "Crude PMMH", "Correlated PMMH")
table <- result$estimates[match(order, result$estimates$method), required, drop = FALSE]
if (anyNA(table$method)) stop("The paper table requires all three methods.")
table$gain <- table$efficiency_per_1000 / table$efficiency_per_1000[1L]
table <- table[, c("method", "final_proposal_scale", "efficiency_per_1000", "gain",
                   "production_acceptance", "zero_estimate_rate")]
names(table) <- c("Method", "h_final", "Efficiency", "Gain", "Acc.", "Zero rate")
print(table, row.names = FALSE, digits = 15)

