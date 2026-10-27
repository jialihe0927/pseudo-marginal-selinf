options(stringsAsFactors = FALSE)

source("functions_nonnull_lasso.R", local = TRUE)
source("samplers_nonnull_lasso.R", local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript inference_nonnull_lasso.R [output.rds] [bank.rds]")
output_path <- if (length(args)) args[[1L]] else "rerun_inference_nonnull_lasso.rds"
if (file.exists(output_path) || dir.exists(output_path)) stop("Output already exists: ", output_path)

bank_path <- if (length(args) >= 2L) args[[2L]] else "inference_nonnull_lasso.rds"
if (!file.exists(bank_path)) stop("Data bank not found: ", bank_path)
input <- readRDS(bank_path)
bank <- input[c("config", "datasets", "rng_assignments")]
if (!is.null(input$raw_attempts)) bank$raw_attempts <- input$raw_attempts
validate_bank_nonnull_lasso(bank)
config <- bank$config
datasets <- bank$datasets
assignments <- bank$rng_assignments
dataset_keys <- vapply(datasets, function(dataset) {
  paste(dataset$delta, dataset$dataset_id, sep = ":")
}, character(1L))

rows <- vector("list", nrow(assignments) + length(datasets))
for (job in seq_len(nrow(assignments))) {
  assignment <- assignments[job, , drop = FALSE]
  index <- match(paste(assignment$delta, assignment$dataset_id, sep = ":"), dataset_keys)
  message("delta=", assignment$delta, ", dataset=", assignment$dataset_id,
          "/", config$datasets_per_delta, ", method=", assignment$method_id)
  rows[[job]] <- run_pmmh_nonnull_lasso(
    datasets[[index]], assignment$method_id, assignment$rng_state[[1L]], config
  )
}
for (index in seq_along(datasets)) {
  dataset <- datasets[[index]]
  p_value <- fission_pvalue_nonnull_lasso(dataset, config)
  rows[[nrow(assignments) + index]] <- result_row_nonnull_lasso(
    dataset, "data_fission", p_value, config
  )
}
data <- do.call(rbind, rows)
data <- data[order(match(data$delta, config$delta_grid), data$dataset_id,
                   match(data$method_id, config$methods)), , drop = FALSE]
row.names(data) <- NULL
result <- c(bank, list(dataset_results = data, summary = summarize_nonnull_lasso(data, config)))
validate_result_nonnull_lasso(result)
saveRDS(result, output_path, version = 3)
message("Saved results to ", normalizePath(output_path, winslash = "/"))
