main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  source("../../../diagnostics.R", local = TRUE)
  if (length(arguments) > 1L) stop("Usage: Rscript table_efficiency_lasso.R [efficiency_lasso.rds]")
  input <- if (length(arguments)) arguments[[1L]] else "efficiency_lasso.rds"

  table <- efficiency_table(readRDS(input))
  print(table[setdiff(names(table), "method_id")], row.names = FALSE, digits = 15)
}

if (sys.nframe() == 0L) {
  main()
}

