args <- commandArgs(TRUE)
.libPaths(c(args[[1]], .libPaths()))
library(glyrepr)
worker <- as.integer(args[[2]])
method <- args[[3]]
x <- readRDS("benchmarks/compact-cpp/results/workloads.rds")$mixed_500
inputs <- list(
  valid_error = x,
  valid_na = x,
  partial_na = c(x, rep("bad", 20), NA_character_)
)
rows <- list()
for (name in names(inputs)) {
  input <- inputs[[name]]
  policy <- if (name == "valid_error") "error" else "na"
  invisible(as_glycan_structure(input[1:3], on_failure = policy))
  for (round in 1:2) {
    gc()
    elapsed <- system.time(
      value <- suppressWarnings(as_glycan_structure(input, on_failure = policy))
    )[["elapsed"]]
    stopifnot(length(value) == length(input))
    rows[[length(rows) + 1L]] <- data.frame(
      worker,
      method,
      workload = name,
      round,
      elapsed
    )
  }
}
write.csv(
  do.call(rbind, rows),
  file.path(
    "benchmarks",
    "compact-cpp",
    "results-api",
    paste0("timings-", method, "-", worker, ".csv")
  ),
  row.names = FALSE
)
