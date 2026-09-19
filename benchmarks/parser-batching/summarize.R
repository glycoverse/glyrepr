output <- "benchmarks/parser-batching/results"
x <- do.call(
  rbind,
  lapply(1:3, function(i) {
    read.csv(file.path(output, paste0("timings-", i, ".csv")))
  })
)
baseline <- x[
  x$variant == "baseline",
  c("worker", "round", "operation", "elapsed")
]
names(baseline)[[4]] <- "baseline_elapsed"
paired <- merge(
  x[x$variant != "baseline", ],
  baseline,
  by = c("worker", "round", "operation")
)
paired$speedup <- paired$baseline_elapsed / paired$elapsed
summary <- do.call(
  rbind,
  lapply(
    split(paired, list(paired$operation, paired$variant), drop = TRUE),
    function(z) {
      data.frame(
        operation = z$operation[[1]],
        variant = z$variant[[1]],
        baseline_median_s = median(z$baseline_elapsed),
        variant_median_s = median(z$elapsed),
        paired_speedup_median = median(z$speedup),
        paired_speedup_min = min(z$speedup),
        paired_speedup_max = max(z$speedup)
      )
    }
  )
)
write.csv(paired, file.path(output, "paired-timings.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 4)
