output <- "benchmarks/rcpp-prototype/results"
x <- do.call(
  rbind,
  lapply(1:3, function(i) {
    read.csv(file.path(output, paste0("timings-", i, ".csv")))
  })
)
r <- x[x$implementation == "R", ]
c <- x[x$implementation == "Rcpp", ]
paired <- merge(
  r,
  c,
  by = c("worker", "round", "operation"),
  suffixes = c("_R", "_Rcpp")
)
paired$speedup <- paired$elapsed_R / paired$elapsed_Rcpp
summary <- do.call(
  rbind,
  lapply(unique(paired$operation), function(operation) {
    z <- paired[paired$operation == operation, ]
    data.frame(
      operation,
      R_median_s = median(z$elapsed_R),
      Rcpp_median_s = median(z$elapsed_Rcpp),
      paired_speedup_median = median(z$speedup),
      paired_speedup_min = min(z$speedup),
      paired_speedup_max = max(z$speedup)
    )
  })
)
write.csv(paired, file.path(output, "paired-timings.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 4)
