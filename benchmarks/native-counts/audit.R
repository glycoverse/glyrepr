library(glyrepr)
ns <- asNamespace("glyrepr")
reference_graph <- get(".graph_to_composition_reference", ns)
new_comp <- get("new_glycan_composition", ns)
reference <- function(x) new_comp(smap(x, reference_graph))
timing <- function(f) {
  elapsed <- system.time(f())[["elapsed"]]
  n <- max(1L, min(100L, ceiling(0.1 / max(elapsed, 0.001))))
  samples <- replicate(
    5,
    system.time(
      for (i in seq_len(n)) {
        f()
      }
    )[["elapsed"]] /
      n
  )
  list(median = median(samples), samples = samples, iterations = n)
}
inputs <- list(
  repeated = rep(c(n_glycan_core(), o_glycan_core_1(), o_glycan_core_2()), 100),
  unique = as_glycan_structure(vapply(
    1:50,
    function(n) {
      paste0(paste(rep("Gal(b1-4)", n), collapse = ""), "Glc(a1-")
    },
    character(1)
  ))
)
rows <- list()
raw <- list()
for (workload in names(inputs)) {
  x <- inputs[[workload]]
  for (query in list("composition", NULL, "Hex", "Gal", "S")) {
    operation <- if (is.null(query)) "total" else query
    old <- if (identical(query, "composition")) {
      function() reference(x)
    } else {
      function() count_mono(reference(x), query)
    }
    current <- if (identical(query, "composition")) {
      function() as_glycan_composition(x)
    } else {
      function() count_mono(x, query)
    }
    stopifnot(identical(old(), current()))
    before <- timing(old)
    after <- timing(current)
    key <- paste(workload, operation)
    raw[[key]] <- list(reference = before, native = after)
    rows[[key]] <- data.frame(
      workload,
      operation,
      reference_seconds = before$median,
      native_seconds = after$median,
      speedup = before$median / after$median
    )
  }
}
write.csv(
  do.call(rbind, rows),
  "benchmarks/native-counts/results.csv",
  row.names = FALSE
)
saveRDS(
  list(inputs = inputs, timings = raw, session = sessionInfo()),
  "benchmarks/native-counts/audit.rds"
)
print(do.call(rbind, rows))
