# Rscript benchmarks/array-backend/benchmark.R /path/to/release/library
args <- commandArgs(TRUE)
library(glyrepr, lib.loc = args[[1]])
ns <- asNamespace("glyrepr")
build <- get(".compact_build_graph", ns)
validate <- get(".validate_structure_arrays", ns)
native <- get(".compact_arrays_native", ns)
bliss <- get(".compact_bliss_labels", ns)
outdir <- "benchmarks/array-backend/results"
records <- readRDS(file.path(outdir, "workload.rds"))
workloads <- list(
  distinct_500 = records,
  repeated_1000 = rep(records[1:20], 50)
)
rows <- list()
for (label in names(workloads)) {
  a <- workloads[[label]]
  graphs <- lapply(a, build)
  checked <- lapply(a, validate)
  methods <- list(
    reference_constructor = function() as_glycan_structure(lapply(a, build)),
    array_constructor = function() structure_from_arrays(a),
    graph_split = function() {
      g <- lapply(graphs, function(g) {
        canonicalize_glycan_graph(validate_glycan_graph(g))
      })
      list(graphs = g, iupac = vapply(g, graph_to_iupac, character(1)))
    },
    graph_fused = function() canonicalize_glycan_graphs(graphs),
    native_arrays_only = function() {
      native(
        checked,
        base::order,
        bliss,
        identical(Sys.getlocale("LC_COLLATE"), "C")
      )
    }
  )
  cat(label, "parity and warmup", format(Sys.time()), "\n")
  flush.console()
  reference <- methods$reference_constructor()
  stopifnot(identical(
    as.character(reference),
    as.character(methods$array_constructor())
  ))
  stopifnot(identical(methods$graph_split()$iupac, methods$graph_fused()$iupac))
  for (f in methods) {
    invisible(f())
  }
  for (round in 1:3) {
    order <- if (round %% 2) names(methods) else rev(names(methods))
    for (method in order) {
      gc()
      cat(label, "round", round, method, "\n")
      flush.console()
      elapsed <- system.time(invisible(methods[[method]]()))[["elapsed"]]
      rows[[length(rows) + 1L]] <- data.frame(
        workload = label,
        method,
        round,
        elapsed
      )
      write.csv(
        do.call(rbind, rows),
        file.path(outdir, "timings.csv"),
        row.names = FALSE
      )
      cat("seconds", elapsed, "\n")
      flush.console()
    }
  }
  cat(label, "complete\n")
}
write.csv(
  do.call(rbind, rows),
  file.path(outdir, "timings.csv"),
  row.names = FALSE
)
writeLines(
  c(capture.output(sessionInfo()), paste("library", args[[1]])),
  file.path(outdir, "benchmark-session.txt")
)
