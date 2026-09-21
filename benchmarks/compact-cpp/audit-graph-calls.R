pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
x <- readRDS(file.path(compact_dir, "results", "workloads.rds"))$fast_500
counts <- new.env(parent = emptyenv())
functions <- c(
  "make_graph",
  "make_empty_graph",
  "add_edges",
  "permute",
  "graph_from_data_frame"
)
for (fn in functions) {
  counts[[fn]] <- 0L
  trace(
    fn,
    where = asNamespace("igraph"),
    print = FALSE,
    tracer = substitute(
      .GlobalEnv$counts[[NAME]] <- .GlobalEnv$counts[[NAME]] + 1L,
      list(NAME = fn)
    )
  )
}
candidate <- compact_structure(c(x, x, NA_character_))
for (fn in functions) {
  untrace(fn, where = asNamespace("igraph"))
}
observed <- vapply(functions, function(fn) counts[[fn]], 0L)
stopifnot(
  observed[["make_graph"]] == length(attr(candidate, "graphs")),
  all(
    observed[c(
      "make_empty_graph",
      "add_edges",
      "permute",
      "graph_from_data_frame"
    )] ==
      0L
  )
)
write.csv(
  data.frame(function_name = names(observed), count = unname(observed)),
  file.path(compact_dir, "results", "graph-calls.csv"),
  row.names = FALSE
)
cat("Graph creation audit passed\n")
