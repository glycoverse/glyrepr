pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
x <- readRDS(file.path(compact_dir, "results", "corpus.rds"))
# All formerly unsupported inputs, plus supported trees, repeats, and NA.
coverage <- read.csv(file.path(compact_dir, "results", "coverage.csv"))
x <- x[c(
  which(coverage$status == "fallback"),
  head(which(coverage$status == "ok"), 50)
)]
counts <- new.env(parent = emptyenv())
functions <- c(
  "make_graph",
  "make_empty_graph",
  "add_edges",
  "permute",
  "graph_from_data_frame",
  "canonical_permutation"
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
.GlobalEnv$final_count <- .GlobalEnv$bliss_count <- 0L
trace(
  "compact_graph",
  print = FALSE,
  tracer = quote(.GlobalEnv$final_count <- .GlobalEnv$final_count + 1L)
)
trace(
  "compact_bliss",
  print = FALSE,
  tracer = quote(.GlobalEnv$bliss_count <- .GlobalEnv$bliss_count + 1L)
)
out <- compact_structure(c(x, x, NA_character_), fallback = FALSE)
for (fn in functions) {
  untrace(fn, where = asNamespace("igraph"))
}
untrace("compact_graph")
untrace("compact_bliss")
observed <- vapply(functions, function(fn) counts[[fn]], 0L)
stopifnot(
  final_count == length(attr(out, "graphs")),
  observed[["make_graph"]] == final_count + bliss_count,
  observed[["canonical_permutation"]] == bliss_count,
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
observed <- c(
  observed,
  final_glycan_graph = final_count,
  augmented_bliss_graph = bliss_count
)
write.csv(
  data.frame(function_name = names(observed), count = unname(observed)),
  file.path(compact_dir, "results-v2", "graph-calls.csv"),
  row.names = FALSE
)
