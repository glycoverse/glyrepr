.libPaths(c(commandArgs(TRUE)[[1]], .libPaths()))
library(glyrepr)
dir <- file.path("benchmarks", "compact-cpp", "results-api")
graph_signature <- function(g) {
  list(
    directed = igraph::is_directed(g),
    edges = unname(igraph::as_edgelist(g, names = FALSE)),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
.GlobalEnv$reference_calls <- 0L
trace(
  ".parse_iupac_condensed_single",
  where = asNamespace("glyrepr"),
  print = FALSE,
  tracer = quote(.GlobalEnv$reference_calls <- .GlobalEnv$reference_calls + 1L)
)
count <- 0L
for (path in list.files(dir, "^baseline-[0-9]+\\.rds$", full.names = TRUE)) {
  baseline <- readRDS(path)
  for (policy in c("error", "na")) {
    out <- as_glycan_structure(baseline$input, on_failure = policy)
    stopifnot(
      identical(as.character(out), baseline$values),
      identical(names(out), baseline$names),
      identical(lapply(attr(out, "graphs"), graph_signature), baseline$graphs)
    )
  }
  count <- count + length(baseline$input)
}
untrace(".parse_iupac_condensed_single", where = asNamespace("glyrepr"))
stopifnot(count == 8573L, reference_calls == 0L)
writeLines(
  c(
    paste("corpus", count),
    "policies error, na",
    paste("reference_parser_calls", reference_calls),
    "PASS: installed package exact graph/key parity under both failure policies"
  ),
  file.path(dir, "installed-parity.txt")
)
