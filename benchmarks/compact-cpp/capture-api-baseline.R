pkgload::load_all(quiet = TRUE)
cat("BASELINE LOADED\n")
dir <- file.path("benchmarks", "compact-cpp", "results-api")
exports <- getNamespaceExports("glyrepr")
api <- lapply(sort(exports), function(n) {
  value <- getExportedValue("glyrepr", n)
  if (is.function(value)) formals(value) else class(value)
})
names(api) <- sort(exports)
saveRDS(api, file.path(dir, "api-baseline.rds"))
x <- readRDS("benchmarks/compact-cpp/results/corpus.rds")
graph_signature <- function(g) {
  list(
    directed = igraph::is_directed(g),
    edges = unname(igraph::as_edgelist(g, names = FALSE)),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
# A graph-free signature snapshot avoids external pointer identity comparisons.
for (start in seq(1L, length(x), by = 500L)) {
  ids <- start:min(start + 499L, length(x))
  out <- as_glycan_structure(x[ids])
  saveRDS(
    list(
      input = x[ids],
      values = as.character(out),
      names = names(out),
      graphs = lapply(attr(out, "graphs"), graph_signature)
    ),
    file.path(dir, paste0("baseline-", start, ".rds"))
  )
  cat("Baseline", max(ids), "\n")
}
