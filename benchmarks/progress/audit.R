# Run from the repository root with a separately installed package library.
# Rscript benchmarks/progress/audit.R /path/to/library baseline /tmp/progress-audit
# Rscript benchmarks/progress/audit.R /path/to/new-library current /tmp/progress-audit
args <- commandArgs(TRUE)
library(glyrepr, lib.loc = args[[1]])
mode <- args[[2]]
root <- args[[3]]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
corpus <- readRDS('benchmarks/native-graph/results/corpus.rds')
strings <- rep(c(corpus, NA_character_), 4L)
names(strings) <- paste0('s', seq_along(strings))
graphs <- as.list(as_glycan_structure(corpus))
records <- lapply(graphs, getFromNamespace('.compact_graph_record', 'glyrepr'))
records <- rep(c(records, list(NULL)), 4L)
names(records) <- names(strings)
inputs <- list(strings = strings, arrays = records, graphs = graphs)
signature <- function(x) {
  list(
    values = as.character(x),
    graphs = lapply(attr(x, 'graphs'), function(g) {
      list(
        vertices = igraph::vertex_attr(g),
        edges = igraph::as_edgelist(g),
        edge_attrs = igraph::edge_attr(g),
        attrs = igraph::graph_attr(g)
      )
    })
  )
}
options(cli.progress_show_after = Inf)
rows <- list()
outputs <- list()
for (kind in names(inputs)) {
  fun <- if (kind == 'arrays') structure_from_arrays else as_glycan_structure
  for (policy in c('error', 'na')) {
    for (p in if (mode == 'baseline') 'off' else c('off', 'on')) {
      call <- list(x = inputs[[kind]], on_failure = policy)
      if (mode != 'baseline') {
        call$progress <- p == 'on'
      }
      key <- paste(kind, policy, p, sep = '-')
      outputs[[key]] <- signature(do.call(fun, call))
      times <- replicate(5, {
        gc()
        system.time(do.call(fun, call))[[3L]]
      })
      rows[[key]] <- data.frame(
        kind,
        policy,
        progress = p,
        n = length(inputs[[kind]]),
        seconds = median(times)
      )
    }
  }
}
saveRDS(outputs, file.path(root, paste0(mode, '-outputs.rds')))
write.csv(
  do.call(rbind, rows),
  file.path(root, paste0(mode, '-timings.csv')),
  row.names = FALSE
)

if (mode == "current") {
  baseline <- readRDS(file.path(root, "baseline-outputs.rds"))
  for (key in names(outputs)) {
    stopifnot(identical(outputs[[key]], baseline[[sub("-on$", "-off", key)]]))
  }
  message("All 12 output signatures match the frozen baseline.")
}
