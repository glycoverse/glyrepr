# Run from the repository root: Rscript benchmarks/rcpp-prototype/run.R verify
# Or: Rscript benchmarks/rcpp-prototype/run.R bench 1 (repeat with 2 and 3).
pkgload::load_all(quiet = TRUE)
source("benchmarks/rcpp-prototype/prototype.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "verify"
output <- "benchmarks/rcpp-prototype/results"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("glyrepr")
all_graphs <- attr(glydb::glydb_structures(), "graphs")
is_ordinary <- function(g) {
  !any(
    c("floating_parts", "floating_substituents") %in%
      igraph::graph_attr_names(g)
  )
}
ordinary <- Filter(is_ordinary, all_graphs)
set.seed(123)
graphs <- unname(ordinary[sample.int(length(ordinary), 500L)])
strings <- vapply(graphs, graph_to_iupac, character(1))
permuted <- lapply(graphs, function(g) {
  g <- igraph::permute(g, sample.int(igraph::vcount(g)))
  ns$.permute_edges(g, sample.int(igraph::ecount(g)))
})
graph_signature <- function(g) {
  list(
    edges = igraph::as_edgelist(g, names = FALSE),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
check_graph <- function(g) {
  if (is_ordinary(g)) {
    cache <- ns$build_seq_cache(g)
    old <- ns$seq_glycan_order_iupac(cache$root, cache)
    new <- prototype_env$sequence_graph(g)
    stopifnot(
      identical(as.numeric(old$vertices), new$vertices),
      identical(as.numeric(old$edges), new$edges),
      identical(old$iupac, new$iupac)
    )
  }
  old <- prototype_env$baseline_canonical(g)
  new <- prototype_env$canonical_cpp(g)
  stopifnot(
    identical(old$iupac, new$iupac),
    identical(graph_signature(old$graph), graph_signature(new$graph)),
    identical(
      prototype_env$baseline_iupac(old$graph),
      prototype_env$iupac_cpp(old$graph)
    )
  )
}
if (mode == "verify") {
  for (i in seq_along(all_graphs)) {
    check_graph(all_graphs[[i]])
    if (i %% 1000L == 0L) cat("Corpus parity:", i, "\n")
  }
  for (g in permuted) {
    check_graph(g)
  }
  fixtures <- as_glycan_structure(c(
    "Gal(b1-4)GlcNAc-ol",
    "Neu5Ac(a2-3/6)Gal",
    "Gal(?1-?)GlcNAc",
    "Gal6S(b1-4)GlcNAc",
    "Gal(?1-?)[Gal(?1-?)]GlcNAc"
  ))
  fixture_graphs <- as.list(fixtures)
  fixture_graphs[[1]] <- igraph::set_graph_attr(
    fixture_graphs[[1]],
    "custom",
    list(a = 1)
  )
  fixture_graphs[[1]] <- igraph::set_vertex_attr(
    fixture_graphs[[1]],
    "custom",
    value = seq_len(igraph::vcount(fixture_graphs[[1]]))
  )
  fixture_graphs[[1]] <- igraph::set_edge_attr(
    fixture_graphs[[1]],
    "custom",
    value = seq_len(igraph::ecount(fixture_graphs[[1]]))
  )
  locales <- unique(c(Sys.getlocale("LC_COLLATE"), "C", "en_US.UTF-8"))
  original_locale <- Sys.getlocale("LC_COLLATE")
  checked_locales <- character()
  for (locale in locales) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) {
      for (g in c(fixture_graphs, permuted[1:50])) {
        check_graph(g)
      }
      checked_locales <- c(checked_locales, locale)
    }
  }
  Sys.setlocale("LC_COLLATE", original_locale)
  named <- c(
    first = strings[[1]],
    missing = NA_character_,
    again = strings[[1]]
  )
  old <- as_glycan_structure(named)
  new <- with_prototype(as_glycan_structure(named))
  stopifnot(
    identical(names(old), names(new)),
    identical(as.character(old), as.character(new))
  )
  cat("Running complete package test suite with prototype bindings\n")
  Sys.setenv(NOT_CRAN = "true")
  tests <- with_prototype(testthat::test_dir(
    "tests/testthat",
    env = new.env(parent = ns),
    reporter = "summary",
    package = "glyrepr",
    load_package = "none",
    stop_on_failure = TRUE
  ))
  results <- as.data.frame(tests)
  write.csv(
    results[c(
      "file",
      "test",
      "nb",
      "failed",
      "skipped",
      "error",
      "warning",
      "passed"
    )],
    file.path(output, "tests.csv"),
    row.names = FALSE
  )
  writeLines(
    c(
      paste("Corpus graphs:", length(all_graphs)),
      paste("Ordinary kernel parity:", length(ordinary)),
      paste("Floating fallback parity:", length(all_graphs) - length(ordinary)),
      paste("Permuted graphs:", length(permuted)),
      paste("Locales:", paste(checked_locales, collapse = ", ")),
      paste("Passed expectations:", sum(results$passed)),
      paste("Failed:", sum(results$failed)),
      paste("Errors:", sum(results$error)),
      paste("Warnings:", sum(results$warning)),
      paste("Skipped:", sum(results$skipped))
    ),
    file.path(output, "validation.txt")
  )
} else if (mode == "bench") {
  worker <- as.integer(args[[2]])
  operations <- list(
    sequence_kernel = function() {
      lapply(graphs, function(g) {
        cache <- ns$build_seq_cache(g)
        ns$seq_glycan_order_iupac(cache$root, cache)
      })
    },
    graph_to_iupac = function() lapply(graphs, ns$graph_to_iupac),
    canonical_graph_constructor = function() as_glycan_structure(graphs),
    permuted_graph_constructor = function() as_glycan_structure(permuted),
    iupac_constructor = function() as_glycan_structure(strings)
  )
  run <- function(name, implementation) {
    if (implementation == "R") {
      return(operations[[name]]())
    }
    if (name == "sequence_kernel") {
      return(lapply(graphs, prototype_env$sequence_graph))
    }
    with_prototype(operations[[name]]())
  }
  # Warm both paths; compilation and corpus preparation are outside timers.
  for (name in names(operations)) {
    for (implementation in c("R", "Rcpp")) {
      invisible(run(name, implementation))
    }
  }
  rows <- list()
  for (round in 1:3) {
    order <- if ((worker + round) %% 2L) c("R", "Rcpp") else c("Rcpp", "R")
    for (name in names(operations)) {
      for (implementation in order) {
        gc()
        elapsed <- system.time(invisible(run(name, implementation)))[[
          "elapsed"
        ]]
        rows[[length(rows) + 1L]] <- data.frame(
          worker,
          round,
          operation = name,
          implementation,
          elapsed
        )
        cat(worker, round, name, implementation, elapsed, "\n")
      }
    }
  }
  write.csv(
    do.call(rbind, rows),
    file.path(output, paste0("timings-", worker, ".csv")),
    row.names = FALSE
  )
  writeLines(strings, file.path(output, "sample-iupac.txt"))
  capture.output(
    sessionInfo(),
    file = file.path(output, paste0("session-", worker, ".txt"))
  )
} else {
  stop("Use verify or bench WORKER_ID")
}
