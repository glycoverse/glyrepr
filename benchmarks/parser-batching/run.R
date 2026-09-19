# Run from repository root; see README.md for commands.
pkgload::load_all(quiet = TRUE)
source("benchmarks/parser-batching/prototype.R")
source("benchmarks/rcpp-prototype/prototype.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
output <- "benchmarks/parser-batching/results"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
ns <- asNamespace("glyrepr")
strings <- readLines("benchmarks/rcpp-prototype/results/sample-iupac.txt")
stopifnot(length(strings) == 500L)
run_variant <- function(variant, code) {
  switch(
    variant,
    baseline = force(code),
    batch = with_batch_parser(force(code)),
    sequence = with_prototype(force(code)),
    combined = with_batch_parser(with_prototype(force(code)))
  )
}
graph_signature <- function(g) {
  list(
    edges = igraph::as_edgelist(g, names = FALSE),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
if (mode == "verify") {
  corpus <- names(attr(glydb::glydb_structures(), "graphs"))
  for (i in seq_along(corpus)) {
    old <- ns$.parse_iupac_condensed_single(corpus[[i]])
    new <- with_batch_parser(ns$.parse_iupac_condensed_single(corpus[[i]]))
    stopifnot(identical(graph_signature(old), graph_signature(new)))
    if (i %% 1000L == 0L) cat("Raw parser parity:", i, "\n")
  }
  fixture <- c(
    strings,
    missing = NA_character_,
    duplicate = strings[[1]],
    alditol = "Gal(b1-4)GlcNAc-ol",
    ambiguous = "Neu5Ac(a2-3/6)Gal",
    substituent = "Gal6S(b1-4)GlcNAc",
    floating = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  expected <- as_glycan_structure(fixture)
  for (variant in c("batch", "combined")) {
    actual <- run_variant(variant, as_glycan_structure(fixture))
    stopifnot(
      identical(as.character(expected), as.character(actual)),
      identical(names(expected), names(actual)),
      identical(
        lapply(attr(expected, "graphs"), graph_signature),
        lapply(attr(actual, "graphs"), graph_signature)
      )
    )
  }
  # Include every corpus entry in public constructor parity, including floating metadata.
  expected <- as_glycan_structure(corpus)
  for (variant in c("batch", "combined")) {
    actual <- run_variant(variant, as_glycan_structure(corpus))
    stopifnot(
      identical(as.character(expected), as.character(actual)),
      identical(
        lapply(attr(expected, "graphs"), graph_signature),
        lapply(attr(actual, "graphs"), graph_signature)
      )
    )
    cat("Corpus constructor parity:", variant, length(corpus), "\n")
  }
  Sys.setenv(NOT_CRAN = "true")
  counts <- list()
  for (variant in c("batch", "combined")) {
    cat("Package tests:", variant, "\n")
    result <- run_variant(
      variant,
      testthat::test_dir(
        "tests/testthat",
        env = new.env(parent = ns),
        package = "glyrepr",
        load_package = "none",
        reporter = "summary",
        stop_on_failure = TRUE
      )
    )
    result <- as.data.frame(result)
    write.csv(
      result[c(
        "file",
        "test",
        "nb",
        "failed",
        "skipped",
        "error",
        "warning",
        "passed"
      )],
      file.path(output, paste0("tests-", variant, ".csv")),
      row.names = FALSE
    )
    counts[[variant]] <- data.frame(
      variant,
      passed = sum(result$passed),
      failed = sum(result$failed),
      errors = sum(result$error),
      warnings = sum(result$warning),
      skipped = sum(result$skipped)
    )
  }
  write.csv(
    do.call(rbind, counts),
    file.path(output, "validation.csv"),
    row.names = FALSE
  )
  writeLines(
    paste(
      "Exact raw-parser and public-constructor corpus parity:",
      length(corpus)
    ),
    file.path(output, "parity.txt")
  )
} else if (mode == "bench") {
  worker <- as.integer(args[[2]])
  # Prepared data excludes all string processing from graph-assembly timing.
  parsed <- lapply(strings, ns$.parse_iupac_tree_single)
  arrays <- lapply(parsed, function(g) {
    list(
      mono = igraph::vertex_attr(g, "mono"),
      sub = igraph::vertex_attr(g, "sub"),
      edges = as.integer(t(igraph::as_edgelist(g, names = FALSE))),
      linkage = igraph::edge_attr(g, "linkage")
    )
  })
  for (a in arrays) {
    stopifnot(identical(
      graph_signature(batch_env$graph_incremental(a)),
      graph_signature(batch_env$graph_from_arrays(a))
    ))
  }
  workloads <- list(
    graph_assembly = c("baseline", "batch"),
    complete_tree_parser = c("baseline", "batch"),
    public_constructor = c("baseline", "batch", "sequence", "combined")
  )
  run <- function(operation, variant) {
    if (operation == "graph_assembly") {
      fun <- if (variant == "baseline") {
        batch_env$graph_incremental
      } else {
        batch_env$graph_from_arrays
      }
      return(lapply(arrays, fun))
    }
    if (operation == "complete_tree_parser") {
      return(run_variant(variant, lapply(strings, ns$.parse_iupac_tree_single)))
    }
    run_variant(variant, as_glycan_structure(strings))
  }
  for (operation in names(workloads)) {
    for (variant in workloads[[operation]]) {
      invisible(run(operation, variant))
    }
  }
  rows <- list()
  for (round in 1:3) {
    for (operation in names(workloads)) {
      variants <- workloads[[operation]]
      # Rotate and reverse order to spread position effects over all four variants.
      shift <- (worker + round - 2L) %% length(variants)
      variants <- variants[
        ((seq_along(variants) + shift - 1L) %% length(variants)) + 1L
      ]
      if ((worker + round) %% 2L) {
        variants <- rev(variants)
      }
      for (variant in variants) {
        gc()
        elapsed <- system.time(invisible(run(operation, variant)))[["elapsed"]]
        rows[[length(rows) + 1L]] <- data.frame(
          worker,
          round,
          operation,
          variant,
          elapsed
        )
        cat(worker, round, operation, variant, elapsed, "\n")
      }
    }
  }
  write.csv(
    do.call(rbind, rows),
    file.path(output, paste0("timings-", worker, ".csv")),
    row.names = FALSE
  )
  writeLines(
    trimws(capture.output(sessionInfo()), which = "right"),
    file.path(output, paste0("session-", worker, ".txt"))
  )
} else {
  stop("Use verify or bench WORKER_ID")
}
