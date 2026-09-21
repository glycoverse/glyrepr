# Source after pkgload::load_all(). This is deliberately not a package backend.
compact_dir <- file.path("benchmarks", "compact-cpp")
Rcpp::sourceCpp(
  file.path(compact_dir, "compact.cpp"),
  cacheDir = file.path(tempdir(), "compact-cpp-cache")
)
compact_residues <- unique(na.omit(c(
  monosaccharides$concrete,
  monosaccharides$generic
)))
compact_positions <- infer_anomer_pos(compact_residues)

compact_bliss <- function(edges, colors) {
  graph <- igraph::make_graph(edges, n = length(colors), directed = FALSE)
  as.integer(igraph::canonical_permutation(graph, colors = colors)$labeling)
}

compact_arrays <- function(x) {
  compact_parse(
    x,
    compact_residues,
    compact_positions,
    names(unusual_configuration_monosaccharides),
    unname(unusual_configuration_monosaccharides),
    base::order,
    compact_bliss,
    byte_order = identical(Sys.getlocale("LC_COLLATE"), "C")
  )
}

compact_graph <- function(a) {
  g <- igraph::make_graph(a$edges, n = length(a$mono), directed = TRUE)
  igraph::vertex_attr(g) <- list(
    name = as.character(seq_along(a$mono)),
    mono = a$mono,
    sub = a$sub
  )
  igraph::edge_attr(g) <- list(linkage = a$linkage)
  igraph::graph_attr(g) <- list(anomer = a$anomer, alditol = a$alditol)
  if (!is.null(a$floating_parts)) {
    g <- igraph::set_graph_attr(g, "floating_parts", a$floating_parts)
  }
  if (!is.null(a$floating_substituents)) {
    g <- igraph::set_graph_attr(
      g,
      "floating_substituents",
      a$floating_substituents
    )
  }
  g
}

compact_structure <- function(x, fallback = TRUE) {
  stopifnot(is.character(x))
  unique_x <- unique(unname(x[!is.na(x)]))
  arrays <- compact_arrays(unique_x)
  keys <- character(length(arrays))
  graphs <- list()
  for (i in seq_along(arrays)) {
    a <- arrays[[i]]
    if (a$status != "ok") {
      if (!fallback) {
        stop(a$reason)
      }
      # Reference fallback supplies public diagnostics for native errors/size limits.
      # Coverage and parity audits disable it for every valid input.
      value <- as_glycan_structure(unique_x[i])
      key <- as.character(value)[[1]]
      graph <- attr(value, "graphs")[[key]]
    } else {
      key <- a$iupac
      graph <- if (is.null(graphs[[key]])) compact_graph(a) else NULL
    }
    keys[i] <- key
    if (is.null(graphs[[key]])) graphs[[key]] <- graph
  }
  iupac <- keys[match(x, unique_x)]
  names(iupac) <- names(x)
  new_glycan_structure(iupac, graphs)
}
