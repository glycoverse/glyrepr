# Source after pkgload::load_all(). No installed/package files are modified.
prototype_env <- new.env(parent = asNamespace("glyrepr"))
Rcpp::sourceCpp(
  "benchmarks/rcpp-prototype/sequence.cpp",
  env = prototype_env,
  cacheDir = file.path(tempdir(), "glyrepr-rcpp")
)
evalq(
  {
    sequence_graph <- function(graph) {
      endpoints <- igraph::as_edgelist(graph, names = FALSE)
      storage.mode(endpoints) <- "integer"
      attrs <- igraph::vertex_attr(graph)
      sequence_cpp(
        endpoints,
        attrs$mono,
        attrs$sub,
        igraph::edge_attr(graph, "linkage")
      )
    }
    baseline_canonical <- canonicalize_graph_with_iupac
    baseline_iupac <- graph_to_iupac
    canonical_cpp <- function(graph) {
      if (
        any(
          c("floating_parts", "floating_substituents") %in%
            igraph::graph_attr_names(graph)
        )
      ) {
        return(baseline_canonical(graph))
      }
      graph <- normalize_alditol_attr(graph)
      graph <- ensure_name_vertex_attr(graph)
      canonical_names <- as.character(seq_len(igraph::vcount(graph)))
      if (!identical(igraph::V(graph)$name, canonical_names)) {
        igraph::V(graph)$name <- canonical_names
      }
      sequence <- sequence_graph(graph)
      iupac <- format_reducing_end_iupac(sequence$iupac, graph)
      if (!is_canonical_sequence_order(sequence)) {
        graph <- .reorder_by_sequence_order(graph, sequence)
      }
      list(graph = graph, iupac = iupac)
    }
    iupac_cpp <- function(graph) {
      checkmate::assert_class(graph, "igraph")
      if (
        any(
          c("floating_parts", "floating_substituents") %in%
            igraph::graph_attr_names(graph)
        )
      ) {
        return(baseline_iupac(graph))
      }
      format_reducing_end_iupac(sequence_graph(graph)$iupac, graph)
    }
  },
  prototype_env
)

with_prototype <- function(code) {
  ns <- asNamespace("glyrepr")
  replacements <- list(
    canonicalize_graph_with_iupac = prototype_env$canonical_cpp,
    graph_to_iupac = prototype_env$iupac_cpp
  )
  originals <- mget(names(replacements), ns)
  set_bindings <- function(values) {
    for (name in names(values)) {
      locked <- bindingIsLocked(name, ns)
      if (locked) {
        unlockBinding(name, ns)
      }
      assign(name, values[[name]], ns)
      if (locked) lockBinding(name, ns)
    }
  }
  on.exit(set_bindings(originals), add = TRUE)
  set_bindings(replacements)
  force(code)
}
