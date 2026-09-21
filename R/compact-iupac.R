# Private character-input backend. Public graph APIs retain their validation,
# arbitrary attributes, and trusted-input contracts.
.compact_bliss_labels <- function(edges, colors) {
  graph <- igraph::make_graph(edges, n = length(colors), directed = FALSE)
  as.integer(igraph::canonical_permutation(graph, colors = colors)$labeling)
}

.compact_iupac_arrays <- function(x) {
  residues <- c(monosaccharides$concrete, monosaccharides$generic)
  residues <- unique(residues[!is.na(residues)])
  .compact_parse_native(
    unname(x),
    residues,
    infer_anomer_pos(residues),
    names(unusual_configuration_monosaccharides),
    unname(unusual_configuration_monosaccharides),
    base::order,
    .compact_bliss_labels,
    byte_order = identical(Sys.getlocale("LC_COLLATE"), "C")
  )
}

.compact_build_graph <- function(a) {
  graph <- igraph::make_graph(a$edges, n = length(a$mono), directed = TRUE)
  igraph::vertex_attr(graph) <- list(
    name = as.character(seq_along(a$mono)),
    mono = a$mono,
    sub = a$sub
  )
  igraph::edge_attr(graph) <- list(linkage = a$linkage)
  igraph::graph_attr(graph) <- list(anomer = a$anomer, alditol = a$alditol)
  if (!is.null(a$floating_parts)) {
    graph <- igraph::set_graph_attr(graph, "floating_parts", a$floating_parts)
  }
  if (!is.null(a$floating_substituents)) {
    graph <- igraph::set_graph_attr(
      graph,
      "floating_substituents",
      a$floating_substituents
    )
  }
  graph
}

.compact_structure_from_arrays <- function(x, unique_x, arrays) {
  keys <- vapply(arrays, `[[`, character(1), "iupac")
  unique_keys <- !duplicated(keys)
  graphs <- lapply(arrays[unique_keys], .compact_build_graph)
  names(graphs) <- keys[unique_keys]
  # This private constructor has historically returned unnamed values; the
  # existing vctrs method restores names, including its empty-input semantics.
  new_glycan_structure(unname(keys[match(x, unique_x)]), graphs)
}

.compact_iupac_outcomes <- function(x) {
  arrays <- .compact_iupac_arrays(x)
  graph_cache <- new.env(hash = TRUE, parent = emptyenv())
  lapply(seq_along(arrays), function(i) {
    tryCatch(
      {
        a <- arrays[[i]]
        if (!identical(a$status, "ok")) {
          # Only the failing element takes the reference path under recovery.
          # This preserves its condition class, message, and normalization rules.
          return(process_glycan_structure_element(.parse_iupac_condensed_single(x[[
            i
          ]])))
        }
        key <- a$iupac
        if (!exists(key, graph_cache, inherits = FALSE)) {
          graph_cache[[key]] <- .compact_build_graph(a)
        }
        list(graph = graph_cache[[key]], iupac = key)
      },
      error = \(cnd) cnd
    )
  })
}
