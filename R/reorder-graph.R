#' Reorder Vertices and Edges
#'
#' @description
#' Reorder the vertices and edges of all underlying graphs of a [glycan_structure()] object,
#' to be in line with the IUPAC-style sequence.
#'
#' @param x A [glycan_structure()] object.
#' @returns A [glycan_structure()] object with the underlying graphs reordered.
#' @noRd
reorder_graphs <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")
  old_structures <- attr(x, "structures")
  new_structures <- purrr::map(old_structures, .reorder_one_graph)
  attr(x, "structures") <- new_structures
  x
}

.reorder_one_graph <- function(graph) {
  root <- which(igraph::degree(graph, mode = "in") == 0)
  seq_cache <- build_seq_cache(graph, root)
  pseudo_seq <- seq_glycan(root, seq_cache)
  .reorder_by_pseudo_seq(graph, pseudo_seq)
}

.reorder_by_pseudo_seq <- function(graph, pseudo_seq) {
  # Reorder the vertices
  vertices <- stringr::str_extract_all(pseudo_seq, "V(\\d+)")[[1]]
  vertices <- stringr::str_sub(vertices, 2L, -1L)
  # Fix: Use the correct permutation vector
  # vertices contains the desired order, we need the inverse permutation
  target_order <- as.numeric(vertices)
  permutation <- order(target_order)
  graph <- igraph::permute(graph, permutation)
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))

  # Reorder the edges
  edges <- stringr::str_extract_all(pseudo_seq, "E(\\d+)")[[1]]
  edges <- stringr::str_sub(edges, 2L, -1L)
  graph <- .permute_edges(graph, as.numeric(edges))

  graph
}

.permute_edges <- function(g, order) {
  edges <- igraph::as_data_frame(g, what = "edges")
  verts <- igraph::as_data_frame(g, what = "vertices")
  edges <- edges[order, , drop = FALSE]
  
  # Ensure 'name' column is first in vertices data frame if it exists
  # This is required by igraph::graph_from_data_frame() to avoid "Duplicate vertex names" error
  if ("name" %in% names(verts)) {
    name_col <- verts[["name"]]
    other_cols <- verts[, !names(verts) %in% "name", drop = FALSE]
    verts <- cbind(name = name_col, other_cols)
  }
  
  new_g <- igraph::graph_from_data_frame(edges, directed = igraph::is_directed(g), vertices = verts)
  new_g$anomer <- g$anomer
  new_g
}
