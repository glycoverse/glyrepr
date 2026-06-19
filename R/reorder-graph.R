#' Reorder Vertices and Edges
#'
#' @description
#' Reorder the vertices and edges of a list of glycan graphs to be in line
#' with the IUPAC-style sequence.
#'
#' @param graphs A list of igraph graph objects.
#' @returns A list of reordered igraph graph objects.
#' @noRd
reorder_graphs <- function(graphs) {
  checkmate::assert_list(graphs, types = "igraph")
  purrr::map(graphs, .reorder_one_graph)
}

#' Reorder graphs and return indices mapping
#'
#' @param graphs A list of igraph graph objects.
#' @returns A list with 'graphs' (reordered graphs) and 'indices' (original indices).
#' @noRd
reorder_graphs_with_indices <- function(graphs) {
  checkmate::assert_list(graphs, types = "igraph")
  n <- length(graphs)
  reordered <- purrr::map(graphs, .reorder_one_graph)
  list(graphs = reordered, indices = seq_len(n))
}

.reorder_one_graph <- function(graph) {
  root <- which(igraph::degree(graph, mode = "in") == 0)
  seq_cache <- build_seq_cache(graph, root)
  order <- seq_glycan_order(root, seq_cache)
  .reorder_by_sequence_order(graph, order)
}

.reorder_by_sequence_order <- function(graph, sequence_order) {
  target_order <- as.numeric(sequence_order$vertices)
  permutation <- order(target_order)
  graph <- igraph::permute(graph, permutation)
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))

  .permute_edges(graph, as.numeric(sequence_order$edges))
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

  new_g <- igraph::graph_from_data_frame(
    edges,
    directed = igraph::is_directed(g),
    vertices = verts
  )
  new_g$anomer <- g$anomer
  new_g
}
