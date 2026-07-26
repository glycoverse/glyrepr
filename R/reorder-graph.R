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
  has_parts_attr <- "floating_parts" %in% igraph::graph_attr_names(graph)
  if (has_parts_attr) {
    parts <- igraph::graph_attr(graph, "floating_parts")
    if (length(parts) > 0) {
      return(canonicalize_floating_graph(graph))
    }
    graph <- delete_floating_parts_attr(graph)
  }

  seq_cache <- build_seq_cache(graph)
  root <- seq_cache$root
  order <- seq_glycan_order(root, seq_cache)
  if (is_canonical_sequence_order(order)) {
    return(graph)
  }

  .reorder_by_sequence_order(graph, order)
}

is_canonical_sequence_order <- function(sequence_order) {
  vertex_order <- as.numeric(sequence_order$vertices)
  edge_order <- as.numeric(sequence_order$edges)

  identical(base::order(vertex_order), seq_along(vertex_order)) &&
    identical(edge_order, as.numeric(seq_along(edge_order)))
}

canonicalize_graph_with_iupac <- function(graph) {
  if ("floating_parts" %in% igraph::graph_attr_names(graph)) {
    graph <- canonicalize_glycan_graph(graph)
    return(list(graph = graph, iupac = graph_to_iupac(graph)))
  }

  graph <- ensure_name_vertex_attr(graph)
  canonical_names <- as.character(seq_len(igraph::vcount(graph)))
  if (!identical(igraph::V(graph)$name, canonical_names)) {
    igraph::V(graph)$name <- canonical_names
  }

  seq_cache <- build_seq_cache(graph)
  root <- seq_cache$root
  sequence <- seq_glycan_order_iupac(root, seq_cache)
  iupac <- paste0(
    sequence$iupac,
    "(",
    graph$anomer,
    "-"
  )
  if (!is_canonical_sequence_order(sequence)) {
    graph <- .reorder_by_sequence_order(graph, sequence)
  }

  list(graph = graph, iupac = iupac)
}

.reorder_by_sequence_order <- function(graph, sequence_order) {
  target_order <- as.numeric(sequence_order$vertices)
  permutation <- order(target_order)
  graph <- igraph::permute(graph, permutation)
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))

  .permute_edges(graph, as.numeric(sequence_order$edges))
}

.permute_edges <- function(g, order) {
  graph_attrs <- igraph::graph_attr(g)
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
  for (attr_name in names(graph_attrs)) {
    new_g <- igraph::set_graph_attr(
      new_g,
      attr_name,
      value = graph_attrs[[attr_name]]
    )
  }
  new_g
}
