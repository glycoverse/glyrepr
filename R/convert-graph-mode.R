#' Convert between NE and DN glycan graphs
#'
#' @description
#' This function converts a Node-Edge (NE) glycan graph into a Dual-Node (DN) glycan graph,
#' or vice versa.
#' In an NE glycan graph, nodes represent monosaccharides and edges represent linkages.
#' In a DN glycan graph, both monosaccharides and linkages are nodes, and they alternate.
#' See [as_glycan_graph()] for details.
#'
#' @param glycan A glycan graph object of class `glycan_graph`.
#'
#' @return The converted glycan graph object.
#'
#' @examples
#' library(igraph)
#' # Create a simple NE glycan graph: GlcNAc(b1-4)GlcNAc
#' ne_graph <- make_graph(~ 1-+2)
#' V(ne_graph)$mono <- c("GlcNAc", "GlcNAc")
#' E(ne_graph)$linkage <- "b1-4"
#' ne_glycan <- as_glycan_graph(ne_graph)
#'
#' # Convert NE glycan graph to DN glycan graph
#' dn_glycan <- convert_ne_to_dn(ne_glycan)
#'
#' # Inspect the DN glycan graph
#' print(dn_glycan, verbose = TRUE)
#'
#' @export
convert_ne_to_dn <- function(glycan) {
  if (!inherits(glycan, "ne_glycan_graph")) {
    rlang::abort("Input must be an NE glycan graph.")
  }

  # Deal with edge case that only one monosaccharide exists
  if (igraph::vcount(glycan) == 1) {
    new_graph <- glycan
    class(new_graph) <- "igraph"
    igraph::V(new_graph)$type <- "mono"
    igraph::V(new_graph)$linkage <- NA_character_
    return(new_dn_glycan_graph(new_graph))
  }

  # Get edge list data frame
  edge_list <- igraph::as_data_frame(glycan, what = "edges")
  n_nodes <- igraph::vcount(glycan)

  # Initiate new graph
  new_graph <- glycan
  class(new_graph) <- "igraph"
  igraph::V(new_graph)$type <- "mono"
  igraph::V(new_graph)$linkage <- NA_character_
  new_graph <- igraph::delete_edge_attr(new_graph, "linkage")

  # Iteratively add new linkage nodes
  for (i in 1:nrow(edge_list)) {
    from = edge_list$from[i]
    to = edge_list$to[i]
    linkage = edge_list$linkage[i]
    name <- as.character(n_nodes + i)
    new_graph <- igraph::add_vertices(
      new_graph,
      nv = 1,
      name = name,
      type = "linkage",
      linkage = linkage,
      mono = NA_character_
    )
    new_graph <- igraph::add_edges(new_graph, c(from, name, name, to))
    new_graph <- igraph::delete_edges(new_graph, paste0(from, "|", to))
  }

  new_dn_glycan_graph(new_graph)
}


#' @rdname convert_ne_to_dn
#' @export
convert_dn_to_ne <- function(glycan) {
  if (!inherits(glycan, "dn_glycan_graph")) {
    rlang::abort("Input must be a DN glycan graph.")
  }

  # Deal with edge case that only one monosaccharide exists
  if (igraph::vcount(glycan) == 1) {
    new_graph <- igraph::make_graph(~ 1)
    igraph::V(new_graph)$mono <- igraph::V(glycan)$mono
    return(new_ne_glycan_graph(new_graph))
  }

  # Get edge list and vertex list
  edge_list <- igraph::as_data_frame(glycan, what = "edges")
  vertex_list <- igraph::as_data_frame(glycan, what = "vertices")
  linkage_vertices <- vertex_list[vertex_list$type == "linkage", c("name", "linkage")]

  # Initiate new graph
  new_graph <- glycan
  class(new_graph) <- "igraph"
  new_graph <- igraph::delete_vertices(new_graph, which(igraph::V(new_graph)$type == "linkage"))
  new_graph <- igraph::delete_vertex_attr(new_graph, "linkage")
  new_graph <- igraph::delete_vertex_attr(new_graph, "type")

  # Iteratively add new edges
  for (i in 1:nrow(linkage_vertices)) {
    name <- as.character(linkage_vertices$name[i])
    linkage <- linkage_vertices$linkage[i]
    from <- edge_list$from[edge_list$to == name]
    to <- edge_list$to[edge_list$from == name]
    new_graph <- igraph::add_edges(new_graph, c(from, to), linkage = linkage)
  }
  igraph::V(new_graph)$name <- as.integer(1:igraph::vcount(new_graph))

  new_ne_glycan_graph(new_graph)
}
