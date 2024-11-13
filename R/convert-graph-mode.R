#' Convert between NE and DN glycan graphs
#'
#' @description
#' This function converts a Node-Edge (NE) glycan graph into a Dual-Node (DN) glycan graph,
#' or vice versa.
#' In an NE glycan graph, nodes represent monosaccharides and edges represent linkages.
#' In a DN glycan graph, both monosaccharides and linkages are nodes, and they alternate.
#' See [as_glycan_graph()] for details.
#'
#' @param glycan A glycan graph.
#' @param to The target graph mode, either "ne" or "dn".
#' @param strict If `TRUE`, an error will be thrown if the glycan is already
#' in the target mode. If `FALSE`, the glycan will be returned as is.
#' Default is `TRUE`.
#'
#' @return The converted glycan graph object.
#'
#' @examples
#' (glycan <- n_glycan_core(mode = "ne"))
#'
#' # Convert to DN mode
#' convert_graph_mode(glycan, to = "dn")
#'
#' # Converting to NE mode will raise an error
#' try(convert_graph_mode(glycan, to = "ne"))
#'
#' # Set `strict` to `FALSE` to return the glycan unchanged
#' convert_graph_mode(glycan, to = "ne", strict = FALSE)
#'
#' @export
convert_graph_mode <- function(glycan, to, strict = TRUE) {
  checkmate::assert_class(glycan, "glycan_graph")
  checkmate::assert_choice(to, c("ne", "dn"))
  checkmate::assert_flag(strict)

  if (to == "ne" && is_dn_glycan(glycan)) {
    return(convert_dn_to_ne(glycan))
  }
  if (to == "dn" && is_ne_glycan(glycan)) {
    return(convert_ne_to_dn(glycan))
  }
  if (strict) {
    cli::cli_abort("The glycan is already in {.val {stringr::str_to_upper(to)}} mode.")
  } else {
    return(glycan)
  }
}


convert_ne_to_dn <- function(glycan) {
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


convert_dn_to_ne <- function(glycan) {
  # Deal with edge case that only one monosaccharide exists
  if (igraph::vcount(glycan) == 1) {
    new_graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
    igraph::V(new_graph)$mono <- igraph::V(glycan)$mono
    igraph::V(new_graph)$sub <- igraph::V(glycan)$sub
    igraph::E(new_graph)$linkage <- character(0)
    new_graph$anomer <- glycan$anomer
    new_graph$alditol <- glycan$alditol
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
