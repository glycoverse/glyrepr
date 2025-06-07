#' Convert igraph Graph to Glycan Graph
#'
#' @description
#' A glycan graph is a subclass of igraph graph with additional constraints.
#' This function checks these constraints and append the S3 class `glycan_graph`
#' to the graph object.
#'
#' Generally, **it is not recommended** to create a glycan graph using this
#' function manually.
#' Use the [glyparse](https://github.com/glycoverse/glyparse) package
#' to generate a glycan graph from a structure string.
#'
#' A glycan graph is a directed modeling of a glycan structure,
#' where nodes are monosaccharides and edges are linkages.
#'
#' @details
#' # S3 Classes
#'
#' A glycan graph is merely an igraph graph with an additional S3 class `glycan_graph`.
#' Therefore, `sloop::s3_class()` of a glycan graph is 
#' `c("glycan_graph", "igraph")`.
#'
#' Constraints:
#' - The graph must be directed and an outward tree (reducing end as root).
#' - The graph must have a graph attribute `anomer`, in the form of "a1".
#' - The graph must have a graph attribute `alditol`, a logical value.
#' - The graph must have a vertex attribute `mono` for monosaccharide names.
#' - The graph must have a vertex attribute `sub` for substituents.
#' - The graph must have an edge attribute `linkage` for linkages.
#' - Monosaccharide name must be known, either generic (Hex, HexNAc, etc.)
#'   or concrete (Glc, Gal, etc.), but not a mixture of both.
#'   NA is not allowed.
#' - Substituent must be valid, in the form of "xY", where x is position
#'   and Y is substituent name, e.g. "2Ac", "3S", etc.
#' - Linkages must be valid, in the form of "a/bX-Y", where X and Y are integers,
#'   e.g. "b1-4", "a2-3", etc.
#'   NA is not allowed.
#'
#' @param graph An igraph graph object.
#'
#' @return A glycan graph object.
#'
#' @examples
#' library(igraph)
#'
#' # A simple glycan graph: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2)  # 1 and 2 are monosaccharides
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")
#' V(graph)$sub <- ""
#' E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' as_glycan_graph(graph)
#'
#' @importFrom magrittr %>%
#' @export
as_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  
  graph %>%
    new_glycan_graph() %>%
    validate_glycan_graph() %>%
    ensure_name_vertex_attr()
}


# Add "name" vertex attributes if missing.
ensure_name_vertex_attr <- function(glycan) {
  if (!("name" %in% igraph::vertex_attr_names(glycan))) {
    names <- as.character(seq_len(igraph::vcount(glycan)))
    glycan <- igraph::set_vertex_attr(glycan, "name", value = names)
  }
  glycan
}
