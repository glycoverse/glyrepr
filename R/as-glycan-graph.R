#' Convert igraph Graph to Glycan Graph
#'
#' @description
#' A glycan graph is a subclass of igraph graph with additional constraints.
#' This function checks these constraints and append the S3 class `glycan_graph`
#' to the graph object.
#'
#' Two types of glycan graphs are supported: dual-node (DN) and node-edge (NE).
#' A NE glycan graph is a directed modeling of a glycan structure,
#' where nodes are monosaccharides and edges are linkages.
#' A DN glycan graph is another way to model a glycan structure,
#' where both monosaccharides and linkages are nodes.
#' Monosaccharides and linkages alternate in a DN glycan graph.
#' For more details, see "details" section.
#'
#' Unknown linkages could either be assigned as NA or "??-?".
#' NA will be converted to "??-?" internally.
#'
#' `as_dn_glycan_graph()` is the same as `as_glycan_graph()` with `type = "dn"`.
#' `as_ne_glycan_graph()` is the same as `as_glycan_graph()` with `type = "ne"`.
#'
#' @details
#' # S3 Classes
#'
#' A glycan graph is merely an igraph graph with an additional S3 class `glycan_graph`.
#' Based on type of the glycan, `ne_glycan_graph` or `dn_glycan_graph` is also
#' appended to the graph.
#' Therefore, `sloop::s3_class()` of a glycan graph is either
#' `c("ne_glycan_graph", "glycan_graph", "igraph")` or
#' `c("dn_glycan_graph", "glycan_graph", "igraph")`.
#'
#' # Node-Edge (NE) Glycan Graph
#'
#' This is the most intuitive way to model a glycan structure.
#' Nodes represent monosaccharides and edges represent linkages.
#'
#' Constraints:
#' - The graph must be directed and an outward tree (reducing end as root).
#' - The graph must have a vertex attribute `mono` for monosaccharide names.
#' - The graph must have an edge attribute `linkage` for linkages.
#' - Monosaccharide name must be known, either generic (Hex, HexNAc, etc.)
#'   or concrete (Glc, Gal, etc.), but not a mixture of both.
#'   NA is not allowed.
#' - Linkages must be valid, in the form of "a/bX-Y", where X and Y are integers,
#'   e.g. "b1-4", "a2-3", etc.
#'   NA is allowed for unknown linkages.
#'
#' # Dual-Node (DN) Glycan Graph
#'
#' This is another way to model a glycan structure.
#' In this mode, both monosaccharides and linkages are nodes,
#' and they alternate in the graph to represent the structure.
#' Python library `glycowork` and `glypy` use this model.
#' It is not that intuitive, but allows a more flexible motif matching.
#'
#' Constraints:
#' - The graph must be directed and an outward tree (reducing end as root).
#' - The graph must have vertex attributes `type`, `mono`, and `linkage`.
#'   `type` should be either "mono" or "linkage".
#'   NA is not allowed for `type`.
#'   `mono` and `linkage` are monosaccharide name and linkage, respectively.
#'   For nodes with `type = "mono"`, `linkage` is always NA.
#'   For nodes with `type = "linkage"`, `mono` is always NA.
#' - Monosaccharide name must be known, either generic or concrete, but not a mixture of both.
#'   For nodes with `type = "mono"`, NA is not allowed for `mono`.
#' - Linkages must be valid.
#'   For nodes with `type = "linkage"`, NA is allowed for `linkage`.
#' - `mono` nodes and `linkage` nodes must alternate in the graph.
#'   In another word, each edge must connect a `mono` node and a `linkage` node.
#' - The outermost nodes (nodes with out-degree 0) must be `mono` nodes,
#'   not `linkage` nodes.
#'
#' @param graph An igraph graph object.
#' @param type A character string specifying the glycan graph type.
#'   It can be "auto", "ne", or "dn".
#'   If "auto" (the default), the function tries to infer the type,
#'   and raises an error if it cannot.
#'   If "ne" or "dn", the function assumes the graph is of the specified type.
#'
#' @return A glycan graph object.
#'
#' @examples
#' library(igraph)
#'
#' # A simple NE glycan graph: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2)  # 1 and 2 are monosaccharides
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")
#' E(graph)$linkage <- "b1-4"
#' as_glycan_graph(graph)
#'
#' # A simple DN glycan graph: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2, 2-+3)  # 1 and 3 are monosaccharides, 2 is linkage
#' V(graph)$type <- c("mono", "linkage", "mono")
#' V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
#' V(graph)$linkage <- c(NA, "b1-4", NA)
#' as_glycan_graph(graph, type = "dn")
#'
#' @export
as_glycan_graph <- function(graph, type = "auto") {
  stopifnot(igraph::is_igraph(graph))
  stopifnot(type %in% c("auto", "ne", "dn"))

  if (type == "auto") {
    glycan <- new_dn_glycan_graph(graph)
    tryCatch(
      validate_dn_glycan_graph(glycan),
      error = function(e) {
        glycan <- new_ne_glycan_graph(graph)
        tryCatch(
          validate_ne_glycan_graph(glycan),
          error = function(e) {
            rlang::abort(paste0(
              "Could not infer glycan graph type. ",
              "Please specify 'type' argument manually to ",
              "see potential problems in the graph."))
          }
        )
      }
    )
  } else if (type == "ne") {
    as_ne_glycan_graph(graph)
  } else {  # type == "dn"
    as_dn_glycan_graph(graph)
  }
}


#' @rdname as_glycan_graph
#' @importFrom magrittr %>%
#' @export
as_dn_glycan_graph <- function(graph) {
  stopifnot(igraph::is_igraph(graph))
  graph %>%
    new_dn_glycan_graph() %>%
    validate_dn_glycan_graph() %>%
    ensure_name_vertex_attr() %>%
    clean_dn_linkages()
}


#' @rdname as_glycan_graph
#' @importFrom magrittr %>%
#' @export
as_ne_glycan_graph <- function(graph) {
  stopifnot(igraph::is_igraph(graph))
  graph %>%
    new_ne_glycan_graph() %>%
    validate_ne_glycan_graph() %>%
    ensure_name_vertex_attr() %>%
    clean_ne_linkages()
}


# Add "name" vertex attributes if missing.
ensure_name_vertex_attr <- function(glycan) {
  if (!("name" %in% igraph::vertex_attr_names(glycan))) {
    names <- as.character(seq_len(igraph::vcount(glycan)))
    glycan <- igraph::set_vertex_attr(glycan, "name", value = names)
  }
  glycan
}


# Replace NA in linkages to "??-?".
clean_ne_linkages <- function(glycan) {
  linkages <- igraph::edge_attr(glycan, "linkage")
  linkages[is.na(linkages)] <- "??-?"
  igraph::set_edge_attr(glycan, "linkage", value = linkages)
}


clean_dn_linkages <- function(glycan) {
  types <- igraph::vertex_attr(glycan, "type")
  linkages <- igraph::vertex_attr(glycan, "linkage")
  linkages[(types == "linkage") & is.na(linkages)] <- "??-?"
  igraph::set_vertex_attr(glycan, "linkage", value = linkages)
}
