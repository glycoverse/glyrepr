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
#' Two types of glycan graphs are supported: dual-node (DN) and node-edge (NE).
#' A NE glycan graph is a directed modeling of a glycan structure,
#' where nodes are monosaccharides and edges are linkages.
#' A DN glycan graph is another way to model a glycan structure,
#' where both monosaccharides and linkages are nodes.
#' Monosaccharides and linkages alternate in a DN glycan graph.
#' For more details, see "details" section.
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
#' - The graph must have a graph attribute `anomer`, in the form of "a1".
#' - The graph must have a graph attribute `alditol`, a logical value.
#' - The graph must have vertex attributes `type`, `mono`, `sub`, and `linkage`.
#'   `type` should be either "mono" or "linkage".
#'   NA is not allowed for `type`.
#'   `mono` and `linkage` are monosaccharide name and linkage, respectively.
#'   For nodes with `type = "mono"`, `linkage` is always NA.
#'   For nodes with `type = "linkage"`, `mono` and `sub` are always NA.
#' - Monosaccharide name must be known, either generic or concrete, but not a mixture of both.
#'   For nodes with `type = "mono"`, NA is not allowed for `mono`.
#' - Substituent must be valid.
#'   For nodes with `type = "mono"`, NA is not allowed for `sub`.
#' - Linkages must be valid.
#'   For nodes with `type = "linkage"`, NA is not allowed for `linkage`.
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
#' V(graph)$sub <- ""
#' E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' as_glycan_graph(graph)
#'
#' # A simple DN glycan graph: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2, 2-+3)  # 1 and 3 are monosaccharides, 2 is linkage
#' V(graph)$type <- c("mono", "linkage", "mono")
#' V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
#' V(graph)$sub <- c("", NA, "")
#' V(graph)$linkage <- c(NA, "b1-4", NA)
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' as_glycan_graph(graph, type = "dn")
#'
#' @export
as_glycan_graph <- function(graph, type = "auto") {
  checkmate::assert_class(graph, "igraph")
  checkmate::assert_choice(type, c("auto", "ne", "dn"))

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
  checkmate::assert_class(graph, "igraph")
  graph %>%
    new_dn_glycan_graph() %>%
    validate_dn_glycan_graph() %>%
    ensure_name_vertex_attr()
}


#' @rdname as_glycan_graph
#' @importFrom magrittr %>%
#' @export
as_ne_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  graph %>%
    new_ne_glycan_graph() %>%
    validate_ne_glycan_graph() %>%
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
