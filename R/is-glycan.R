#' Check if an Object is a Glycan Graph
#'
#' These functions check if an object is a glycan graph
#' (an object of class `glycan_graph`),
#' a node-edge glycan graph (an object of class `ne_glycan_graph`),
#' or a dual-node glycan graph (an object of class `dn_glycan_graph`).
#' See docs for [as_glycan_graph()] for more details.
#'
#' @param x An object to be checked.
#'
#' @return A logical value.
#'
#' @seealso [as_glycan_graph()], [glycan_graph_mode()]
#'
#' @export
is_glycan <- function(x) {
  inherits(x, "glycan_graph")
}


#' @rdname is_glycan
#' @export
is_ne_glycan <- function(x) {
  inherits(x, "ne_glycan_graph")
}


#' @rdname is_glycan
#' @export
is_dn_glycan <- function(x) {
  inherits(x, "dn_glycan_graph")
}


#' Determine the Glycan Graph Mode
#'
#' These functions determine the mode of a glycan graph.
#' The mode can be either "ne" (node-edge) or "dn" (dual-node).
#' See docs for [as_glycan_graph()] for more details.
#'
#' @param x A glycan graph.
#'
#' @return A character string.
#'
#' @examples
#' glycan <- o_glycan_core_1(mode = "ne")
#' glycan_graph_mode(glycan)
#'
#' glycan <- o_glycan_core_1(mode = "dn")
#' glycan_graph_mode(glycan)
#'
#' @seealso [as_glycan_graph()], [is_glycan()], [is_ne_glycan()], [is_dn_glycan()]
#'
#' @export
glycan_graph_mode <- function(x) {
  checkmate::assert_class(x, "glycan_graph")
  if (is_ne_glycan(x)) "ne" else "dn"
}
