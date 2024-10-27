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
#' @export
glycan_graph_mode <- function(x) {
  if (!is_glycan(x)) {
    stop("Input must be a glycan graph.")
  }
  if (is_ne_glycan(x)) "ne" else "dn"
}
