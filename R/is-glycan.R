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
