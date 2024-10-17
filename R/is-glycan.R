#' Check if an Object is a Glycan Graph
#'
#' This function checks if an object is a glycan graph.
#' That is an object of class `glycan_graph`.
#'
#' @param x An object to be checked.
#'
#' @return A logical value.
#' @export
is_glycan <- function(x) {
  inherits(x, "glycan_graph")
}
