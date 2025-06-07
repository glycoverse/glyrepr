#' Check if an Object is a Glycan Graph
#'
#' This function checks if an object is a glycan graph
#' (an object of class `glycan_graph`).
#' See docs for [as_glycan_graph()] for more details.
#'
#' @param x An object to be checked.
#'
#' @return A logical value.
#'
#' @seealso [as_glycan_graph()]
#'
#' @export
is_glycan <- function(x) {
  inherits(x, "glycan_graph")
}
