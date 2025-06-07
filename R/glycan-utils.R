#' Get Composition of a Glycan
#'
#' Get the composition of a glycan graph.
#' The composition is returned as a glyrepr_composition object.
#'
#' @param glycan A glycan graph.
#'
#' @return A glyrepr_composition object.
#'
#' @examples
#' glycan <- n_glycan_core()
#' get_composition(glycan)
#' 
#' @export
get_composition <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_graph")
  monos <- igraph::V(glycan)$mono
  result_tb <- table(monos)
  result <- as.integer(result_tb)
  names(result) <- names(result_tb)
  composition(result)
}
