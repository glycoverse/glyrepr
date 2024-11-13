#' Get Composition of a Glycan
#'
#' Get the composition of a glycan graph.
#' The composition is a named integer vector.
#'
#' @param glycan A glycan graph.
#'
#' @return A named integer vector.
#'
#' @examples
#' glycan <- n_glycan_core()
#' get_composition(glycan)
#'
#' @seealso [count_monos()]
#'
#' @export
get_composition <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_graph")
  monos <- igraph::V(glycan)$mono
  result_tb <- table(monos)
  result <- as.integer(result_tb)
  names(result) <- names(result_tb)
  result
}


#' Count the Number of Monosaccharides
#'
#' This function returns the total number of monosaccharides of a
#' glycan graph.
#'
#' @param glycan A glycan graph.
#'
#' @return An integer.
#'
#' @examples
#' glycan <- n_glycan_core()
#' count_monos(glycan)
#'
#' @seealso [get_composition()]
#'
#' @export
count_monos <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_graph")
  sum(get_composition(glycan))
}
