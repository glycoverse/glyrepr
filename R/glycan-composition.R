#' Get Compositions of Glycans
#'
#' Get the compositions of a glycan graph or a list of glycan graphs.
#' The composition is a named integer vector.
#'
#' @param glycan A glycan graph.
#' @param glycans A list of glycan graphs.
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
  stopifnot(is_glycan(glycan))
  monos <- igraph::V(glycan)$mono
  result_tb <- table(monos)
  result <- as.integer(result_tb)
  names(result) <- names(result_tb)
  result
}


#' @rdname get_composition
#' @export
get_compositions <- function(glycans) {
  purrr::map(glycans, get_composition)
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
#' @seealso [get_composition()], [get_compositions()]
#'
#' @export
count_monos <- function(glycan) {
  stopifnot(is_glycan(glycan))
  sum(get_composition(glycan))
}
