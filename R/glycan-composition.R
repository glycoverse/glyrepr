#' Get Compositions of Glycans
#'
#' Get the compositions of a glycan graph or a list of glycan graphs.
#' The composition is a named integer vector.
#'
#' @param glycan A glycan graph.
#' @return A named integer vector.
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
