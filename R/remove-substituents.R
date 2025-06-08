#' Remove All Substituents from a Glycan
#'
#' This function replaces all substituents in a glycan structure with empty strings.
#'
#' @param glycan A glycan structure.
#'
#' @return A glycan structure.
#'
#' @examples
#' glycan <- n_glycan_core()
#' igraph::V(glycan)$sub[[1]] <- "3Me"
#' glycan
#' remove_substituents(glycan)
#'
#' @export
remove_substituents <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_structure")
  igraph::set_vertex_attr(glycan, "sub", value = "")
}
