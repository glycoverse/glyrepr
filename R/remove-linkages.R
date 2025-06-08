#' Remove All Linkages from a Glycan Structure
#'
#' @description
#' An unkonwn linkage in a glycan structure is represented by "??-?".
#' This function replaces all linkages in a glycan structure with "??-?",
#' so that `has_linkages()` will return `FALSE`.
#'
#' @param glycan A glycan structure.
#'
#' @return A glycan structure with all linkages removed.
#'
#' @examples
#' glycan <- n_glycan_core()
#' print(glycan)
#'
#' glycan <- remove_linkages(glycan)
#' print(glycan)
#'
#' @seealso [has_linkages()], [possible_linkages()]
#'
#' @export
remove_linkages <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_structure")
  igraph::set_edge_attr(glycan, "linkage", value = "??-?")
}
