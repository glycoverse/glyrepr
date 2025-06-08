#' Determine if a Glycan Structure has Linkages
#'
#' @description
#' Unknown linkages in a glycan structure are represented by "??-?".
#' This function checks if all linkages in a glycan structure are unknown.
#' Note that even only one linkage is partial known (e.g. "a?-?"),
#' this function will return `TRUE`.
#'
#' @param glycan A glycan structure.
#'
#' @return A logical value indicating if the glycan structure has linkages.
#'
#' @examples
#' glycan <- o_glycan_core_1(linkage = TRUE)
#' has_linkages(glycan)
#' print(glycan)
#'
#' glycan <- remove_linkages(glycan)
#' has_linkages(glycan)
#' print(glycan)
#'
#' @seealso [remove_linkages()], [possible_linkages()]
#'
#' @export
has_linkages <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_structure")
  any(igraph::E(glycan)$linkage != "??-?")
}
