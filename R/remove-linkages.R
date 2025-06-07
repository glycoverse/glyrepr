#' Remove All Linkages from a Glycan Graph
#'
#' @description
#' An unkonwn linkage in a glycan graph is represented by "??-?".
#' This function replaces all linkages in a glycan graph with "??-?",
#' so that `has_linkages()` will return `FALSE`.
#'
#' @param glycan A glycan graph.
#'
#' @return A glycan graph with all linkages removed.
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
  checkmate::assert_class(glycan, "glycan_graph")
  igraph::set_edge_attr(glycan, "linkage", value = "??-?")
}
