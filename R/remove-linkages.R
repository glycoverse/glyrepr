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
  if (is_ne_glycan(glycan)) {
    remove_linkages_ne(glycan)
  } else {
    remove_linkages_dn(glycan)
  }
}


remove_linkages_ne <- function(glycan) {
  igraph::set_edge_attr(glycan, "linkage", value = "??-?")
}


remove_linkages_dn <- function(glycan) {
  linkages <- igraph::V(glycan)$linkage
  linkages[igraph::V(glycan)$type == "linkage"] <- "??-?"
  igraph::set_vertex_attr(glycan, "linkage", value = linkages)
}
