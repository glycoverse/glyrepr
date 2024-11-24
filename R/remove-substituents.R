#' Remove All Substituents from a Glycan
#'
#' This function replaces all substituents in a glycan graph with empty strings.
#'
#' @param glycan A glycan graph.
#'
#' @return A glycan graph.
#'
#' @examples
#' glycan <- n_glycan_core()
#' igraph::V(glycan)$sub[[1]] <- "3Me"
#' glycan
#' remove_substituents(glycan)
#'
#' @export
remove_substituents <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_graph")
  if (is_ne_glycan(glycan)) {
    remove_substituents_ne(glycan)
  } else {
    remove_substituents_dn(glycan)
  }
}


remove_substituents_ne <- function(glycan) {
  igraph::set_vertex_attr(glycan, "sub", value = "")
}


remove_substituents_dn <- function(glycan) {
  subs <- igraph::V(glycan)$sub
  subs[igraph::V(glycan)$type == "mono"] <- ""
  igraph::set_vertex_attr(glycan, "sub", value = subs)
}
