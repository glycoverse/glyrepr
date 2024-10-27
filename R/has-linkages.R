#' Determine if a Glycan Graph has Linkages
#'
#' @description
#' Unknown linkages in a glycan graph are represented by "??-?".
#' This function checks if all linkages in a glycan graph are unknown.
#' Note that even only one linkage is partial known (e.g. "a?-?"),
#' this function will return `TRUE`.
#'
#' @param glycan A glycan graph.
#'
#' @return A logical value indicating if the glycan graph has linkages.
#' @export
has_linkages <- function(glycan) {
  if (!is_glycan(glycan)) {
    stop("Input must be a glycan graph.")
  }
  if (is_ne_glycan(glycan)) {
    has_linkages_ne(glycan)
  } else {
    has_linkages_dn(glycan)
  }
}


has_linkages_ne <- function(glycan) {
  any(igraph::E(glycan)$linkage != "??-?")
}


has_linkages_dn <- function(glycan) {
  any(igraph::V(glycan)$linkage != "??-?", na.rm = TRUE)
}
