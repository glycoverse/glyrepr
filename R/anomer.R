#' Get the Anomeric information
#'
#' @param x A glycan structure vector or a glycan `igraph`.
#' @returns A character vector for structure-vector input, or a character scalar
#'   for graph input.
#' @export
#'
#' @examples
#' x <- n_glycan_core()
#' get_anomer(x)
get_anomer <- function(x) {
  if (inherits(x, "igraph")) {
    return(igraph::graph_attr(x, "anomer"))
  }
  smap_chr(x, function(g) g$anomer)
}
