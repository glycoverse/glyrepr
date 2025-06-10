#' Get the Anomeric information
#'
#' @param x A glycan structure vector (glyrepr_structure).
#' @returns a character vector of the anomeric information.
#' @export
#'
#' @examples
#' x <- n_glycan_core()
#' get_anomer(x)
get_anomer <- function(x) {
  structure_map_chr(x, function(g) g$anomer)
}
