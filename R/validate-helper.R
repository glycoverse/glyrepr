# Is the graph directed?
is_directed_graph <- function(graph) {
  igraph::is_directed(graph)
}


# Is the graph an outward tree?
is_out_tree <- function(graph) {
  igraph::is_tree(graph, mode = "out")
}


# Does the graph have these vertex attributes?
has_vertex_attrs <- function(graph, attrs) {
  all(attrs %in% igraph::vertex_attr_names(graph))
}


# Does the graph have these edge attributes?
has_edge_attrs <- function(graph, attrs) {
  all(attrs %in% igraph::edge_attr_names(graph))
}


.unique_no_na <- function(x) unique(x[!is.na(x)])


# Are all monosaaccharides known?
is_known_mono <- function(monos) {
  known_monos <- c(
    .unique_no_na(monosaccharides$simple),
    .unique_no_na(monosaccharides$generic),
    monosaccharides$concrete
  )
  monos %in% known_monos
}


# Are generic and concrete monosaccharides mixed?
mix_generic_concrete <- function(monos) {
  has_simple <- any(monos %in% .unique_no_na(monosaccharides$simple))
  has_generic <- any(monos %in% .unique_no_na(monosaccharides$generic))
  has_concrete <- any(monos %in% monosaccharides$concrete)
  sum(as.integer(c(has_simple, has_generic, has_concrete))) > 1
}


#' Check if Linkages are Valid
#'
#' Valid linkages are in the form of "a1-2", "b1-4", "a?-1", etc.
#' Specifically, the pattern is `xy-z`:
#' - `x`: the anomer, either "a", "b", or "?".
#' - `y`: the first position, either "1", "2" or "?".
#' - `z`: the second position, either a 1-9 digit or "?".
#' Can also be multiple positions separated by "|", e.g. "1|2|3".
#' "?" could not be used with "|".
#'
#' @param linkages A character vector of linkages.
#'
#' @return A logical vector.
#'
#' @examples
#' # Valid linkages
#' valid_linkages(c("a1-2", "?1-4", "a?-1", "b?-?", "??-?", "a1|2-3"))
#'
#' # Invalid linkages
#' valid_linkages(c("a1-2|?", "1-4", "a|b1-2", "c1-2", "a9-1"))
#'
#' @export
valid_linkages <- function(linkages) {
  if (!is.character(linkages)) {
    rlang::abort("Linkages must be a character vector.")
  }
  anomer_p <- "[ab\\?]"
  pos1_p <- "([12]|\\?)"
  pos2_p <- "([1-9](\\|[1-9])*|\\?)"
  linkage_pattern <- stringr::str_glue("^{anomer_p}{pos1_p}-{pos2_p}$")
  stringr::str_detect(linkages, linkage_pattern)
}
