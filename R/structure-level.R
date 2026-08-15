#' Get the Structure Resolution Levels
#'
#' @description
#' Glycan structures can have three possible levels of resolution, determined
#' only by linkage and anomer information:
#' - "intact": No linkage or anomer is unknown or ambiguous.
#' - "partial": At least one linkage or anomer is unknown or ambiguous, and
#'   at least one linkage or anomer contains known information.
#' - "topological": All linkages and anomers are completely unknown
#'   ("??-?"/"??").
#'
#' Floating metadata does not by itself affect the structure level. A floating
#' part's attachment linkage is evaluated like an ordinary linkage, while
#' candidate-parent ambiguity for floating parts and substituents is independent
#' of linkage resolution. Consequently, a glycan with floating metadata is
#' regarded as "intact" when every graph-edge and floating-part attachment
#' linkage, as well as the reducing-end anomer, is fully specified, even though
#' a parent residue remains unlocalized.
#'
#' @param x A [glycan_structure()] vector or a glycan `igraph`.
#'
#' @returns For vector input, a character vector containing one structure level
#'   per element. Missing elements return `NA_character_`, and an empty input
#'   returns `character()`. For graph input, a character scalar.
#'
#' @examples
#' glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
#' get_structure_level(glycan)
#'
#' floating <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' get_structure_level(floating)
#'
#' @seealso [has_linkages()], [has_floating_parts()],
#'   [has_floating_substituents()], [get_mono_type()]
#' @export
get_structure_level <- function(x) {
  if (inherits(x, "igraph")) {
    return(get_graph_structure_level(x))
  }
  checkmate::assert_class(x, "glyrepr_structure")

  if (length(x) == 0) {
    return(character())
  }

  smap_chr(x, get_graph_structure_level)
}

get_graph_structure_level <- function(graph) {
  has_linkages_strict <- .has_linkages_single(graph, strict = TRUE)
  has_linkages_lenient <- .has_linkages_single(graph, strict = FALSE)

  if (has_linkages_strict) {
    return("intact")
  }

  if (has_linkages_lenient) {
    return("partial")
  }

  "topological"
}
