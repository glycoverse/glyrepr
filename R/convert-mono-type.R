#' Convert the Type of Monosaacharides in a Glycan Graph
#'
#' @description
#' This function converts the all monosaccharides in a glycan graph
#' to a different type.
#' The types are: concrete, generic, and simple.
#' The conversion can only be done from "concrete" to "generic" or "simple",
#' and from "generic" to "simple".
#' Conversion in other orders is not allowed.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param glycan A glycan graph.
#' @param to A character string specifying the target monosaccharide type.
#'  It can be "concrete", "generic", or "simple".
#'
#' @return A glycan graph with monosaccharides converted to the target type.
#'
#' @examples
#' concrete_glycan <- n_glycan_core(mono_type = "concrete")
#' convert_glycan_mono_type(concrete_glycan, to = "generic")
#' convert_glycan_mono_type(concrete_glycan, to = "simple")
#' generic_glycan <- n_glycan_core(mono_type = "generic")
#' convert_glycan_mono_type(generic_glycan, to = "simple")
#'
#' @export
convert_glycan_mono_type <- function(glycan, to) {
  if (!(to %in% c("concrete", "generic", "simple"))) {
    rlang::abort("Must be one of: concrete, generic, simple.")
  }
  from <- decide_glycan_mono_type(glycan)
  valid_from_to(from, to)
  if (inherits(glycan, "ne_glycan_graph")) {
    convert_mono_type_ne(glycan, from, to)
  } else {  # "dn_glycan_graph"
    convert_mono_type_dn(glycan, from, to)
  }
}


#' Decide the Type of Monosaacharides in a Glycan Graph
#'
#' This function trys to decide the type of monosaccharides in a glycan graph.
#'
#' @details
#' By saying "trys" in the description, it means that the function only
#' checks the first monosaccharide in the graph to decide the type.
#' This is reasonable because the monosaccharides in a glycan graph are always
#' of the same type if it was created with the functions in this package.
#' The validation functions will ensure this.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param glycan A glycan graph.
#'
#' @return A character string specifying the monosaccharide type.
#' @export
decide_glycan_mono_type <- function(glycan) {
  stopifnot(is_glycan(glycan))
  if (inherits(glycan, "ne_glycan_graph")) {
    decide_glycan_mono_type_ne(glycan)
  } else {  # "dn_glycan_graph"
    decide_glycan_mono_type_dn(glycan)
  }
}


#' Decide the Type of a Monosaacharide
#'
#' This function decides the type of a monosaccharide name.
#'
#' @details
#' # Three types of monosaccharides
#' There are three types of monosaccharides:
#' - concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.
#' - generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.
#' - simple: e.g. "H", "N", "S", "F".
#'
#' For the full list of monosaccharides, see `glyrepr::monosaccharides`.
#'
#' @param mono A character string specifying a monosaccharide name.
#'
#' @return A character string specifying the monosaccharide type.
#' @export
decide_mono_type <- function(mono) {
  stopifnot(is.character(mono))
  if (mono %in% monosaccharides$concrete) {
    "concrete"
  } else if (mono %in% monosaccharides$generic) {
    "generic"
  } else if (mono %in% monosaccharides$simple) {
    "simple"
  } else {
    rlang::abort("Unknown monosaccharide: {mono}")
  }
}


valid_from_to <- function(from, to) {
  if ((from == "generic" && to == "concrete") || (from == "simple" && to %in% c("concrete", "generic"))) {
    cli::cli_abort(c(
      "Cannot convert from {.val {from}} to {.val {to}}.",
      "i" = "Can only convert in this order: concrete -> generic -> simple."
    ), call = rlang::caller_call())
  }
  if (from == to) {
    cli::cli_abort("It is already {.val {to}}.", call = rlang::caller_call())
  }
}


decide_glycan_mono_type_ne <- function(glycan) {
  first_mono <- igraph::vertex_attr(glycan, "mono")[[1]]
  decide_mono_type(first_mono)
}


decide_glycan_mono_type_dn <- function(glycan) {
  mono_nodes <- which(igraph::V(glycan)$type == "mono")
  mono_names <- igraph::vertex_attr(glycan, "mono")[mono_nodes]
  first_mono <- mono_names[[1]]
  decide_mono_type(first_mono)
}


convert_mono_type_ne <- function(glycan, from, to) {
  from_ <- monosaccharides[[from]]
  to_ <- monosaccharides[[to]]
  new_names <- to_[match(igraph::V(glycan)$mono, from_)]
  igraph::set_vertex_attr(glycan, "mono", value = new_names)
}


convert_mono_type_dn <- function(glycan, from, to) {
  from_ <- monosaccharides[[from]]
  to_ <- monosaccharides[[to]]
  mono_nodes <- which(igraph::V(glycan)$type == "mono")
  old_names <- igraph::vertex_attr(glycan, "mono")[mono_nodes]
  new_names <- to_[match(old_names, from_)]
  igraph::set_vertex_attr(glycan, "mono", value = new_names, index = mono_nodes)
}
