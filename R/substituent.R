#' Available Substituents
#'
#' Get the available substituents for monosaccharides.
#'
#' @return A character vector.
#' @export
available_substituents <- function() {
  c("Me", "Ac", "NAc", "P", "S", "Pyr", "PC", "PPEtn", "PEtn", "N")
}

#' Remove All Substituents from a Glycan
#'
#' This function replaces all substituents in a glycan structure with empty strings.
#'
#' @param glycan A glyrepr_structure vector.
#'
#' @return A glyrepr_structure vector with all substituents removed.
#' @export
remove_substituents <- function(glycan) {
  if (!is_glycan_structure(glycan)) {
    rlang::abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }
  
  structure_map_structure(glycan, .remove_substituents_single)
}

# Internal function to remove substituents from a single igraph
.remove_substituents_single <- function(glycan) {
  igraph::set_vertex_attr(glycan, "sub", value = "")
}


