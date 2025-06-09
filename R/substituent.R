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
#' @param glycan An igraph object representing a glycan structure, or a glyrepr_structure vector.
#'
#' @return An igraph object representing a glycan structure, or a glyrepr_structure vector.
#' @export
remove_substituents <- function(glycan) {
  UseMethod("remove_substituents")
}

#' @export
remove_substituents.igraph <- function(glycan) {
  igraph::set_vertex_attr(glycan, "sub", value = "")
}

#' @export
remove_substituents.glyrepr_structure <- function(glycan) {
  # Extract individual igraph objects and apply remove_substituents to each
  data <- vctrs::vec_data(glycan)
  codes <- vctrs::field(data, "codes")
  structures <- attr(glycan, "structures")
  
  # Remove substituents from each unique structure
  modified_structures <- purrr::map(structures, ~ remove_substituents.igraph(.x))
  names(modified_structures) <- names(structures)
  
  # Update the stored structures
  attr(glycan, "structures") <- modified_structures
  
  # Regenerate IUPAC codes for the modified structures
  new_iupacs <- purrr::map_chr(modified_structures, structure_to_iupac.igraph)
  names(new_iupacs) <- names(modified_structures)
  attr(glycan, "iupacs") <- new_iupacs
  
  glycan
}

#' @export
remove_substituents.default <- function(glycan) {
  # Try to handle legacy glycan_structure objects (single igraphs with glycan_structure class)
  if (inherits(glycan, "glycan_structure") && inherits(glycan, "igraph")) {
    return(remove_substituents.igraph(glycan))
  }
  
  rlang::abort(c(
    "Cannot remove substituents for object of class {.cls {class(glycan)}}.",
    "i" = "Supported types: igraph object or glyrepr_structure vector."
  ))
}
