#' Available Substituents
#'
#' Get the available substituents for monosaccharides.
#'
#' @returns A character vector.
#'
#' @examples
#' available_substituents()
#'
#' @export
available_substituents <- function() {
  c("Me", "Ac", "NAc", "P", "S", "Pyr", "PC", "PPEtn", "PEtn", "N")
}

#' Normalize Substituent String
#'
#' Takes a substituent string (potentially with multiple substituents) and 
#' returns a normalized string with substituents sorted by position.
#'
#' @param sub A character string representing substituents, e.g., "4Ac,3Me" or "6S"
#' @returns A character string with substituents sorted by position, e.g., "3Me,4Ac"
#'
#' @examples
#' normalize_substituents("4Ac,3Me")  # Returns "3Me,4Ac"
#' normalize_substituents("6S")       # Returns "6S"
#' normalize_substituents("")         # Returns ""
#' @noRd
normalize_substituents <- function(sub) {
  checkmate::assert_character(sub, len = 1)
  
  if (sub == "") {
    return("")
  }
  
  # Split by commas
  individual_subs <- stringr::str_split(sub, ",")[[1]]
  
  # Remove any empty strings (in case of double commas)
  individual_subs <- individual_subs[individual_subs != ""]
  
  if (length(individual_subs) == 0) {
    return("")
  }
  
  # Extract positions for sorting
  positions <- purrr::map_chr(individual_subs, ~ stringr::str_extract(.x, "^[\\d\\?]"))
  
  # Convert to numeric for sorting (? becomes Inf)
  numeric_positions <- purrr::map_dbl(positions, function(pos) {
    if (pos == "?") Inf else as.numeric(pos)
  })
  
  # Sort and combine
  sorted_indices <- order(numeric_positions)
  sorted_subs <- individual_subs[sorted_indices]
  
  stringr::str_c(sorted_subs, collapse = ",")
}

#' Remove All Substituents from a Glycan
#'
#' This function replaces all substituents in a glycan structure with empty strings.
#'
#' @param glycan A glyrepr_structure vector.
#'
#' @returns A glyrepr_structure vector with all substituents removed.
#'
#' @examples
#' (glycan <- glycan_structure(o_glycan_core_1()))
#' remove_substituents(glycan)
#'
#' @export
remove_substituents <- function(glycan) {
  if (!is_glycan_structure(glycan)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }

  smap_structure(glycan, .remove_substituents_single)
}

# Internal function to remove substituents from a single igraph
.remove_substituents_single <- function(glycan) {
  igraph::set_vertex_attr(glycan, "sub", value = "")
}