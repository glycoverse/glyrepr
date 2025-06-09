#' Determine if a Glycan Structure has Linkages
#'
#' @description
#' Unknown linkages in a glycan structure are represented by "??-?".
#' This function checks if all linkages in a glycan structure are unknown.
#' Note that even only one linkage is partial known (e.g. "a?-?"),
#' this function will return `TRUE`.
#'
#' @param glycan An igraph object representing a glycan structure, or a glyrepr_structure vector.
#'
#' @return A logical value indicating if the glycan structure has linkages, or a vector of logical values if input is vectorized.
#'
#' @examples
#' glycan <- o_glycan_core_1(linkage = TRUE)
#' has_linkages(glycan)
#' print(glycan)
#'
#' glycan <- remove_linkages(glycan)
#' has_linkages(glycan)
#' print(glycan)
#'
#' @seealso [remove_linkages()], [possible_linkages()]
#'
#' @export
has_linkages <- function(glycan) {
  UseMethod("has_linkages")
}

#' @export
has_linkages.igraph <- function(glycan) {
  any(igraph::E(glycan)$linkage != "??-?")
}

#' @export
has_linkages.glyrepr_structure <- function(glycan) {
  # Extract individual igraph objects and apply has_linkages to each
  data <- vctrs::vec_data(glycan)
  codes <- vctrs::field(data, "codes")
  structures <- attr(glycan, "structures")
  
  # Get individual structures and apply has_linkages
  individual_graphs <- purrr::map(codes, ~ structures[[.x]])
  purrr::map_lgl(individual_graphs, has_linkages.igraph)
}

#' @export
has_linkages.default <- function(glycan) {
  # Try to handle legacy glycan_structure objects (single igraphs with glycan_structure class)
  if (inherits(glycan, "glycan_structure") && inherits(glycan, "igraph")) {
    return(has_linkages.igraph(glycan))
  }
  
  rlang::abort(c(
    "Cannot check linkages for object of class {.cls {class(glycan)}}.",
    "i" = "Supported types: igraph object or glyrepr_structure vector."
  ))
}


#' Generate Possible Linkages
#'
#' @description
#' Given an obscure linkage format (having "?", e.g. "a2-?"),
#' this function generates all possible linkages based on the format.
#' See [valid_linkages()] for details.
#'
#' The ranges of possible anomers, first positions, and second positions
#' can be specified using `anomer_range`, `pos1_range`, and `pos2_range`.
#'
#' @param linkage A linkage string.
#' @param anomer_range A character vector of possible anomers.
#' Default is `c("a", "b")`.
#' @param pos1_range A numeric vector of possible first positions.
#' Default is `1:2`.
#' @param pos2_range A numeric vector of possible second positions.
#' Default is `1:9`.
#' @param include_unknown A logical value. If `TRUE`, "?" will be included.
#' Default is `FALSE`.
#'
#' @returns A character vector of possible linkages.
#'
#' @examples
#' possible_linkages("a2-?")
#' possible_linkages("??-2")
#' possible_linkages("a1-3")
#' possible_linkages("a?-?", pos1_range = 2, pos2_range = c(2, 3))
#' possible_linkages("?1-6", include_unknown = TRUE)
#'
#' @seealso [has_linkages()], [remove_linkages()], [valid_linkages()]
#'
#' @export
possible_linkages <- function(
  linkage,
  anomer_range = c("a", "b"),
  pos1_range = 1:2,
  pos2_range = 1:9,
  include_unknown = FALSE
) {
  # Input checks
  checkmate::assert_character(linkage, len = 1)
  checkmate::assert_character(anomer_range, pattern = "^[ab]$", unique = TRUE)
  checkmate::assert_numeric(pos1_range, lower = 1, upper = 2, unique = TRUE)
  checkmate::assert_numeric(pos2_range, lower = 1, upper = 9, unique = TRUE)
  checkmate::assert_flag(include_unknown)

  # Check if the linkage is valid
  if (!valid_linkages(linkage)) {
    rlang::abort("Invalid linkage format.")
  }

  # Add unknown elements
  if (include_unknown) {
    anomer_range <- c(anomer_range, "?")
    pos1_range <- c(pos1_range, "?")
    pos2_range <- c(pos2_range, "?")
  }

  # Possible linkage elements
  current_anomer <- stringr::str_sub(linkage, 1, 1)
  current_pos1 <- stringr::str_sub(linkage, 2, 2)
  current_pos2 <- stringr::str_sub(linkage, 4, -1)

  anomers <- if (current_anomer == "?") anomer_range else current_anomer
  first_positions <- if (current_pos1 == "?") pos1_range else current_pos1
  if (current_pos2 == "?") {
    second_positions <- pos2_range
  } else if (stringr::str_detect(current_pos2, "/")) {
    second_positions <- stringr::str_split(current_pos2, "/")[[1]]
  } else {
    second_positions <- current_pos2
  }

  # Generate possible linkages
  purrr::pmap_chr(
    expand.grid(anomers, first_positions, second_positions),
    ~ paste0(..1, ..2, "-", ..3)
  )
}


#' Remove All Linkages from a Glycan
#'
#' This function replaces all linkages in a glycan structure with "??-?".
#'
#' @param glycan An igraph object representing a glycan structure, or a glyrepr_structure vector.
#'
#' @return An igraph object representing a glycan structure, or a glyrepr_structure vector.
#'
#' @examples
#' glycan <- o_glycan_core_1(linkage = TRUE)
#' glycan
#' remove_linkages(glycan)
#'
#' @export
remove_linkages <- function(glycan) {
  UseMethod("remove_linkages")
}

#' @export
remove_linkages.igraph <- function(glycan) {
  igraph::set_edge_attr(glycan, "linkage", value = "??-?")
}

#' @export
remove_linkages.glyrepr_structure <- function(glycan) {
  # Extract individual igraph objects and apply remove_linkages to each
  data <- vctrs::vec_data(glycan)
  codes <- vctrs::field(data, "codes")
  structures <- attr(glycan, "structures")
  
  # Remove linkages from each unique structure
  modified_structures <- purrr::map(structures, ~ remove_linkages.igraph(.x))
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
remove_linkages.default <- function(glycan) {
  # Try to handle legacy glycan_structure objects (single igraphs with glycan_structure class)
  if (inherits(glycan, "glycan_structure") && inherits(glycan, "igraph")) {
    return(remove_linkages.igraph(glycan))
  }
  
  rlang::abort(c(
    "Cannot remove linkages for object of class {.cls {class(glycan)}}.",
    "i" = "Supported types: igraph object or glyrepr_structure vector."
  ))
}
