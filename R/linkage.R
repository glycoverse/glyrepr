#' Determine if a Glycan Structure has Linkages
#'
#' @description
#' Unknown linkages in a glycan structure are represented by "??-?".
#' This function checks if all linkages in a glycan structure are unknown.
#' Note that even only one linkage is partial known (e.g. "a?-?"),
#' this function will return `TRUE`.
#'
#' @param glycan A glyrepr_structure vector.
#'
#' @return A logical vector indicating if each glycan structure has linkages.
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
  if (!is_glycan_structure(glycan)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }
  
  smap_lgl(glycan, .has_linkages_single)
}

# Internal function to check linkages in a single igraph
.has_linkages_single <- function(glycan) {
  any(igraph::E(glycan)$linkage != "??-?")
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
    cli::cli_abort("Invalid linkage format.")
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
#' @param glycan A glyrepr_structure vector.
#'
#' @return A glyrepr_structure vector with all linkages removed.
#'
#' @examples
#' glycan <- o_glycan_core_1(linkage = TRUE)
#' glycan
#' remove_linkages(glycan)
#'
#' @export
remove_linkages <- function(glycan) {
  if (!is_glycan_structure(glycan)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }
  
  smap_structure(glycan, .remove_linkages_single)
}

# Internal function to remove linkages from a single igraph
.remove_linkages_single <- function(glycan) {
  igraph::set_edge_attr(glycan, "linkage", value = "??-?")
}
