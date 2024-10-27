#' Generate Possible Linkages
#'
#' @description
#' Given an obscure linkage format (having "?", e.g. "a2-?"),
#' this function generates all possible linkages based on the format.
#'
#' A linkage character string has the following format: `xy-z`.
#' - `x` is the anomer, either "a" or "b".
#' - `y` is the first position, an integer from 1 to 2.
#' - `z` is the second position, an integer from 1 to 9.
#'
#' The ranges of possible anomers, first positions, and second positions
#' can be specified using `anomer_range`, `pos1_range`, and `pos2_range`.
#'
#' @param linkage A linkage string.
#' @param ... Not used.
#' @param anomer_range A character vector of possible anomers.
#' Default is `c("a", "b")`.
#' @param pos1_range A numeric vector of possible first positions.
#' Default is `1:2`.
#' @param pos2_range A numeric vector of possible second positions.
#' Default is `1:9`.
#'
#' @returns A character vector of possible linkages.
#'
#' @examples
#' possible_linkages("a2-?")
#' possible_linkages("??-2")
#' possible_linkages("a1-3")
#' possible_linkages("a?-?", pos1_range = 2, pos2_range = c(2, 3))
#'
#' @export
possible_linkages <- function(
  linkage,
  ...,
  anomer_range = c("a", "b"),
  pos1_range = 1:2,
  pos2_range = 1:9
) {
  # Input checks
  if (!is.character(linkage)) {
    rlang::abort("Linkage must be a character.")
  }
  if (length(linkage) != 1) {
    rlang::abort("Linkage must be a single character.")
  }
  if (!is.character(anomer_range)) {
    rlang::abort("Anomer range must be a character vector.")
  }
  if (!is.numeric(pos1_range) || !is.numeric(pos2_range)) {
    rlang::abort("Position ranges must be a numeric vector.")
  }

  # Check if the linkage is valid
  if (!valid_linkages(linkage)) {
    rlang::abort("Invalid linkage format.")
  }

  # Possible linkage elements
  current_anomer <- stringr::str_sub(linkage, 1, 1)
  current_pos1 <- stringr::str_sub(linkage, 2, 2)
  current_pos2 <- stringr::str_sub(linkage, 4, 4)

  anomers <- if (current_anomer == "?") anomer_range else current_anomer
  first_positions <- if (current_pos1 == "?") pos1_range else current_pos1
  second_positions <- if (current_pos2 == "?") pos2_range else current_pos2

  # Generate possible linkages
  purrr::pmap_chr(
    expand.grid(anomers, first_positions, second_positions),
    ~ paste0(..1, ..2, "-", ..3)
  )
}
