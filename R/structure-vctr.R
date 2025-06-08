#' Create a Structure Vector
#'
#' A vector of glycan structures.
#'
#' @details
#' The underlying implementation uses hash values of IUPAC codes of the glycan structures.
#' This prevents redundant storage and computation,
#' which is very useful for glycan structures.
#'
#' @param x A list of glycan structures.
#'
#' @return A structure_vctr object.
#' @export
structure_vctr <- function(x = list()) {
  # Check arguments
  checkmate::assert_list(x)

  # Use hash values of IUPAC codes of the glycan structures.
  iupacs <- purrr::map_chr(x, structure_to_iupac)
  codes <- purrr::map_chr(iupacs, rlang::hash)

  # Create a unique list of x based on uniqueness of codes.
  unique_indices <- which(!duplicated(codes))
  x <- x[unique_indices]
  names(x) <- codes[unique_indices]
  iupacs <- iupacs[unique_indices]
  names(iupacs) <- codes[unique_indices]

  new_structure_vctr(codes, iupacs, x)
}

new_structure_vctr <- function(codes = character(), iupacs = character(), structures = list()) {
  vctrs::new_vctr(
    codes,
    iupacs = iupacs,
    structures = structures,
    class = "glyrepr_structure_vctr"
  )
}

#' @export
#' @rdname structure_vctr
is_structure_vctr <- function(x) {
  inherits(x, "glyrepr_structure_vctr")
}

#' @export
format.glyrepr_structure_vctr <- function(x, ...) {
  unname(attr(x, "iupacs")[vctrs::vec_data(x)])
}

#' @export
vec_ptype_abbr.glyrepr_structure_vctr <- function(x, ...) "struc_vctr"

#' @export
vec_ptype_full.glyrepr_structure_vctr <- function(x, ...) "structure_vctr"

#' @export
obj_print_footer.glyrepr_structure_vctr <- function(x, ...) {
  cat("# Unique structures: ", format(length(attr(x, "iupacs"))), "\n", sep = "")
}
