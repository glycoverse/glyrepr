#' Get the Number of Monosaccharides
#'
#' When `mono` is:
#' - `NULL` (default), returns the total number of monosaccharides and substituents.
#' - A string, returns the number of the specified monosaccharide or substituent.
#'
#' @details
#' When `mono` is "generic" (e.g. "Hex", "HexNAc"),
#' it counts all "concrete" monosaccharides that match.
#' For example, "Hex" will count all Glc, Man, Gal, etc.
#' When `mono` is "concrete" (e.g. "Gal", "GalNAc"),
#' NA is returned when the composition is "generic".
#'
#' @param x A glycan composition (`glyrepr_composition`) or
#'   a glycan structure (`glyrepr_structure`) vector.
#' @param mono The monosaccharide or substituent to count. A character scalar.
#'   If `NULL` (default), return the total number of monosaccharides.
#' @param include_subs Whether to include substituents when `mono` is `NULL`.
#'   Default is `FALSE`.
#'
#' @returns A numeric vector of the same length as `x`.
#'
#' @examples
#' comp <- glycan_composition(c(Hex = 5, HexNAc = 2), c(Gal = 1, Man = 1, GalNAc = 1))
#' count_mono(comp, "Hex")
#' count_mono(comp, "Gal")
#'
#' struct <- as_glycan_structure("Gal(b1-3)GlcNAc(b1-4)Glc(a1-")
#' count_mono(struct, "Gal")
#'
#' # Total number of monosaccharides
#' count_mono(comp)
#'
#' @export
count_mono <- function(x, mono = NULL, include_subs = FALSE) {
  UseMethod("count_mono")
}

.check_count_mono_args <- function(mono, include_subs) {
  checkmate::assert_flag(include_subs)
  if (is.null(mono)) {
    return()
  }
  checkmate::assert_string(mono)
  if (!mono %in% c(available_monosaccharides(), available_substituents())) {
    cli::cli_abort("{.arg mono} must be a known monosaccharide or substituent.")
  }
}

#' @rdname count_mono
#' @export
count_mono.glyrepr_composition <- function(x, mono = NULL, include_subs = FALSE) {
  .check_count_mono_args(mono, include_subs)

  # Special behavior when `mono` is NULL: count all monosaccharides and substituents
  if (is.null(mono)) {
    data <- vctrs::field(vctrs::vec_data(x), "data")
    if (!include_subs) {
      data <- purrr::map(data, ~ .x[names(.x) %in% available_monosaccharides()])
    }
    return(purrr::map_int(data, sum))
  }

  # Get the type of the monosaccharide
  if (mono %in% available_substituents()) {
    mono_type <- "substituent"  # the value is not used anywhere, just for readability
  } else {
    mono_type <- get_mono_type(mono)
  }
  # special monosaccharides are those having the same name for both generic and concrete types
  if (is.na(mono_type)) {
    mono_type <- "special"  # the value is not used anywhere, just for readability
  }

  # Convert to generic if needed
  if (mono_type == "generic") {
    x <- convert_to_generic(x)
  }
  data <- vctrs::field(vctrs::vec_data(x), "data")

  # Count the number of monosaccharides or substituents
  count_one <- function(one_mono, mono) {
    if (mono %in% names(one_mono)) {
      n <- one_mono[[mono]]
      if (is.na(n)) {
        n <- 0L
      }
    } else {
      n <- 0L
    }
    n
  }
  res <- purrr::map_int(data, count_one, mono = mono)
  if (mono_type == "concrete") {
    x_mono_type <- get_mono_type(x)
    res[x_mono_type == "generic"] <- NA_integer_
  }
  res
}

#' @rdname count_mono
#' @export
count_mono.glyrepr_structure <- function(x, mono = NULL, include_subs = FALSE) {
  .check_count_mono_args(mono, include_subs)
  comps <- as_glycan_composition(x)
  count_mono.glyrepr_composition(comps, mono, include_subs)
}
