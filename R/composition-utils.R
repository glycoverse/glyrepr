#' Get the Number of Monosaccharides
#'
#' Get the number of monosaccharides in a glycan composition or glycan structure.
#' When `mono` is "generic" (e.g. "Hex", "HexNAc"), 
#' it counts all "concrete" monosaccharides that match.
#' For example, "Hex" will count all Glc, Man, Gal, etc.
#' When `mono` is "concrete" (e.g. "Gal", "GalNAc"),
#' NA is returned when the composition is "generic".
#'
#' @param x A glycan composition (`glyrepr_composition`) or
#'   a glycan structure (`glyrepr_structure`) vector
#' @param mono The monosaccharide to count. A character scalar.
#'
#' @returns A numeric vector of the same length as `x`.
#'
#' @examples
#' comp <- glycan_composition(c(Hex = 5, HexNAc = 2), c(Gal = 1, Man = 1,GalNAc = 1))
#' count_mono(comp, "Hex")
#' count_mono(comp, "Gal")
#'
#' struct <- glycan_structure("Gal(b1-3)GlcNAc(b1-4)Glc(a1-)")
#' count_mono(struct, "Gal")
#'
#' @export
count_mono <- function(x, mono) {
  UseMethod("count_mono")
}

.check_mono_arg <- function(mono) {
  checkmate::assert_string(mono)
  if (!is_known_mono(mono)) {
    cli::cli_abort("{.arg mono} must be a known monosaccharide.")
  }
}

#' @rdname count_mono
#' @export
count_mono.glyrepr_composition <- function(x, mono) {
  .check_mono_arg(mono)
  mono_type <- get_mono_type(mono)
  if (mono_type == "generic") {
    x <- convert_mono_type(x, "generic")
  }
  data <- vctrs::field(vctrs::vec_data(x), "data")
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
count_mono.glyrepr_structure <- function(x, mono) {
  .check_mono_arg(mono)
  comps <- as_glycan_composition(x)
  count_mono(comps, mono)
}
