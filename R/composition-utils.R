#' Get the Number of Monosaccharides
#'
#' Get the number of monosaccharides in a glycan composition or glycan structure.
#' When `mono` is "generic" (e.g. "Hex", "HexNAc"), 
#' it counts all "concrete" monosaccharides that match.
#' For example, "Hex" will count all Glc, Man, Gal, etc.
#'
#' @param x A glycan composition (`glyrepr_composition`) or
#'   a glycan structure (`glyrepr_structure`) vector
#' @param mono The monosaccharide to count. A character scalar.
#'
#' @returns A numeric vector of the same length as `x`.
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
  if (mono %in% monosaccharides$generic) {
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
  purrr::map_int(data, count_one, mono = mono)
}

#' @rdname count_mono
#' @export
count_mono.glyrepr_structure <- function(x, mono) {
  .check_mono_arg(mono)
  comps <- as_glycan_composition(x)
  count_mono(comps, mono)
}
