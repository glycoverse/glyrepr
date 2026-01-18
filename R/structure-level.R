#' Get the Structure Resolution Levels
#'
#' @description
#' Glycan structures can have four possible levels of resolution:
#' - "intact": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   and all linkages are fully determined (e.g. "a2-3", "b1-4").
#' - "partial": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   but some linkage information is missing (e.g. "a2-?").
#' - "topological": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   but the linkage information is completely unknown ("??-?").
#' - "basic": All monosaccharides are generic (e.g. "Hex", "HexNAc"),
#'   and the linkage information is completely unknown ("??-?").
#'
#' Note that in theory you can have a glycan with generic monosaccharides with all linkages determined.
#' For example, "Hex(b1-3)HexNAc(a1-" is a valid glycan structure.
#' But in reality, this is almost impossible,
#' because linkage information is far more difficult to acquire than monosaccharide information.
#' This kind of glycan structure is also assigned to "basic" level.
#'
#' @param x A [glycan_structure()] vector.
#'
#' @returns A character vector of the same length as `x`,
#'   containing the structure level for each element.
#'
#' @examples
#' glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
#' get_structure_level(glycan)
#'
#' @seealso [has_linkages()], [get_mono_type()]
#' @export
get_structure_level <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")

  # Capture input names for preservation
  input_names <- names(x)

  codes <- vctrs::vec_data(x)
  mono_type <- get_mono_type.glyrepr_structure(x)
  has_linkages_strict <- has_linkages(x, strict = TRUE)
  has_linkages_lenient <- has_linkages(x, strict = FALSE)

  result <- dplyr::case_when(
    mono_type == "concrete" & has_linkages_strict ~ "intact",
    mono_type == "concrete" & (!has_linkages_strict) & has_linkages_lenient ~ "partial",
    mono_type == "concrete" & (!has_linkages_lenient) ~ "topological",
    mono_type == "generic" & (!has_linkages_strict) ~ "basic",
    .default = "basic"
  )

  # Restore names
  names(result) <- input_names

  result
}

#' Reduce a Glycan Structure to a Lower Resolution Level
#'
#' This function reduces a glycan structure from a higher resolution level to a lower resolution level
#' (see [get_structure_level()] for four possible levels of resolution).
#' For example, it can reduce an "intact" structure to a "topological" structure,
#' or a "partial" structure to a "basic" structure.
#' One exception is that you can never reduce an "intact" structure to "partial" level,
#' because the "partial" level is not deterministic.
#'
#' @details
#' The logic is as follows:
#' - If `to_level` is "topological", this function calls [remove_linkages()] to remove all linkages.
#' - If `to_level` is "basic", this function calls [remove_linkages()] to remove all linkages,
#'   and [convert_to_generic()] to convert all monosaccharides to generic.
#'
#' @param x A [glycan_structure()] vector.
#' @param to_level The resolution level to reduce to. Can be "basic" or "topological".
#'   Must be a lower resolution level than any structure in `x`
#'   ("intact" > "partial" > "topological" > "basic").
#'   If `to_level` is the same as some structure in `x`, the result will be the same as the input.
#'   You can use [get_structure_level()] to check the structure levels of `x`.
#'
#' @returns A [glycan_structure()] vector reduced to the given resolution level.
#' @examples
#' glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
#' reduce_structure_level(glycan, to_level = "topological")
#'
#' @seealso [get_structure_level()]
#' @export
reduce_structure_level <- function(x, to_level) {
  checkmate::assert_class(x, "glyrepr_structure")
  checkmate::assert_choice(to_level, c("basic", "topological"))

  # Check if the target level is lower than any structure in `x`
  struc_levels <- get_structure_level(x)
  level_ranks <- c("basic" = 1, "topological" = 2, "partial" = 3, "intact" = 4)
  from_level_ranks <- level_ranks[struc_levels]
  to_level_rank <- level_ranks[to_level]
  if (any(from_level_ranks < to_level_rank)) {
    larger_levels <- unique(struc_levels[from_level_ranks < to_level_rank])
    cli::cli_abort(c(
      "Cannot reduce a structure to a higher resolution level.",
      "x" = "Some structures in {.arg x} have levels: {.val {larger_levels}}.",
      "i" = "Target level: {.val {to_level}} (> {.val {larger_levels}}).",
      "i" = "You can use {.fn get_structure_level} to check the structure levels of {.arg x}."
    ))
  }

  # Reduce the structure level
  if (to_level == "basic") {
    x <- remove_linkages(x)
    x <- convert_to_generic(x)
  } else { # `to_level` is "topological"
    x <- remove_linkages(x)
  }
  x
}