#' Get the Structure Resolution Levels
#'
#' @description
#' Glycan structures can have four possible levels of resolution:
#' - "intact": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   and no linkage or anomer contains "?".
#' - "partial": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   at least one linkage or anomer contains "?",
#'   and at least one linkage or anomer has a non-"?" annotation.
#' - "topological": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
#'   and all linkages and anomers are completely unknown ("??-?"/"??").
#' - "basic": All monosaccharides are generic (e.g. "Hex", "HexNAc").
#'
#' Note that in theory you can have a glycan with generic monosaccharides with all linkages determined.
#' For example, "Hex(b1-3)HexNAc(a1-" is a valid glycan structure.
#' But in reality, this is almost impossible,
#' because linkage information is far more difficult to acquire than monosaccharide information.
#' This kind of glycan structure is also assigned to "basic" level.
#'
#' @param x A [glycan_structure()] vector.
#'
#' @returns A character scalar containing the structure level for `x`.
#'   If `x` is empty or all structures in `x` are NA, returns NA_character_.
#'
#' @examples
#' glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
#' get_structure_level(glycan)
#'
#' @seealso [has_linkages()], [get_mono_type()]
#' @export
get_structure_level <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")

  non_na <- !structure_na_mask(x)

  if (!any(non_na)) {
    if (length(x) == 0) {
      return(NA_character_)
    }
    return(NA_character_)
  }

  x_valid <- x[non_na]
  mono_type <- get_mono_type.glyrepr_structure(x_valid)
  has_linkages_strict <- has_linkages(x_valid, strict = TRUE)
  has_linkages_lenient <- has_linkages(x_valid, strict = FALSE)

  if (mono_type == "generic") {
    .warn_generic_linkage_structure_level(has_linkages_lenient)
    return("basic")
  }

  if (all(has_linkages_strict)) {
    return("intact")
  }

  if (any(has_linkages_lenient)) {
    return("partial")
  }

  "topological"
}

#' Warn About Generic Structures With Linkage Annotation
#'
#' Generic structures are always treated as basic resolution,
#' even if they contain linkage or anomer annotations.
#'
#' @param has_linkages_lenient A logical vector returned by [has_linkages()]
#'   with `strict = FALSE`.
#' @returns Nothing. Called for its warning side effect.
#' @noRd
.warn_generic_linkage_structure_level <- function(has_linkages_lenient) {
  if (any(has_linkages_lenient)) {
    cli::cli_warn(c(
      "Generic glycan structures with linkage annotations are treated as {.val basic}.",
      "i" = "Linkage information is ignored when residues are generic."
    ), class = "glyrepr_warning_generic_structure_linkages")
  }
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
#'   Must be a lower resolution level than `x`
#'   ("intact" > "partial" > "topological" > "basic").
#'   If `to_level` is the same as the structure level of `x`,
#'   the result will be the same as the input.
#'   You can use [get_structure_level()] to check the structure level of `x`.
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

  from_level <- get_structure_level(x)

  if (is.na(from_level)) {
    # Two situations can lead to NA structure level:
    #. 1. x is empty.
    #. 2. All structures in x are NA.
    # In both cases, we can just return x without any modification.
    return(x)
  }

  level_ranks <- c("basic" = 1, "topological" = 2, "partial" = 3, "intact" = 4)
  if (level_ranks[[from_level]] < level_ranks[[to_level]]) {
    cli::cli_abort(c(
      "Cannot reduce a structure to a higher resolution level.",
      "x" = "Structure level of {.arg x}: {.val {from_level}}.",
      "i" = "Target level: {.val {to_level}} (> {.val {from_level}}).",
      "i" = "You can use {.fn get_structure_level} to check the structure level of {.arg x}."
    ))
  }

  x <- remove_linkages(x)
  if (to_level == "basic") {
    x <- convert_to_generic(x)
  }
  x
}
