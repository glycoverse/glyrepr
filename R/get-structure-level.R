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
#' structures <- as_glycan_structure(c(
#'   "Gal(b1-3)GalNAc(a1-",
#'   "Gal(b1-?)GalNAc(a1-",
#'   "Gal(??-?)GalNAc(??-",
#'   "Hex(??-?)HexNAc(??-",
#'   "Hex(b1-3)HexNAc(a1-"
#' ))
#' get_structure_level(structures)
#'
#' @export
get_structure_level <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")
  data <- vctrs::vec_data(x)
  mono_types <- vctrs::field(data, "mono_type")
  has_linkages_strict <- has_linkages(x, strict = TRUE)
  has_linkages_lenient <- has_linkages(x, strict = FALSE)
  dplyr::case_when(
    mono_types == "concrete" & has_linkages_strict ~ "intact",
    mono_types == "concrete" & (!has_linkages_strict) & has_linkages_lenient ~ "partial",
    mono_types == "concrete" & (!has_linkages_lenient) ~ "topological",
    mono_types == "generic" & (!has_linkages_strict) ~ "basic",
    .default = "basic"
  )
}
