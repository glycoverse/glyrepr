#' Fill Anomer Positions
#'
#' Add anomer positions to glycan structures with missing anomer position
#' information. For example, `"Gal(??-?)GalNAc(??-"` is converted to
#' `"Gal(?1-?)GalNAc(?1-"`.
#'
#' For anomer positions that are already specified in the input structures,
#' this function does not modify them.
#'
#' @param strucs A [glycan_structure()] vector with concrete monosaccharides.
#'
#' @returns A [glycan_structure()] vector with anomer positions added where
#'   missing.
#'
#' @examples
#' glycans <- as_glycan_structure(c(
#'   "Gal(??-?)GalNAc(??-",
#'   "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
#' ))
#' fill_anomer_pos(glycans)
#'
#' @export
fill_anomer_pos <- function(strucs) {
  checkmate::assert_class(strucs, "glyrepr_structure")
  .assert_concrete_structure(strucs)
  smap_structure(strucs, .fill_anomer_pos_single)
}


#' Assert That a Glycan Structure Has Concrete Monosaccharides
#'
#' @param strucs A [glycan_structure()] vector.
#'
#' @returns `strucs`, invisibly. Throws an error if non-missing structures are
#'   not concrete.
#' @noRd
.assert_concrete_structure <- function(strucs) {
  mono_type <- get_mono_type(strucs)
  if (length(mono_type) == 0 || is.na(mono_type)) {
    return(invisible(strucs))
  }

  if (mono_type != "concrete") {
    cli::cli_abort("{.arg strucs} must have concrete monosaccharides.")
  }

  invisible(strucs)
}


#' Fill One Glycan Graph's Anomer Positions
#'
#' @param struc An igraph glycan structure.
#'
#' @returns An igraph glycan structure with missing anomer positions filled.
#' @noRd
.fill_anomer_pos_single <- function(struc) {
  root <- which(igraph::degree(struc, mode = "in") == 0)
  root_mono <- igraph::vertex_attr(struc, "mono", index = root)
  root_anomer <- igraph::graph_attr(struc, "anomer")
  struc <- igraph::set_graph_attr(
    struc,
    "anomer",
    value = .fill_anomer_pos_value(root_anomer, root_mono)
  )

  linkages <- igraph::edge_attr(struc, "linkage")
  if (length(linkages) == 0) {
    return(struc)
  }

  edges <- igraph::ends(struc, igraph::E(struc), names = FALSE)
  donor_monos <- igraph::vertex_attr(struc, "mono", index = edges[, 2])
  linkages <- purrr::map2_chr(linkages, donor_monos, .fill_anomer_pos_value)
  igraph::set_edge_attr(struc, "linkage", value = linkages)
}


#' Fill the Anomer Position in One Linkage or Reducing-End Anomer
#'
#' @param anomer A linkage string like `"??-?"` or reducing-end anomer like
#'   `"??"`.
#' @param mono A concrete monosaccharide name.
#'
#' @returns `anomer` with the second character filled when missing.
#' @noRd
.fill_anomer_pos_value <- function(anomer, mono) {
  if (stringr::str_sub(anomer, 2, 2) != "?") {
    return(anomer)
  }

  stringr::str_c(
    stringr::str_sub(anomer, 1, 1),
    infer_anomer_pos(mono),
    stringr::str_sub(anomer, 3)
  )
}
