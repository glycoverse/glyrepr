#' Fill Anomer Positions
#'
#' Add anomer positions to glycan structures with missing anomer position
#' information. For example, `"Gal(??-?)GalNAc(??-"` is converted to
#' `"Gal(?1-?)GalNAc(?1-"`.
#'
#' For anomer positions that are already specified in the input structures,
#' this function does not modify them.
#'
#' For a structure with floating parts, the reducing-end position is inferred
#' from the root of the main tree. Positions in virtual floating-part
#' attachment linkages are inferred from each floating part's root residue.
#'
#' @param strucs A [glycan_structure()] vector with concrete or generic
#'   monosaccharides.
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
  smap_structure(strucs, .fill_anomer_pos_single)
}


#' Fill One Glycan Graph's Anomer Positions
#'
#' @param struc An igraph glycan structure.
#'
#' @returns An igraph glycan structure with missing anomer positions filled.
#' @noRd
.fill_anomer_pos_single <- function(struc) {
  parts <- normalize_floating_parts(struc)
  main_vertices <- if (length(parts) == 0) {
    seq_len(igraph::vcount(struc))
  } else {
    floating_graph_info(struc, parts)$main_vertices
  }
  root <- intersect(
    which(igraph::degree(struc, mode = "in") == 0),
    main_vertices
  )
  root_mono <- igraph::vertex_attr(struc, "mono", index = root)
  root_anomer <- igraph::graph_attr(struc, "anomer")
  struc <- igraph::set_graph_attr(
    struc,
    "anomer",
    value = .fill_anomer_pos_value(root_anomer, root_mono)
  )

  linkages <- igraph::edge_attr(struc, "linkage")
  if (length(linkages) > 0) {
    edges <- igraph::ends(struc, igraph::E(struc), names = FALSE)
    donor_monos <- igraph::vertex_attr(struc, "mono", index = edges[, 2])
    linkages <- purrr::map2_chr(linkages, donor_monos, .fill_anomer_pos_value)
    struc <- igraph::set_edge_attr(struc, "linkage", value = linkages)
  }

  parts <- purrr::map(parts, function(part) {
    root_mono <- igraph::vertex_attr(struc, "mono", index = part$root)
    part$linkage <- .fill_anomer_pos_value(part$linkage, root_mono)
    part
  })
  set_floating_parts_attr(struc, parts)
}


#' Fill the Anomer Position in One Linkage or Reducing-End Anomer
#'
#' @param anomer A linkage string like `"??-?"` or reducing-end anomer like
#'   `"??"`.
#' @param mono A concrete or generic monosaccharide name.
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
