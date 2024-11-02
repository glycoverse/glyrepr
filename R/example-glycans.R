#' Example Glycan Graphs
#'
#' @description
#' Create example glycan graphs for testing and demonstration.
#' Includes **N-glycan core** and **O-glycan core 1** and **core 2**.
#'
#' @details
#' # N-Glycan Core
#'
#' **N-Glycans** are branched oligosaccharides that are bound, most commonly,
#' via GlcNAc to an Asn residue of the protein backbone.
#' A common motif of all N-glycans is the **chitobiose core**,
#' composed of three mannose and two GlcNAc moieties,
#' which is commonly attached to the protein backbone via GlcNAc.
#' The mannose residue is branched and connected via α1,3- and α1,6-glycosidic
#' linkages to the two other mannose building blocks.
#'
#' ```
#'     Man
#'   a1-6 \   b1-4      b1-4
#'         Man -- GlcNAc -- GlcNAc
#'   a1-3 /
#'     Man
#' ```
#'
#' # O-Glycan Core
#'
#' **O-Glycans** are highly abundant in extracellular proteins.
#' Generally, O-glycans are extended following four major core structures:
#' **core 1**, **core 2**, core 3, and core 4.
#' The first two are by far the most common core structures in O-glycosylation
#' and are found throughout the body.
#'
#' **core 1**:
#' ```
#'     GalNAc
#'    / b1-3
#' Gal
#' ```
#'
#' **core 2**:
#' ```
#' GlcNAc
#'       \ b1-6
#'        GalNAc
#'       / b1-3
#'    Gal
#' ```
#'
#' @param mode A character string, either "ne" or "dn".
#' See [as_glycan_graph()] for details.
#' @param linkage A logical indicating whether to include linkages (e.g. "b1-4").
#'   Default is `TRUE`.
#' @param mono_type A character string specifying the type of monosaccharides.
#'   Can be "simple" (H, N, F, S), "generic" (Hex, HexNAc, dHex, NeuAc, etc.),
#'   or "concrete" (Man, Gal, HexNAc, Fuc, etc.). Default is "concrete".
#'
#' @return A glycan graph object.
#'
#' @examples
#' print(n_glycan_core(mode = "ne"), verbose = TRUE)
#' print(o_glycan_core_1(mode = "dn"), verbose = TRUE)
#'
#' @export
n_glycan_core <- function(mode = "ne", linkage = TRUE, mono_type = "concrete") {
  build_example_graph(mode, linkage, mono_type, n_glycan_core_ne)
}


#' @rdname n_glycan_core
#' @export
o_glycan_core_1 <- function(mode = "ne", linkage = TRUE, mono_type = "concrete") {
  build_example_graph(mode, linkage, mono_type, o_glycan_core_1_ne)
}


#' @rdname n_glycan_core
#' @export
o_glycan_core_2 <- function(mode = "ne", linkage = TRUE, mono_type = "concrete") {
  build_example_graph(mode, linkage, mono_type, o_glycan_core_2_ne)
}


build_example_graph <- function(mode, linkage, mono_type, builder) {
  if (!mode %in% c("ne", "dn")) {
    rlang::abort("Mode must be 'ne' or 'dn'.")
  }
  if (!mono_type %in% c("simple", "generic", "concrete")) {
    rlang::abort("Mono type must be 'simple', 'generic', or 'concrete'.")
  }
  if (!is.logical(linkage) && length(linkage) != 1) {
    rlang::abort("Linkage must be a single logical.")
  }
  glycan <- builder()
  if (!linkage) {
    igraph::E(glycan)$linkage <- "??-?"
  }
  if (mode == "dn") {
    glycan <- convert_ne_to_dn(glycan)
  }
  if (mono_type %in% c("simple", "generic")) {
    glycan <- convert_glycan_mono_type(glycan, mono_type)
  }
  glycan
}


n_glycan_core_ne <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  graph$anomer <- "?1"
  graph$alditol <- FALSE
  new_ne_glycan_graph(graph)
}


o_glycan_core_1_ne <- function() {
  graph <- igraph::make_graph(~ 1-+2)
  igraph::V(graph)$mono <- c("GalNAc", "Gal")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-3"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  new_ne_glycan_graph(graph)
}


o_glycan_core_2_ne <- function() {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("GalNAc", "Gal", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-3", "b1-6")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  new_ne_glycan_graph(graph)
}
