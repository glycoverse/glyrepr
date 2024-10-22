#' N-glycan Core
#'
#' @description
#' Create a N-glycan core graph, either in NE or DN mode.
#'
#' ```
#'     Man
#'   a1-6 \   b1-4      b1-4
#'         Man -- GlcNAc -- GlcNAc
#'   a1-3 /
#'     Man
#' ```
#'
#' @param mode A character string, either "ne" or "dn".
#' @param linkage A logical indicating whether to include linkages.
#'   Default is `TRUE`.
#' @param mono_type A character string specifying the type of monosaccharides.
#'   Can be "simple" (H, N, F, S), "generic" (Hex, HexNAc, dHex, NeuAc, etc.),
#'   or "concrete" (Man, Gal, HexNAc, Fuc, etc.). Default is "concrete".
#'
#' @return A glycan graph object.
#'
#' @examples
#' print(n_glycan_core(mode = "ne"), verbose = TRUE)
#' print(n_glycan_core(mode = "dn"), verbose = TRUE)
#'
#' @export
n_glycan_core <- function(mode = "ne", linkage = TRUE, mono_type = "concrete") {
  if (!mode %in% c("ne", "dn")) {
    rlang::abort("Mode must be 'ne' or 'dn'.")
  }
  if (!mono_type %in% c("simple", "generic", "concrete")) {
    rlang::abort("Mono type must be 'simple', 'generic', or 'concrete'.")
  }
  if (!is.logical(linkage) && length(linkage) != 1) {
    rlang::abort("Linkage must be a single logical.")
  }
  glycan <- n_glycan_core_ne()
  if (!linkage) {
    igraph::E(glycan)$linkage <- NA_character_
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
  igraph::E(graph)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  new_ne_glycan_graph(graph)
}
