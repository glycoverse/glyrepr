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
  switch(
    mode,
    ne = n_glycan_core_ne(linkage, mono_type),
    dn = n_glycan_core_dn(linkage, mono_type),
  )
}


n_glycan_core_ne <- function(linkage, mono_type) {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph)$mono <- switch(
    mono_type,
    simple = c("N", "N", "H", "H", "H"),
    generic = c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"),
    concrete = c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  )
  linkages <- c("b1-4", "b1-4", "a1-3", "a1-6")
  igraph::E(graph)$linkage <- if (linkage) linkages else NA_character_
  new_ne_glycan_graph(graph)
}


n_glycan_core_dn <- function(linkage, mono_type) {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 4-+5, 5-+6, 6-+7, 5-+8, 8-+9)
  igraph::V(graph)$type <- c("mono", rep(c("linkage", "mono"), 4))
  igraph::V(graph)$mono <- switch(
    mono_type,
    simple = c("N", NA, "N", NA, "H", NA, "H", NA, "H"),
    generic = c("HexNAc", NA, "HexNAc", NA, "Hex", NA, "Hex", NA, "Hex"),
    concrete = c("GlcNAc", NA, "GlcNAc", NA, "Man", NA, "Man", NA, "Man")
  )
  linkages <- c(NA, "b1-4", NA, "b1-4", NA, "a1-3", NA, "a1-6", NA)
  igraph::V(graph)$linkage <- if (linkage) linkages else NA_character_
  new_dn_glycan_graph(graph)
}
