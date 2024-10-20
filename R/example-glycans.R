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
#'
#' @return A glycan graph object.
#' @export
n_glycan_core <- function(mode = "ne", linkage = TRUE) {
  if (!mode %in% c("ne", "dn")) rlang::abort("Mode must be 'ne' or 'dn'.")
  switch(
    mode,
    ne = n_glycan_core_ne(linkage),
    dn = n_glycan_core_dn(linkage),
  )
}


#' @rdname n_glycan_core
#' @export
n_glycan_core_ne <- function(linkage = TRUE) {
  if (!is.logical(linkage) && length(linkage) != 1) {
    rlang::abort("Linkage must be a single logical.")
  }
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  linkages <- c("b1-4", "b1-4", "a1-3", "a1-6")
  igraph::E(graph)$linkage <- if (linkage) linkages else NA_character_
  new_ne_glycan_graph(graph)
}


#' @rdname n_glycan_core
#' @export
n_glycan_core_dn <- function(linkage = TRUE) {
  if (!is.logical(linkage) && length(linkage) != 1) {
    rlang::abort("Linkage must be a single logical.")
  }
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 4-+5, 5-+6, 6-+7, 5-+8, 8-+9)
  igraph::V(graph)$type <- c("mono", rep(c("linkage", "mono"), 4))
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc", NA, "Man", NA, "Man", NA, "Man")
  linkages <- c(NA, "b1-4", NA, "b1-4", NA, "a1-3", NA, "a1-6", NA)
  igraph::V(graph)$linkage <- if (linkage) linkages else NA_character_
  new_dn_glycan_graph(graph)
}
