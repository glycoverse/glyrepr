#' Get Alditol Status
#'
#' Determine whether the reducing-end residue of a glycan is an alditol.
#' Alditol status is stored as the graph-level logical attribute `alditol`.
#' Legacy graphs without this attribute are treated as non-alditols.
#'
#' @param x A glycan structure vector or one glycan `igraph`.
#'
#' @returns A logical vector for structure-vector input, or one logical scalar
#'   for graph input. Missing structure-vector elements return `NA`.
#'
#' @examples
#' glycans <- as_glycan_structure(c(
#'   reduced = "Gal(b1-4)GlcNAc-ol(a1-",
#'   ordinary = "Gal(b1-4)GlcNAc(a1-"
#' ))
#' get_alditol(glycans)
#'
#' @export
get_alditol <- function(x) {
  if (inherits(x, "igraph")) {
    return(isTRUE(igraph::graph_attr(x, "alditol")))
  }

  checkmate::assert_class(x, "glyrepr_structure")
  smap_lgl(x, function(graph) {
    isTRUE(igraph::graph_attr(graph, "alditol"))
  })
}


normalize_alditol_attr <- function(graph) {
  if (is.null(igraph::graph_attr(graph, "alditol"))) {
    graph <- igraph::set_graph_attr(graph, "alditol", value = FALSE)
  }
  graph
}


format_reducing_end_iupac <- function(sequence, graph) {
  suffix <- if (isTRUE(igraph::graph_attr(graph, "alditol"))) {
    "-ol"
  } else {
    ""
  }
  paste0(sequence, suffix, "(", graph$anomer, "-")
}
