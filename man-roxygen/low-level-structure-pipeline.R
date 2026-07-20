#' @section Name-preserving manual construction:
#' The five low-level functions can reproduce strict graph-based construction
#' while preserving the names of the input graph list:
#'
#' ```r
#' input_names <- names(graphs)
#' graphs <- unname(graphs)
#'
#' graphs <- purrr::map(graphs, validate_glycan_graph)
#' graphs <- purrr::map(graphs, canonicalize_glycan_graph)
#' validate_glycan_graph_vector(graphs)
#'
#' iupacs <- purrr::map_chr(graphs, graph_to_iupac)
#' names(iupacs) <- input_names
#'
#' unique <- !duplicated(unname(iupacs))
#' unique_graphs <- graphs[unique]
#' names(unique_graphs) <- unname(iupacs[unique])
#'
#' new_glycan_structure(iupacs, unique_graphs)
#' ```
#'
#' Unlike `as_glycan_structure(graphs, on_failure = "na")`, this strict
#' pipeline stops at the first invalid graph.
