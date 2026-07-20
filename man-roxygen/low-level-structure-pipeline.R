#' @section Low-level API warning:
#' These functions are low-level, developer-facing APIs. Calling them directly
#' is usually not a good idea unless you understand and can guarantee all glycan
#' graph and `glyrepr_structure` invariants. Prefer [as_glycan_structure()] for
#' ordinary construction. Incorrect use of these functions can create invalid
#' structure vectors that fail in later operations.
#'
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
