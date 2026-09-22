#' Validate and Canonicalize a Batch of Glycan Graphs
#'
#' Return canonical graphs and their IUPAC keys together. For ordinary trees,
#' ordering and key generation share one traversal. Graph, vertex, and edge
#' attributes are retained with their corresponding objects. Unlike structure
#' vector construction, this function does not deduplicate graphs: equal
#' structures can retain different source attributes.
#'
#' Strict failures have class `glyrepr_error_structure_failure` with
#' `position`, `input_name`, and `reason` fields. Recovery warnings have class
#' `glyrepr_warning_structure_failure` with `positions` and `reasons` fields.
#'
#' @param graphs A list of glycan `igraph` objects. `NULL` elements represent
#'   missing structures. List names and positions are preserved.
#' @param validate Whether to validate each graph before canonicalization.
#'   Set to `FALSE` only for graphs already validated with
#'   [validate_glycan_graph()]. This skips semantic validation, not
#'   canonicalization. Array records should use [structure_from_arrays()].
#' @param on_failure Either `"error"` (default) or `"na"`. With `"na"`, invalid
#'   graphs produce a warning and missing output, with details in `reason`.
#' @returns A list containing aligned, named vectors `iupac`, `status`, and
#'   `reason`, and an aligned named list `graphs`. Status is `"ok"`, `"missing"`,
#'   or `"invalid"`. Missing and invalid entries have a `NULL` graph and `NA`
#'   key. Reasons are `NA` except for invalid entries.
#' @examples
#' graphs <- as.list(n_glycan_core())
#' result <- canonicalize_glycan_graphs(graphs)
#' result$iupac
#' new_glycan_structure(result$iupac, stats::setNames(result$graphs, result$iupac))
#' @export
canonicalize_glycan_graphs <- function(
  graphs,
  validate = TRUE,
  on_failure = c("error", "na")
) {
  checkmate::assert_list(graphs)
  checkmate::assert_flag(validate)
  on_failure <- match.arg(on_failure)
  outcomes <- lapply(graphs, function(graph) {
    if (is.null(graph)) {
      return(NULL)
    }
    tryCatch(
      {
        checkmate::assert_class(graph, "igraph")
        if (validate) {
          graph <- validate_glycan_graph(graph)
        }
        canonicalize_graph_with_iupac(graph)
      },
      error = identity
    )
  })
  failed <- vapply(outcomes, inherits, logical(1), "error")
  reasons <- rep(NA_character_, length(graphs))
  reasons[failed] <- vapply(
    outcomes[failed],
    normalize_structure_failure_reason,
    character(1)
  )
  if (any(failed)) {
    if (on_failure == "error") {
      i <- which(failed)[[1]]
      .abort_structure_failure(
        outcomes[[i]],
        i,
        names(graphs),
        rlang::current_env()
      )
    }
    warn_structure_failures(which(failed), reasons[failed], names(graphs))
  }
  result <- list(
    graphs = lapply(outcomes, function(x) {
      if (is.null(x) || inherits(x, "error")) NULL else x$graph
    }),
    iupac = vapply(
      outcomes,
      function(x) {
        if (is.null(x) || inherits(x, "error")) NA_character_ else x$iupac
      },
      character(1)
    ),
    status = as.character(ifelse(
      failed,
      "invalid",
      ifelse(vapply(outcomes, is.null, logical(1)), "missing", "ok")
    )),
    reason = reasons
  )
  lapply(result, stats::setNames, names(graphs))
}
