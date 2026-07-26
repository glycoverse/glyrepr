#' Localize Floating Glycan Parts
#'
#' @description
#' `localize_floating_parts()` attaches selected floating parts to
#' caller-supplied parent nodes. The assignments are interpreted against the
#' canonical node and part identifiers returned by
#' [structure_floating_candidates()].
#'
#' Each selected parent must belong to the floating part's declared candidate
#' domain. For an unrestricted `{<floating>}` part, that domain contains every
#' node in the original main tree. Assignments must also be simultaneously
#' compatible with occupied and potential acceptor linkage positions.
#'
#' Selected virtual attachments become ordinary graph edges. Unassigned
#' floating parts remain floating. The resulting structures are canonicalized,
#' and candidate-parent indices for remaining parts are remapped to the new
#' canonical main-tree order.
#'
#' Missing values, vector positions, and names in `x` are preserved.
#'
#' @param x A glycan structure vector.
#' @param assignments A data frame with integer columns `glycan_id`, `part_id`,
#'   and `parent_node`. Each `(glycan_id, part_id)` pair may occur at most once.
#'
#' @returns A glycan structure vector with the same length and names as `x`.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
#' )
#' assignments <- tibble::tibble(
#'   glycan_id = 1L,
#'   part_id = 1L,
#'   parent_node = 1L
#' )
#' localize_floating_parts(glycan, assignments)
#'
#' @export
localize_floating_parts <- function(x, assignments) {
  checkmate::assert_class(x, "glyrepr_structure")
  assignments <- validate_floating_assignments(assignments, length(x))

  if (nrow(assignments) == 0) {
    return(x)
  }

  graphs <- as.list(x)
  for (glycan_id in unique(assignments$glycan_id)) {
    if (is.null(graphs[[glycan_id]])) {
      cli::cli_abort(c(
        "Cannot localize a missing glycan structure.",
        "x" = "Glycan {.val {glycan_id}} is missing."
      ))
    }

    graph_assignments <- assignments[
      assignments$glycan_id == glycan_id,
      ,
      drop = FALSE
    ]
    graphs[[glycan_id]] <- localize_floating_graph(
      graphs[[glycan_id]],
      graph_assignments,
      glycan_id
    )
  }

  out <- do.call(glycan_structure, unname(graphs))
  names(out) <- names(x)
  out
}


#' Validate a floating localization assignment table
#'
#' @param assignments A candidate assignment table.
#' @param input_size Length of the structure vector.
#' @returns A tibble containing normalized assignment columns.
#' @noRd
validate_floating_assignments <- function(assignments, input_size) {
  required_cols <- c("glycan_id", "part_id", "parent_node")
  assignments <- validate_structure_table(
    assignments,
    required_cols,
    "assignments"
  )

  for (column in required_cols) {
    validate_integerish_structure_column(
      assignments,
      column,
      "assignments"
    )
    assignments[[column]] <- as.integer(assignments[[column]])
  }

  if (any(assignments$glycan_id > input_size)) {
    cli::cli_abort(c(
      "{.arg assignments} contains out-of-range glycan identifiers.",
      "x" = "{.field glycan_id} must be between 1 and {input_size}."
    ))
  }

  assignment_keys <- paste(
    assignments$glycan_id,
    assignments$part_id,
    sep = "\r"
  )
  duplicated_rows <- duplicated(assignment_keys) |
    duplicated(assignment_keys, fromLast = TRUE)
  if (any(duplicated_rows)) {
    duplicated_pairs <- unique(assignments[
      duplicated_rows,
      c("glycan_id", "part_id"),
      drop = FALSE
    ])
    pair_labels <- paste0(
      "(",
      duplicated_pairs$glycan_id,
      ", ",
      duplicated_pairs$part_id,
      ")"
    )
    cli::cli_abort(c(
      "Each floating part can be assigned at most once.",
      "x" = "Duplicated ({.field glycan_id}, {.field part_id}) pair{?s}: {.val {pair_labels}}."
    ))
  }

  assignments
}


#' Localize selected parts in one glycan graph
#'
#' @param graph A valid glycan graph.
#' @param assignments Assignment rows for this graph.
#' @param glycan_id Input glycan identifier used in errors.
#' @returns A valid, canonicalized glycan graph.
#' @noRd
localize_floating_graph <- function(graph, assignments, glycan_id = 1L) {
  parts <- normalize_floating_parts(graph)
  info <- floating_graph_info(graph, parts)

  nonexistent <- assignments$part_id > length(parts)
  if (any(nonexistent)) {
    part_ids <- unique(assignments$part_id[nonexistent])
    cli::cli_abort(c(
      "Assignment refers to a nonexistent floating part.",
      "x" = "Glycan {.val {glycan_id}} has no floating part with {.field part_id} {.val {part_ids}}."
    ))
  }

  for (row_id in seq_len(nrow(assignments))) {
    part_id <- assignments$part_id[[row_id]]
    parent_node <- assignments$parent_node[[row_id]]
    part <- parts[[part_id]]
    candidate_parents <- if (length(part$parents) == 0) {
      info$main_vertices
    } else {
      part$parents
    }

    if (!parent_node %in% candidate_parents) {
      cli::cli_abort(c(
        "Assignment selects a parent outside the candidate domain.",
        "x" = "Glycan {.val {glycan_id}} floating part {.val {part_id}} cannot attach to parent node {.val {parent_node}}.",
        "i" = "Candidate parent domain: {.val {candidate_parents}}."
      ))
    }
  }

  selected_parts <- parts
  for (row_id in seq_len(nrow(assignments))) {
    part_id <- assignments$part_id[[row_id]]
    selected_parts[[part_id]]$parents <- assignments$parent_node[[row_id]]
  }
  validation_graph <- set_floating_parts_attr(graph, selected_parts)
  validate_glycan_graph(validation_graph)

  for (row_id in seq_len(nrow(assignments))) {
    part_id <- assignments$part_id[[row_id]]
    part <- parts[[part_id]]
    graph <- igraph::add_edges(
      graph,
      c(assignments$parent_node[[row_id]], part$root),
      linkage = part$linkage
    )
  }

  remaining <- parts[-assignments$part_id]
  remaining <- purrr::map(remaining, function(part) {
    if (length(part$parents) == 0) {
      part$parents <- as.integer(info$main_vertices)
    }
    part
  })
  graph <- set_floating_parts_attr(graph, remaining)
  graph <- validate_glycan_graph(graph)
  canonicalize_glycan_graph(graph)
}
