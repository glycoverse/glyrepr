#' Localize Floating Glycan Parts
#'
#' @description
#' `localize_floating_parts()` attaches selected floating parts to
#' caller-supplied parent nodes. The assignments are interpreted against the
#' floating-part rows and canonical node identifiers returned by
#' [structure_floating_candidates()].
#'
#' Each selected parent must belong to the floating part's declared candidate
#' domain. For an unrestricted `{<floating>}` part, that domain contains every
#' feasible node outside its own component, including nodes in other floating
#' components. Assignments must also be simultaneously compatible with occupied
#' and potential acceptor linkage positions, acyclic, and ultimately connected
#' to the main tree.
#'
#' Selected virtual attachments become ordinary graph edges. Unassigned
#' floating parts remain floating. The resulting structures are canonicalized,
#' and candidate-parent indices for remaining parts are remapped to the new
#' canonical complete-sequence order. Attaching one floating component to
#' another merges their component metadata and can iteratively resolve newly
#' singleton domains.
#'
#' Missing values, vector positions, and names in structure-vector input are
#' preserved. For graph input, `glycan_id` must be `1L`, and selected edges are
#' appended without canonicalizing or renumbering vertices.
#'
#' @param x A glycan structure vector or a glycan `igraph`.
#' @param assignments A data frame with integer columns `glycan_id`, `part_id`,
#'   and `parent_node`. Each `(glycan_id, part_id)` pair may occur at most once.
#'
#' @returns An object of the same representation as `x`. Structure-vector
#'   output has the same length and names as `x`; graph output retains its
#'   vertex IDs and order.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' assignments <- tibble::tibble(
#'   glycan_id = 1L,
#'   part_id = 1L,
#'   parent_node = 2L
#' )
#' localize_floating_parts(glycan, assignments)
#'
#' @export
localize_floating_parts <- function(x, assignments) {
  is_graph <- inherits(x, "igraph")
  if (!is_graph) {
    checkmate::assert_class(x, "glyrepr_structure")
  }
  input_size <- if (is_graph) 1L else length(x)
  assignments <- validate_floating_assignments(assignments, input_size)

  if (nrow(assignments) == 0) {
    return(x)
  }

  if (is_graph) {
    return(localize_floating_graph(
      x,
      assignments,
      glycan_id = 1L,
      canonicalize = FALSE
    ))
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


#' Enumerate Floating Graph Localizations
#'
#' @description
#' `enumerate_floating_graph_localizations()` generates every conflict-free,
#' fully localized graph permitted by the candidate-parent domains in `graph`.
#' It is the graph-level counterpart to [enumerate_floating_localizations()].
#'
#' The returned graphs are validated but are not canonicalized. Adding the
#' selected attachment edges does not reorder vertices, so every vertex keeps
#' the same integer ID, name, and attributes as in `graph`. Assignment
#' `parent_node` values therefore refer directly to vertex IDs in both the
#' input and localized graphs.
#'
#' Floating parts are localized as ordinary edges, while floating substituents
#' are localized into the selected parent vertex's `sub` attribute. Every
#' conflict-free acyclic assignment that connects all components to the main
#' tree is retained, even when multiple assignments would produce the same
#' canonical IUPAC-condensed structure. A graph without floating metadata
#' produces one identity row with an empty assignment table.
#'
#' `max_variants` is a conservative safeguard. It limits the raw Cartesian
#' product before conflict filtering, so the function may ask for a higher
#' bound even when fewer variants would ultimately remain.
#'
#' @param graph A valid glycan `igraph`, optionally containing unresolved
#'   floating parts or substituents.
#' @param max_variants A positive integer giving the maximum raw candidate
#'   combinations allowed.
#'
#' @returns A tibble with columns:
#'
#' - `variant_id`: the sequential localization identifier.
#' - `graph`: a list-column of fully localized `igraph` objects whose vertex
#'   IDs are identical to those in `graph`.
#' - `assignments`: a list-column of tibbles with `glycan_id`, `part_id`,
#'   `parent_node`, and `substituent_id`. Exactly one of `part_id` and
#'   `substituent_id` is non-missing in each row. `glycan_id` is always `1L`.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' graph <- get_structure_graphs(glycan, return_list = FALSE)
#' localizations <- enumerate_floating_graph_localizations(graph)
#' localizations$graph
#'
#' @template low-level-structure-pipeline
#'
#' @seealso [enumerate_floating_localizations()]
#' @export
enumerate_floating_graph_localizations <- function(
  graph,
  max_variants = 256
) {
  error_call <- rlang::current_call()
  graph <- validate_glycan_graph(graph)
  max_variants <- validate_floating_max_variants(max_variants)

  variants <- enumerate_floating_graph_localizations_one(
    graph,
    max_variants,
    error_call = error_call
  )
  tibble::tibble(
    variant_id = seq_along(variants$graphs),
    graph = variants$graphs,
    assignments = variants$assignments
  )
}


#' Enumerate Floating Localizations
#'
#' @description
#' `enumerate_floating_localizations()` generates every conflict-free,
#' fully localized structure permitted by the candidate-parent domains in `x`.
#' Each variant records the complete assignment that produced it.
#'
#' Candidate combinations for floating parts and substituents are validated
#' simultaneously, including linkages or substituents with multiple possible
#' carbon positions. Floating-component dependencies must be acyclic and
#' ultimately connect to the main tree. Variants are canonicalized and, by
#' default, deduplicated by structure. When multiple assignments produce the same canonical
#' structure, the first assignment in deterministic candidate order is
#' retained. Set `deduplicate = FALSE` to retain every valid assignment and its
#' original-node provenance, including assignments that produce identical
#' canonical structures.
#'
#' `max_variants` is a conservative per-input safeguard. It limits the raw
#' Cartesian product before conflict filtering or canonical deduplication, so
#' the function may ask for a higher bound even when fewer variants would
#' ultimately remain.
#'
#' Missing inputs and structures without floating metadata each produce one row
#' with `variant_id = 1L`, the original structure, and an empty assignment
#' table. Empty structure vectors produce a zero-row result.
#'
#' @param x A glycan structure vector.
#' @param max_variants A positive integer giving the maximum raw candidate
#'   combinations allowed for each input structure.
#' @param deduplicate A logical value. If `TRUE` (default), retain only the
#'   first assignment for each canonical structure. If `FALSE`, retain every
#'   conflict-free assignment.
#'
#' @returns A tibble with columns:
#'
#' - `input_id`: the integer position in `x`.
#' - `variant_id`: the sequential identifier after optional canonical
#'   deduplication.
#' - `structure`: a `glyrepr_structure` vector column containing fully
#'   localized variants.
#' - `assignments`: a list-column of tibbles with `glycan_id`, `part_id`,
#'   `parent_node`, and `substituent_id`. Exactly one of `part_id` and
#'   `substituent_id` is non-missing in each row. Here `glycan_id` equals
#'   `input_id`.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' enumerate_floating_localizations(glycan)
#'
#' @export
enumerate_floating_localizations <- function(
  x,
  max_variants = 256,
  deduplicate = TRUE
) {
  checkmate::assert_class(x, "glyrepr_structure")
  max_variants <- validate_floating_max_variants(max_variants)
  checkmate::assert_flag(deduplicate)

  if (length(x) == 0) {
    return(tibble::tibble(
      input_id = integer(),
      variant_id = integer(),
      structure = x,
      assignments = list()
    ))
  }

  result_tables <- lapply(
    seq_along(x),
    function(input_id) {
      enumerate_floating_localizations_one(
        x[input_id],
        input_id,
        max_variants,
        deduplicate
      )
    }
  )
  dplyr::bind_rows(result_tables)
}


#' Create an empty floating assignment tibble
#'
#' @returns A zero-row tibble with localization assignment columns.
#' @noRd
empty_floating_assignments <- function() {
  tibble::tibble(
    glycan_id = integer(),
    part_id = integer(),
    parent_node = integer(),
    substituent_id = integer()
  )
}


#' Enumerate localizations for one structure
#'
#' @param x A length-one glycan structure vector.
#' @param input_id Original input position.
#' @param max_variants Maximum raw candidate combinations.
#' @param deduplicate Whether to deduplicate canonical structures.
#' @returns A localization result tibble for one input.
#' @noRd
enumerate_floating_localizations_one <- function(
  x,
  input_id,
  max_variants,
  deduplicate
) {
  if (is.na(x)) {
    return(single_identity_localization(x, input_id))
  }

  graph <- as.list(x)[[1]]
  if (!has_floating_metadata(graph)) {
    return(single_identity_localization(x, input_id))
  }

  error_call <- rlang::current_call()
  variants <- enumerate_floating_graph_localizations_one(
    graph,
    max_variants,
    input_id,
    error_call
  )
  localized_graphs <- purrr::map(
    variants$graphs,
    canonicalize_glycan_graph
  )
  assignment_tables <- variants$assignments

  iupacs <- purrr::map_chr(localized_graphs, graph_to_iupac)
  if (deduplicate) {
    unique_variants <- !duplicated(iupacs)
    iupacs <- iupacs[unique_variants]
    localized_graphs <- localized_graphs[unique_variants]
    assignment_tables <- assignment_tables[unique_variants]
  }
  unique_graphs <- !duplicated(iupacs)
  structure_graphs <- localized_graphs[unique_graphs]
  names(structure_graphs) <- iupacs[unique_graphs]
  structures <- new_glycan_structure(iupacs, structure_graphs)
  if (!is.null(names(x))) {
    names(structures) <- rep(names(x), length(structures))
  }

  tibble::tibble(
    input_id = rep(as.integer(input_id), length(structures)),
    variant_id = seq_along(structures),
    structure = structures,
    assignments = assignment_tables
  )
}


#' Validate the floating localization variant bound
#'
#' @param max_variants Maximum raw candidate combinations.
#' @returns `max_variants` as an integer.
#' @noRd
validate_floating_max_variants <- function(max_variants) {
  checkmate::assert_count(max_variants, positive = TRUE)
  if (max_variants > .Machine$integer.max) {
    cli::cli_abort(
      "{.arg max_variants} must not exceed {.val {(.Machine$integer.max)}}."
    )
  }
  as.integer(max_variants)
}


#' Enumerate graph localizations while preserving vertex IDs
#'
#' @param graph A validated glycan graph.
#' @param max_variants Maximum raw candidate combinations.
#' @param input_id Optional structure input position used in assignments and
#'   diagnostics.
#' @param error_call The call to report in localization errors.
#' @returns A list with `graphs` and `assignments` elements.
#' @noRd
enumerate_floating_graph_localizations_one <- function(
  graph,
  max_variants,
  input_id = 1L,
  error_call = rlang::caller_call()
) {
  parts <- normalize_floating_parts(graph)
  substituents <- normalize_floating_substituents(graph)
  if (length(parts) == 0 && length(substituents) == 0) {
    return(list(
      graphs = list(graph),
      assignments = list(empty_floating_assignments())
    ))
  }

  part_domains <- purrr::map(
    parts,
    ~ floating_part_candidate_parents(graph, .x)
  )
  substituent_domains <- purrr::map(
    substituents,
    ~ floating_substituent_candidate_parents(graph, .x)
  )
  candidate_domains <- c(part_domains, substituent_domains)
  combination_count <- prod(as.double(lengths(candidate_domains)))
  if (!is.finite(combination_count) || combination_count > max_variants) {
    combination_label <- format(combination_count, scientific = FALSE)
    cli::cli_abort(
      c(
        "Floating localization count exceeds {.arg max_variants}.",
        "x" = "Input {.val {input_id}} has {combination_label} raw candidate combinations; the limit is {.val {max_variants}}.",
        "i" = "Increase {.arg max_variants} to enumerate this input."
      ),
      call = error_call
    )
  }

  names(candidate_domains) <- c(
    if (length(part_domains) > 0) {
      paste0("part_", seq_along(part_domains))
    } else {
      character()
    },
    if (length(substituent_domains) > 0) {
      paste0("substituent_", seq_along(substituent_domains))
    } else {
      character()
    }
  )
  combinations <- do.call(
    expand.grid,
    c(
      candidate_domains,
      list(
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    )
  )

  localized_graphs <- list()
  assignment_tables <- list()
  for (combination_id in seq_len(nrow(combinations))) {
    parent_nodes <- as.integer(
      unlist(combinations[combination_id, , drop = FALSE])
    )
    part_assignments <- tibble::tibble(
      glycan_id = rep(as.integer(input_id), length(parts)),
      part_id = as.integer(seq_along(parts)),
      parent_node = parent_nodes[seq_along(parts)],
      substituent_id = rep(NA_integer_, length(parts))
    )
    substituent_positions <- length(parts) + seq_along(substituents)
    substituent_assignments <- tibble::tibble(
      glycan_id = rep(as.integer(input_id), length(substituents)),
      part_id = rep(NA_integer_, length(substituents)),
      parent_node = parent_nodes[substituent_positions],
      substituent_id = as.integer(seq_along(substituents))
    )
    assignments <- dplyr::bind_rows(
      part_assignments,
      substituent_assignments
    )

    is_valid <- tryCatch(
      {
        validate_complete_floating_assignment_compatibility(
          graph,
          parts,
          substituents,
          assignments
        )
        TRUE
      },
      error = function(cnd) FALSE
    )
    if (!is_valid) {
      next
    }

    localized_graphs <- append(
      localized_graphs,
      list(localize_floating_graph_preserve_ids(
        graph,
        parts,
        substituents,
        assignments
      ))
    )
    assignment_tables <- append(assignment_tables, list(assignments))
  }

  if (length(localized_graphs) == 0) {
    cli::cli_abort(
      "Input {.val {input_id}} has no conflict-free floating localization.",
      call = error_call
    )
  }

  list(
    graphs = localized_graphs,
    assignments = assignment_tables
  )
}


#' Apply complete floating assignments without canonicalizing
#'
#' @param graph A validated floating glycan graph.
#' @param parts Normalized floating-part metadata.
#' @param substituents Normalized floating-substituent metadata.
#' @param assignments A complete compatible assignment table.
#' @returns A validated graph with its original vertex order.
#' @noRd
localize_floating_graph_preserve_ids <- function(
  graph,
  parts,
  substituents,
  assignments
) {
  part_assignments <- assignments[!is.na(assignments$part_id), , drop = FALSE]
  for (row_id in seq_len(nrow(part_assignments))) {
    part_id <- part_assignments$part_id[[row_id]]
    part <- parts[[part_id]]
    graph <- igraph::add_edges(
      graph,
      c(part_assignments$parent_node[[row_id]], part$root),
      linkage = part$linkage
    )
  }

  substituent_assignments <- assignments[
    !is.na(assignments$substituent_id),
    ,
    drop = FALSE
  ]
  for (row_id in seq_len(nrow(substituent_assignments))) {
    substituent_id <- substituent_assignments$substituent_id[[row_id]]
    graph <- append_vertex_substituent(
      graph,
      substituent_assignments$parent_node[[row_id]],
      substituents[[substituent_id]]$substituent
    )
  }

  graph <- delete_floating_parts_attr(graph)
  graph <- delete_floating_substituents_attr(graph)
  validate_glycan_graph(graph)
}


#' Validate a complete mixed floating assignment
#'
#' @param graph A valid floating glycan graph.
#' @param parts Normalized floating-part metadata.
#' @param substituents Normalized floating-substituent metadata.
#' @param assignments A complete assignment table.
#' @returns `NULL`, invisibly.
#' @noRd
validate_complete_floating_assignment_compatibility <- function(
  graph,
  parts,
  substituents,
  assignments
) {
  selected_parts <- parts
  part_assignments <- assignments[!is.na(assignments$part_id), , drop = FALSE]
  for (row_id in seq_len(nrow(part_assignments))) {
    part_id <- part_assignments$part_id[[row_id]]
    selected_parts[[part_id]]$parents <- part_assignments$parent_node[[row_id]]
  }

  selected_substituents <- substituents
  substituent_assignments <- assignments[
    !is.na(assignments$substituent_id),
    ,
    drop = FALSE
  ]
  for (row_id in seq_len(nrow(substituent_assignments))) {
    substituent_id <- substituent_assignments$substituent_id[[row_id]]
    selected_substituents[[substituent_id]]$parents <-
      substituent_assignments$parent_node[[row_id]]
  }

  validation_graph <- set_floating_parts_attr(graph, selected_parts)
  validation_graph <- set_floating_substituents_attr(
    validation_graph,
    selected_substituents
  )
  validate_glycan_graph(validation_graph)
  invisible(NULL)
}


#' Return one identity localization row
#'
#' @param x A length-one glycan structure vector.
#' @param input_id Original input position.
#' @returns A one-row localization result tibble.
#' @noRd
single_identity_localization <- function(x, input_id) {
  tibble::tibble(
    input_id = as.integer(input_id),
    variant_id = 1L,
    structure = x,
    assignments = list(empty_floating_assignments())
  )
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
#' @param canonicalize Whether to canonicalize the localized graph.
#' @returns A valid glycan graph, canonicalized when `canonicalize` is `TRUE`.
#' @noRd
localize_floating_graph <- function(
  graph,
  assignments,
  glycan_id = 1L,
  canonicalize = TRUE
) {
  parts <- normalize_floating_parts(graph)

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
    candidate_parents <- floating_part_candidate_parents(graph, part)

    if (!parent_node %in% candidate_parents) {
      cli::cli_abort(c(
        "Assignment selects a parent outside the candidate domain.",
        "x" = "Glycan {.val {glycan_id}} floating part {.val {part_id}} cannot attach to parent node {.val {parent_node}}.",
        "i" = "Candidate parent domain: {.val {candidate_parents}}."
      ))
    }
  }

  validate_floating_assignment_compatibility(
    graph,
    parts,
    assignments
  )

  graph <- attach_floating_parts(
    graph,
    parts,
    assignments$part_id,
    assignments$parent_node
  )
  graph <- validate_glycan_graph(graph)
  if (canonicalize) {
    graph <- canonicalize_glycan_graph(graph)
  }
  graph
}


#' Validate simultaneous floating-part assignments
#'
#' @param graph A valid floating glycan graph.
#' @param parts Normalized floating-part metadata.
#' @param assignments Assignment rows for this graph.
#' @returns `NULL`, invisibly.
#' @noRd
validate_floating_assignment_compatibility <- function(
  graph,
  parts,
  assignments
) {
  selected_parts <- parts
  for (row_id in seq_len(nrow(assignments))) {
    part_id <- assignments$part_id[[row_id]]
    selected_parts[[part_id]]$parents <- assignments$parent_node[[row_id]]
  }
  validation_graph <- set_floating_parts_attr(graph, selected_parts)
  validate_glycan_graph(validation_graph)
  invisible(NULL)
}
