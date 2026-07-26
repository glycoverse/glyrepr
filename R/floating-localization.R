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
#'   parent_node = 2L
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


#' Enumerate Floating-Part Localizations
#'
#' @description
#' `enumerate_floating_localizations()` generates every conflict-free,
#' fully localized structure permitted by the candidate-parent domains in `x`.
#' Each variant records the complete assignment that produced it.
#'
#' Candidate combinations are validated simultaneously, including linkages
#' with multiple possible acceptor positions such as `"a2-3/6"`. Variants are
#' canonicalized and then deduplicated by structure; when multiple assignments
#' produce the same canonical structure, the first assignment in deterministic
#' candidate order is retained.
#'
#' `max_variants` is a conservative per-input safeguard. It limits the raw
#' Cartesian product before conflict filtering or canonical deduplication, so
#' the function may ask for a higher bound even when fewer variants would
#' ultimately remain.
#'
#' Missing inputs and structures without floating parts each produce one row
#' with `variant_id = 1L`, the original structure, and an empty assignment
#' table. Empty structure vectors produce a zero-row result.
#'
#' @param x A glycan structure vector.
#' @param max_variants A positive integer giving the maximum raw candidate
#'   combinations allowed for each input structure.
#'
#' @returns A tibble with columns:
#'
#' - `input_id`: the integer position in `x`.
#' - `variant_id`: the sequential identifier after canonical deduplication.
#' - `structure`: a `glyrepr_structure` vector column containing fully
#'   localized variants.
#' - `assignments`: a list-column of tibbles with `glycan_id`, `part_id`, and
#'   `parent_node`. Here `glycan_id` equals `input_id`.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
#' )
#' enumerate_floating_localizations(glycan)
#'
#' @export
enumerate_floating_localizations <- function(
  x,
  max_variants = 256
) {
  checkmate::assert_class(x, "glyrepr_structure")
  checkmate::assert_count(max_variants, positive = TRUE)
  if (max_variants > .Machine$integer.max) {
    cli::cli_abort(
      "{.arg max_variants} must not exceed {.val {(.Machine$integer.max)}}."
    )
  }
  max_variants <- as.integer(max_variants)

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
        max_variants
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
    parent_node = integer()
  )
}


#' Enumerate localizations for one structure
#'
#' @param x A length-one glycan structure vector.
#' @param input_id Original input position.
#' @param max_variants Maximum raw candidate combinations.
#' @returns A localization result tibble for one input.
#' @noRd
enumerate_floating_localizations_one <- function(
  x,
  input_id,
  max_variants
) {
  if (is.na(x)) {
    return(single_identity_localization(x, input_id))
  }

  graph <- as.list(x)[[1]]
  parts <- normalize_floating_parts(graph)
  if (length(parts) == 0) {
    return(single_identity_localization(x, input_id))
  }

  main_vertices <- floating_main_vertices(graph, parts)
  candidate_domains <- purrr::map(parts, function(part) {
    if (length(part$parents) == 0) {
      main_vertices
    } else {
      as.integer(part$parents)
    }
  })
  combination_count <- prod(as.double(lengths(candidate_domains)))
  if (!is.finite(combination_count) || combination_count > max_variants) {
    combination_label <- format(combination_count, scientific = FALSE)
    cli::cli_abort(c(
      "Floating localization count exceeds {.arg max_variants}.",
      "x" = "Input {.val {input_id}} has {combination_label} raw candidate combinations; the limit is {.val {max_variants}}.",
      "i" = "Increase {.arg max_variants} to enumerate this input."
    ))
  }

  names(candidate_domains) <- paste0(
    "part_",
    seq_along(candidate_domains)
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
    assignments <- tibble::tibble(
      glycan_id = rep(as.integer(input_id), length(parts)),
      part_id = seq_along(parts),
      parent_node = as.integer(
        unlist(combinations[combination_id, , drop = FALSE])
      )
    )

    is_valid <- tryCatch(
      {
        validate_floating_assignment_compatibility(
          graph,
          parts,
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
      list(localize_floating_graph(graph, assignments, input_id))
    )
    assignment_tables <- append(assignment_tables, list(assignments))
  }

  if (length(localized_graphs) == 0) {
    cli::cli_abort(
      "Input {.val {input_id}} has no conflict-free floating localization."
    )
  }

  iupacs <- purrr::map_chr(localized_graphs, graph_to_iupac)
  unique_variants <- !duplicated(iupacs)
  iupacs <- iupacs[unique_variants]
  localized_graphs <- localized_graphs[unique_variants]
  assignment_tables <- assignment_tables[unique_variants]
  names(localized_graphs) <- iupacs
  structures <- new_glycan_structure(iupacs, localized_graphs)
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
#' @returns A valid, canonicalized glycan graph.
#' @noRd
localize_floating_graph <- function(graph, assignments, glycan_id = 1L) {
  parts <- normalize_floating_parts(graph)
  main_vertices <- floating_main_vertices(graph, parts)

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
      main_vertices
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

  validate_floating_assignment_compatibility(
    graph,
    parts,
    assignments
  )

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
      part$parents <- main_vertices
    }
    part
  })
  graph <- set_floating_parts_attr(graph, remaining)
  graph <- validate_glycan_graph(graph)
  canonicalize_glycan_graph(graph)
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
