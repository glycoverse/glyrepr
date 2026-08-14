# Internal helpers for substituents with unresolved parent residues.

#' Detect Floating Substituents
#'
#' Test whether each glycan structure contains one or more substituents whose
#' parent residue is unresolved.
#'
#' In the `glyrepr` IUPAC extension, floating substituents appear in braces
#' before the main glycan. `{6S}<main>` allows every feasible residue in the
#' complete structure to carry the substituent, `{6S|1,2}<main>` restricts it
#' to complete-sequence nodes 1 and 2, and `{?S}<main>` also leaves the carbon
#' position unknown. Residue nodes in floating blocks are counted before the
#' main glycan; substituent blocks contribute no node indices.
#'
#' Internally, unresolved substituents are stored in the
#' `floating_substituents` graph attribute. It is a list with one entry per
#' substituent. Each entry contains a canonical `substituent` token and an
#' integer `parents` vector of candidate graph vertex indices. An empty vector
#' means every feasible residue vertex is a candidate.
#'
#' A singleton candidate set is normalized into that vertex's ordinary `sub`
#' attribute. This also happens whenever chemistry leaves only one feasible
#' residue in the complete structure.
#' [structure_floating_substituents()] exposes the metadata as a normalized
#' table.
#'
#' @param x A [glycan_structure()] vector or a glycan `igraph`.
#'
#' @returns A logical vector with the same length and names as vector input, or
#'   a logical scalar for graph input. Missing structures produce `NA`.
#'
#' @examples
#' glycans <- as_glycan_structure(c(
#'   floating = "{6S}Gal(a1-3)Gal(a1-",
#'   restricted = "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-",
#'   ordinary = "Gal6S(a1-"
#' ))
#' has_floating_substituents(glycans)
#'
#' @seealso [structure_floating_substituents()], [has_floating_parts()],
#'   [as_glycan_structure()]
#' @export
has_floating_substituents <- function(x) {
  if (inherits(x, "igraph")) {
    return(has_floating_substituents_graph(x))
  }
  checkmate::assert_class(x, "glyrepr_structure")
  smap_lgl(x, has_floating_substituents_graph)
}

floating_substituents_attr <- function(graph) {
  substituents <- igraph::graph_attr(graph, "floating_substituents")
  if (is.null(substituents)) {
    return(list())
  }
  substituents
}

has_floating_substituents_graph <- function(graph) {
  length(floating_substituents_attr(graph)) > 0
}

has_floating_metadata <- function(graph) {
  has_floating_parts_graph(graph) || has_floating_substituents_graph(graph)
}

delete_floating_substituents_attr <- function(graph) {
  if ("floating_substituents" %in% igraph::graph_attr_names(graph)) {
    graph <- igraph::delete_graph_attr(graph, "floating_substituents")
  }
  graph
}

set_floating_substituents_attr <- function(graph, substituents) {
  graph <- delete_floating_substituents_attr(graph)
  if (length(substituents) > 0) {
    graph <- igraph::set_graph_attr(
      graph,
      "floating_substituents",
      value = substituents
    )
  }
  graph
}

normalize_floating_substituent <- function(substituent, graph_size) {
  if (!is.list(substituent)) {
    cli::cli_abort("Each floating substituent must be a list.")
  }

  allowed <- c("substituent", "parents")
  fields <- names(substituent)
  if (is.null(fields)) {
    fields <- rep("", length(substituent))
  }
  unsupported <- unique(fields[
    is.na(fields) | fields == "" | !fields %in% allowed
  ])
  unsupported[is.na(unsupported) | unsupported == ""] <- "<unnamed>"
  if (length(unsupported) > 0) {
    cli::cli_abort(c(
      "Floating substituent metadata contains unsupported fields.",
      "x" = "Unsupported field{?s}: {.field {unsupported}}."
    ))
  }

  duplicated_fields <- unique(fields[duplicated(fields)])
  if (length(duplicated_fields) > 0) {
    cli::cli_abort(c(
      "Floating substituent metadata contains duplicated fields.",
      "x" = "Duplicated field{?s}: {.field {duplicated_fields}}."
    ))
  }

  missing <- setdiff(allowed, names(substituent))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Floating substituent metadata is incomplete.",
      "x" = "Missing field{?s}: {.field {missing}}."
    ))
  }

  token <- substituent$substituent
  pattern <- substituent_token_pattern(anchored = TRUE)
  if (
    !is.character(token) ||
      length(token) != 1 ||
      is.na(token) ||
      !stringr::str_detect(token, pattern)
  ) {
    cli::cli_abort(
      "Each floating substituent must contain one valid token in {.field substituent}."
    )
  }
  if (!identical(token, normalize_substituent_token(token))) {
    cli::cli_abort(
      "Floating substituent field {.field substituent} must use canonical position order."
    )
  }

  parents <- substituent$parents
  if (
    !is.numeric(parents) ||
      anyNA(parents) ||
      any(!is.finite(parents)) ||
      any(parents > .Machine$integer.max) ||
      any(parents != as.integer(parents)) ||
      any(parents < 1) ||
      any(parents > graph_size)
  ) {
    cli::cli_abort(
      "Floating substituent parent indices must refer to valid graph vertices."
    )
  }
  parents <- as.integer(parents)
  if (anyDuplicated(parents) > 0) {
    cli::cli_abort("Floating substituent parent indices must be unique.")
  }

  list(substituent = token, parents = parents)
}

normalize_floating_substituents <- function(graph) {
  substituents <- floating_substituents_attr(graph)
  if (!is.list(substituents)) {
    cli::cli_abort(
      "Graph attribute {.field floating_substituents} must be a list."
    )
  }

  purrr::map(
    substituents,
    normalize_floating_substituent,
    graph_size = igraph::vcount(graph)
  )
}

floating_metadata_main_vertices <- function(
  graph,
  parts = normalize_floating_parts(graph)
) {
  if (length(parts) == 0) {
    return(as.integer(seq_len(igraph::vcount(graph))))
  }
  floating_main_vertices(graph, parts)
}

validate_floating_substituent_parents <- function(graph) {
  raw_substituents <- igraph::graph_attr(graph, "floating_substituents")
  if (is.null(raw_substituents)) {
    return(invisible(NULL))
  }

  substituents <- normalize_floating_substituents(graph)
  if (length(substituents) == 0) {
    return(invisible(NULL))
  }

  invisible(NULL)
}

floating_substituent_candidate_parents <- function(graph, substituent) {
  if (length(substituent$parents) > 0) {
    return(as.integer(substituent$parents))
  }

  as.integer(seq_len(igraph::vcount(graph)))
}

main_substituent_domains <- function(
  graph,
  vertices = seq_len(igraph::vcount(graph))
) {
  domains <- list()
  for (vertex in vertices) {
    sub <- igraph::vertex_attr(graph, "sub", index = vertex)
    if (identical(sub, "")) {
      next
    }
    tokens <- stringr::str_split(sub, stringr::fixed(","))[[1]]
    for (token in tokens) {
      positions <- substituent_position_tokens(token)
      if (identical(positions, "?")) {
        next
      }
      values <- stringr::str_split(positions, stringr::fixed("/"))[[1]]
      domains[[length(domains) + 1L]] <- list(
        known = TRUE,
        slots = paste(vertex, values, sep = "\r")
      )
    }
  }
  domains
}

floating_substituent_domain <- function(
  substituent,
  substituent_id,
  candidate_parents,
  occupied_slots
) {
  position <- substituent_position_tokens(substituent$substituent)
  if (identical(position, "?")) {
    return(list(known = FALSE, slots = character()))
  }

  positions <- stringr::str_split(position, stringr::fixed("/"))[[1]]
  explicit <- length(substituent$parents) > 0
  parents <- candidate_parents
  slots_by_parent <- lapply(parents, function(parent) {
    candidate_slots <- paste(parent, positions, sep = "\r")
    setdiff(candidate_slots, occupied_slots)
  })

  if (explicit) {
    impossible <- parents[lengths(slots_by_parent) == 0]
    if (length(impossible) > 0) {
      cli::cli_abort(c(
        "Floating substituent has impossible explicit parent metadata.",
        "x" = "Floating substituent {substituent_id} cannot use explicit parent node{?s} {.val {impossible}} because every declared carbon position is already occupied."
      ))
    }
  }

  list(
    known = TRUE,
    slots = unique(unlist(slots_by_parent, use.names = FALSE))
  )
}

validate_floating_substituent_slots <- function(graph) {
  substituents <- normalize_floating_substituents(graph)
  if (length(substituents) == 0) {
    return(invisible(NULL))
  }

  parts <- normalize_floating_parts(graph)
  main_edge_domains <- main_attachment_domains(graph)
  main_sub_domains <- main_substituent_domains(graph)
  occupied_slots <- definitely_occupied_main_slots(c(
    main_edge_domains,
    main_sub_domains
  ))
  floating_part_domains <- lapply(seq_along(parts), function(part_id) {
    part <- parts[[part_id]]
    floating_attachment_domain(
      part,
      part_id,
      floating_part_candidate_parents(graph, part),
      occupied_slots
    )
  })
  floating_sub_domains <- lapply(
    seq_along(substituents),
    function(substituent_id) {
      substituent <- substituents[[substituent_id]]
      floating_substituent_domain(
        substituent,
        substituent_id,
        floating_substituent_candidate_parents(graph, substituent),
        occupied_slots
      )
    }
  )
  domains <- c(
    main_edge_domains,
    main_sub_domains,
    floating_part_domains,
    floating_sub_domains
  )

  if (!has_conflict_free_attachment_assignment(domains)) {
    cli::cli_abort(c(
      "Floating substituents cannot be localized simultaneously.",
      "x" = "No conflict-free assignment exists for the declared parent residues and carbon positions."
    ))
  }

  invisible(NULL)
}

append_vertex_substituent <- function(graph, vertex, substituent) {
  current <- igraph::vertex_attr(graph, "sub", index = vertex)
  tokens <- if (identical(current, "")) {
    substituent
  } else {
    c(stringr::str_split(current, stringr::fixed(","))[[1]], substituent)
  }
  igraph::set_vertex_attr(
    graph,
    "sub",
    index = vertex,
    value = collapse_substituent_tokens(tokens)
  )
}

resolve_single_parent_floating_substituents <- function(graph) {
  substituents <- normalize_floating_substituents(graph)
  if (length(substituents) == 0) {
    return(graph)
  }

  candidate_parents <- purrr::map(
    substituents,
    ~ floating_substituent_candidate_parents(graph, .x)
  )
  resolved <- lengths(candidate_parents) == 1
  if (!any(resolved)) {
    return(graph)
  }

  for (substituent_id in which(resolved)) {
    substituent <- substituents[[substituent_id]]
    graph <- append_vertex_substituent(
      graph,
      candidate_parents[[substituent_id]][[1]],
      substituent$substituent
    )
  }

  set_floating_substituents_attr(graph, substituents[!resolved])
}

materialize_unrestricted_floating_substituents <- function(
  graph,
  vertices = seq_len(igraph::vcount(graph))
) {
  substituents <- normalize_floating_substituents(graph)
  substituents <- purrr::map(substituents, function(substituent) {
    if (length(substituent$parents) == 0) {
      substituent$parents <- as.integer(vertices)
    }
    substituent
  })
  set_floating_substituents_attr(graph, substituents)
}

floating_substituent_iupac <- function(
  substituent
) {
  parents <- if (length(substituent$parents) == 0) {
    ""
  } else {
    paste0("|", paste(substituent$parents, collapse = ","))
  }
  paste0("{", substituent$substituent, parents, "}")
}
