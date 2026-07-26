# Internal helpers for floating glycan substructures.

#' Detect Floating Glycan Parts
#'
#' @description
#' Test whether each glycan structure contains one or more unresolved floating
#' parts.
#'
#' @details
#' A floating part is a known glycan residue or substructure whose parent
#' residue on the main glycan tree is not fully localized. For example, a
#' bi-antennary N-glycan may contain one sialic acid while the available
#' evidence cannot determine which of its two terminal galactoses carries that
#' residue. The sialic acid can then be represented as a floating part with
#' both galactoses as candidate parents.
#'
#' In the `glyrepr` IUPAC extension, floating parts appear in braces before the
#' main glycan:
#'
#' - `{<floating>}<main>` means every feasible main-tree node is a candidate
#'   parent.
#' - `{<floating>|<parents>}<main>` restricts the candidates to the
#'   comma-separated main-tree node indices in `<parents>`.
#'
#' Main-tree indices follow residue order in the IUPAC-condensed sequence. A
#' floating part may contain one residue or an entire subtree, and its virtual
#' attachment linkage may be fully known, partially known, or unknown.
#'
#' Internally, an unresolved structure is an annotated forest containing one
#' main tree and one disconnected tree per floating part. The virtual linkage
#' and candidate parents are stored as graph metadata rather than as an edge;
#' [structure_floating_parts()] exposes this metadata in tabular form.
#'
#' An explicit singleton parent list fully localizes the attachment, as does an
#' omitted parent list when the main tree has only one node. Such a part is
#' normalized to an ordinary graph edge, so `has_floating_parts()` returns
#' `FALSE` for the normalized structure.
#'
#' @param x A [glycan_structure()] vector.
#'
#' @returns A logical vector with the same length and names as `x`. Missing
#'   structures produce `NA`.
#'
#' @examples
#' main <- paste0(
#'   "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
#'   "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
#'   "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' )
#' ambiguous <- as_glycan_structure(
#'   paste0("{Neu5Ac(a2-3)|1,4}", main)
#' )
#' glycans <- c(ambiguous = ambiguous, ordinary = as_glycan_structure(main))
#' has_floating_parts(glycans)
#' structure_floating_parts(ambiguous)
#'
#' localized <- as_glycan_structure(
#'   "{Neu5Ac(a2-3)|1}Gal(b1-4)GlcNAc(b1-"
#' )
#' has_floating_parts(localized)
#'
#' @seealso [structure_floating_parts()], [as_glycan_structure()]
#' @export
has_floating_parts <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")
  smap_lgl(x, has_floating_parts_graph)
}

floating_parts_attr <- function(graph) {
  parts <- igraph::graph_attr(graph, "floating_parts")
  if (is.null(parts)) {
    return(list())
  }
  parts
}

has_floating_parts_graph <- function(graph) {
  length(floating_parts_attr(graph)) > 0
}

delete_floating_parts_attr <- function(graph) {
  if ("floating_parts" %in% igraph::graph_attr_names(graph)) {
    graph <- igraph::delete_graph_attr(graph, "floating_parts")
  }
  graph
}

set_floating_parts_attr <- function(graph, parts) {
  graph <- delete_floating_parts_attr(graph)
  if (length(parts) > 0) {
    graph <- igraph::set_graph_attr(graph, "floating_parts", value = parts)
  }
  graph
}

normalize_floating_part <- function(part, graph_size) {
  if (!is.list(part)) {
    cli::cli_abort("Each floating part must be a list.")
  }

  allowed <- c("root", "linkage", "parents")
  fields <- names(part)
  if (is.null(fields)) {
    fields <- rep("", length(part))
  }
  unsupported <- unique(fields[
    is.na(fields) | fields == "" | !fields %in% allowed
  ])
  unsupported[is.na(unsupported) | unsupported == ""] <- "<unnamed>"
  if (length(unsupported) > 0) {
    cli::cli_abort(c(
      "Floating part metadata contains unsupported fields.",
      "x" = "Unsupported field{?s}: {.field {unsupported}}."
    ))
  }

  duplicated_fields <- unique(fields[duplicated(fields)])
  if (length(duplicated_fields) > 0) {
    cli::cli_abort(c(
      "Floating part metadata contains duplicated fields.",
      "x" = "Duplicated field{?s}: {.field {duplicated_fields}}."
    ))
  }

  required <- c("root", "linkage", "parents")
  missing <- setdiff(required, names(part))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Floating part metadata is incomplete.",
      "x" = "Missing field{?s}: {.field {missing}}."
    ))
  }

  root <- part$root
  if (
    length(root) != 1 ||
      is.na(root) ||
      !is.numeric(root) ||
      !is.finite(root) ||
      root > .Machine$integer.max ||
      root != as.integer(root) ||
      root < 1 ||
      root > graph_size
  ) {
    cli::cli_abort(
      "Floating part {.field root} must be one valid vertex index."
    )
  }

  linkage <- part$linkage
  if (
    !is.character(linkage) ||
      length(linkage) != 1 ||
      is.na(linkage) ||
      !valid_linkages(linkage)
  ) {
    cli::cli_abort(
      "Floating part {.field linkage} must be one valid linkage."
    )
  }

  parents <- part$parents
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
      "Floating part {.field parents} must contain valid vertex indices."
    )
  }
  parents <- as.integer(parents)
  if (anyDuplicated(parents) > 0) {
    cli::cli_abort("Floating part parent indices must be unique.")
  }

  list(
    root = as.integer(root),
    linkage = linkage,
    parents = parents
  )
}

normalize_floating_parts <- function(graph) {
  parts <- floating_parts_attr(graph)
  if (!is.list(parts)) {
    cli::cli_abort(
      "Graph attribute {.field floating_parts} must be a list."
    )
  }

  purrr::map(parts, normalize_floating_part, graph_size = igraph::vcount(graph))
}

floating_graph_info <- function(
  graph,
  parts = normalize_floating_parts(graph)
) {
  components <- igraph::components(graph, mode = "weak")
  membership <- components$membership
  floating_components <- membership[purrr::map_int(parts, "root")]

  if (anyDuplicated(floating_components) > 0) {
    cli::cli_abort(
      "Each floating part must identify a different graph component."
    )
  }

  main_components <- setdiff(unique(membership), floating_components)
  if (length(main_components) != 1) {
    cli::cli_abort(
      "A floating glycan structure must contain exactly one main tree."
    )
  }

  main_component <- main_components[[1]]
  main_vertices <- which(membership == main_component)

  list(
    membership = membership,
    main_component = main_component,
    main_vertices = main_vertices,
    floating_components = floating_components,
    parts = parts
  )
}

validate_floating_graph_shape <- function(graph) {
  parts <- normalize_floating_parts(graph)

  if (length(parts) == 0) {
    if (!is_out_tree(graph)) {
      cli::cli_abort("Glycan structure must be an out tree.")
    }
    return(invisible(NULL))
  }

  info <- floating_graph_info(graph, parts)
  component_ids <- unique(info$membership)
  if (length(component_ids) != length(parts) + 1) {
    cli::cli_abort(
      "Every non-main graph component must be a declared floating part."
    )
  }

  for (component in component_ids) {
    vertices <- which(info$membership == component)
    subgraph <- igraph::induced_subgraph(graph, vertices)
    if (!is_out_tree(subgraph)) {
      cli::cli_abort(
        "The main glycan and every floating part must be an out tree."
      )
    }
  }

  roots <- purrr::map_int(parts, "root")
  indegree <- igraph::degree(graph, mode = "in")
  if (any(indegree[roots] != 0)) {
    cli::cli_abort(
      "Each floating part root must be the root of its graph component."
    )
  }

  for (part in parts) {
    if (
      length(part$parents) > 0 &&
        any(!part$parents %in% info$main_vertices)
    ) {
      cli::cli_abort(
        "Floating part parent indices must refer to vertices in the main tree."
      )
    }
  }

  validate_floating_attachment_slots(graph, info)

  invisible(NULL)
}

floating_linkage_acceptor_positions <- function(linkage) {
  acceptor <- sub(".*-", "", linkage)
  if (identical(acceptor, "?")) {
    return(integer())
  }

  unique(as.integer(strsplit(acceptor, "/", fixed = TRUE)[[1]]))
}

definitely_occupied_main_slots <- function(graph, main_vertices) {
  if (igraph::ecount(graph) == 0) {
    return(character())
  }

  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  linkages <- igraph::edge_attr(graph, "linkage")
  is_main_edge <- endpoints[, 1] %in%
    main_vertices &
    endpoints[, 2] %in% main_vertices

  slots <- character()
  for (edge_id in which(is_main_edge)) {
    positions <- floating_linkage_acceptor_positions(linkages[[edge_id]])
    if (length(positions) == 1) {
      slots <- c(
        slots,
        paste(endpoints[edge_id, 1], positions, sep = "\r")
      )
    }
  }

  unique(slots)
}

floating_attachment_domain <- function(
  part,
  part_id,
  main_vertices,
  occupied_slots
) {
  positions <- floating_linkage_acceptor_positions(part$linkage)
  if (length(positions) == 0) {
    return(list(known = FALSE, slots = character()))
  }

  explicit <- length(part$parents) > 0
  parents <- if (explicit) part$parents else main_vertices
  slots_by_parent <- lapply(parents, function(parent) {
    candidate_slots <- paste(parent, positions, sep = "\r")
    setdiff(candidate_slots, occupied_slots)
  })

  if (explicit) {
    impossible <- parents[lengths(slots_by_parent) == 0]
    if (length(impossible) > 0) {
      cli::cli_abort(c(
        "Floating part has impossible explicit parent metadata.",
        "x" = "Floating part {part_id} cannot use explicit parent node{?s} {.val {impossible}} because every acceptor position declared by linkage {.val {part$linkage}} is already occupied."
      ))
    }
  }

  list(
    known = TRUE,
    slots = unique(unlist(slots_by_parent, use.names = FALSE))
  )
}

has_conflict_free_attachment_assignment <- function(domains) {
  domains <- domains[vapply(domains, `[[`, logical(1), "known")]
  if (length(domains) == 0) {
    return(TRUE)
  }

  slot_sets <- lapply(domains, `[[`, "slots")
  if (any(lengths(slot_sets) == 0)) {
    return(FALSE)
  }
  slot_sets <- slot_sets[order(lengths(slot_sets))]

  assign_slot <- function(index, used) {
    if (index > length(slot_sets)) {
      return(TRUE)
    }

    available <- setdiff(slot_sets[[index]], used)
    for (slot in available) {
      if (assign_slot(index + 1, c(used, slot))) {
        return(TRUE)
      }
    }
    FALSE
  }

  assign_slot(1, character())
}

validate_floating_attachment_slots <- function(graph, info) {
  if (
    !has_edge_attrs(graph, "linkage") ||
      !is.character(igraph::edge_attr(graph, "linkage"))
  ) {
    return(invisible(NULL))
  }

  linkages <- igraph::edge_attr(graph, "linkage")
  if (
    length(linkages) != igraph::ecount(graph) ||
      anyNA(linkages) ||
      !all(valid_linkages(linkages))
  ) {
    return(invisible(NULL))
  }

  occupied_slots <- definitely_occupied_main_slots(
    graph,
    info$main_vertices
  )
  domains <- lapply(seq_along(info$parts), function(part_id) {
    floating_attachment_domain(
      info$parts[[part_id]],
      part_id,
      info$main_vertices,
      occupied_slots
    )
  })

  if (!has_conflict_free_attachment_assignment(domains)) {
    cli::cli_abort(c(
      "Floating parts cannot be attached simultaneously.",
      "x" = "No conflict-free assignment exists for the declared parent and acceptor positions."
    ))
  }

  invisible(NULL)
}

component_sequence_order <- function(graph, vertices) {
  subgraph <- igraph::induced_subgraph(graph, vertices)
  subgraph <- delete_floating_parts_attr(subgraph)
  cache <- build_seq_cache(subgraph)
  order <- seq_glycan_order(cache$root, cache)

  vertex_names <- igraph::V(subgraph)$name[order$vertices]
  edge_table <- igraph::as_data_frame(subgraph, what = "edges")
  edge_keys <- paste(edge_table$from, edge_table$to, sep = "\r")
  ordered_edge_keys <- edge_keys[order$edges]

  graph_edges <- igraph::as_data_frame(graph, what = "edges")
  graph_edge_keys <- paste(graph_edges$from, graph_edges$to, sep = "\r")

  list(
    vertices = match(vertex_names, igraph::V(graph)$name),
    edges = match(ordered_edge_keys, graph_edge_keys),
    iupac = seq_glycan_iupac(cache$root, cache),
    root_name = igraph::V(subgraph)$name[cache$root]
  )
}

resolve_single_parent_floating_parts <- function(graph) {
  parts <- normalize_floating_parts(graph)
  if (length(parts) == 0) {
    return(graph)
  }

  info <- floating_graph_info(graph, parts)
  candidate_parents <- purrr::map(parts, function(part) {
    if (length(part$parents) > 0) {
      part$parents
    } else {
      info$main_vertices
    }
  })
  resolved <- lengths(candidate_parents) == 1
  if (!any(resolved)) {
    return(graph)
  }

  for (part_id in which(resolved)) {
    part <- parts[[part_id]]
    graph <- igraph::add_edges(
      graph,
      c(candidate_parents[[part_id]][[1]], part$root),
      linkage = part$linkage
    )
  }

  remaining <- parts[!resolved]
  if (length(remaining) == 0) {
    return(delete_floating_parts_attr(graph))
  }

  # An unrestricted part refers to every node in the original main tree.
  # Materialize that domain before resolved components enlarge the main tree.
  remaining <- purrr::map(remaining, function(part) {
    if (length(part$parents) == 0) {
      part$parents <- as.integer(info$main_vertices)
    }
    part
  })
  set_floating_parts_attr(graph, remaining)
}

canonicalize_floating_graph <- function(graph) {
  graph <- resolve_single_parent_floating_parts(graph)
  if (!has_floating_parts_graph(graph)) {
    seq_cache <- build_seq_cache(graph)
    order <- seq_glycan_order(seq_cache$root, seq_cache)
    return(.reorder_by_sequence_order(graph, order))
  }

  parts <- normalize_floating_parts(graph)
  info <- floating_graph_info(graph, parts)
  original_names <- igraph::V(graph)$name

  main_order <- floating_symmetry_main_order(graph, info)
  main_names <- original_names[main_order$vertices]
  main_index <- stats::setNames(seq_along(main_names), main_names)

  floating_orders <- purrr::map2(
    parts,
    info$floating_components,
    function(part, component) {
      order <- component_sequence_order(
        graph,
        which(info$membership == component)
      )
      parent_names <- original_names[part$parents]
      parents <- unname(main_index[parent_names])
      parents <- sort(as.integer(parents))
      suffix <- if (length(parents) == 0) {
        ""
      } else {
        paste0("|", paste(parents, collapse = ","))
      }

      list(
        order = order,
        linkage = part$linkage,
        parents = parents,
        key = paste0(
          "{",
          order$iupac,
          "(",
          part$linkage,
          ")",
          suffix,
          "}"
        )
      )
    }
  )

  floating_order <- order(purrr::map_chr(floating_orders, "key"))
  floating_orders <- floating_orders[floating_order]

  vertex_order <- c(
    main_order$vertices,
    unlist(
      purrr::map(floating_orders, c("order", "vertices")),
      use.names = FALSE
    )
  )
  edge_order <- c(
    main_order$edges,
    unlist(
      purrr::map(floating_orders, c("order", "edges")),
      use.names = FALSE
    )
  )

  graph <- delete_floating_parts_attr(graph)
  graph <- .reorder_by_sequence_order(
    graph,
    list(vertices = vertex_order, edges = edge_order)
  )

  main_size <- length(main_order$vertices)
  offset <- main_size
  canonical_parts <- vector("list", length(floating_orders))
  for (i in seq_along(floating_orders)) {
    floating <- floating_orders[[i]]
    component_size <- length(floating$order$vertices)
    root_position <- match(
      floating$order$root_name,
      original_names[floating$order$vertices]
    )
    canonical_parts[[i]] <- list(
      root = as.integer(offset + root_position),
      linkage = floating$linkage,
      parents = floating$parents
    )
    offset <- offset + component_size
  }

  set_floating_parts_attr(graph, canonical_parts)
}

floating_part_iupac <- function(graph, part, membership) {
  component <- membership[[part$root]]
  vertices <- which(membership == component)
  subgraph <- igraph::induced_subgraph(graph, vertices)
  subgraph <- delete_floating_parts_attr(subgraph)
  cache <- build_seq_cache(subgraph)
  parents <- if (length(part$parents) == 0) {
    ""
  } else {
    paste0("|", paste(part$parents, collapse = ","))
  }

  paste0(
    "{",
    seq_glycan_iupac(cache$root, cache),
    "(",
    part$linkage,
    ")",
    parents,
    "}"
  )
}
