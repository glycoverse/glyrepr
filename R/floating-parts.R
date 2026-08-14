# Internal helpers for floating glycan substructures.

#' Detect Floating Glycan Parts
#'
#' @description
#' Test whether each glycan structure contains one or more unresolved floating
#' parts.
#'
#' @details
#' A floating part is a known glycan residue or substructure whose parent
#' residue in the complete glycan structure is not fully localized. For example, a
#' bi-antennary N-glycan may contain one sialic acid while the available
#' evidence cannot determine which of its two terminal galactoses carries that
#' residue. The sialic acid can then be represented as a floating part with
#' both galactoses as candidate parents.
#'
#' In the `glyrepr` IUPAC extension, floating parts appear in braces before the
#' main glycan:
#'
#' - `{<floating>}<main>` means every feasible node outside that floating
#'   component is a candidate parent.
#' - `{<floating>|<parents>}<main>` restricts the candidates to the
#'   comma-separated complete-sequence node indices in `<parents>`.
#'
#' Node indices follow residue order in the complete IUPAC-condensed sequence:
#' residue nodes in floating blocks are counted from left to right before the
#' main glycan, while substituent blocks contribute no node indices. A floating
#' part may target a node in another floating component or the main tree, but
#' not a node in its own component. It may contain one residue or an entire
#' subtree, and its virtual attachment linkage may be fully known, partially
#' known, or unknown.
#'
#' Internally, an unresolved structure is an annotated forest containing one
#' main tree and one disconnected tree per floating part. The virtual linkage
#' and candidate parents are stored as graph metadata rather than as an edge.
#' Each floating-part metadata entry also stores `nodes`, the complete vector of
#' vertex indices in that floating component, so downstream graph operations do
#' not need to rediscover its membership by traversal.
#' [structure_floating_parts()] exposes attachment metadata in tabular form, and
#' [structure_component_membership()] exposes component membership.
#'
#' An effective singleton parent domain fully localizes the attachment. When
#' one floating component attaches to another, their component metadata is
#' merged and newly singleton domains are resolved iteratively. Such parts are
#' normalized to ordinary graph edges, so `has_floating_parts()` returns
#' `FALSE` once every attachment is localized.
#'
#' @param x A [glycan_structure()] vector or a glycan `igraph`.
#'
#' @returns A logical vector with the same length and names as vector input, or
#'   a logical scalar for graph input. Missing structures produce `NA`.
#'
#' @examples
#' main <- paste0(
#'   "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
#'   "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
#'   "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' )
#' ambiguous <- as_glycan_structure(
#'   paste0("{Neu5Ac(a2-3)|2,5}", main)
#' )
#' glycans <- c(ambiguous = ambiguous, ordinary = as_glycan_structure(main))
#' has_floating_parts(glycans)
#' structure_floating_parts(ambiguous)
#'
#' localized <- as_glycan_structure(
#'   "{Neu5Ac(a2-3)|2}Gal(b1-4)GlcNAc(b1-"
#' )
#' has_floating_parts(localized)
#'
#' @seealso [structure_floating_parts()], [has_floating_substituents()],
#'   [as_glycan_structure()]
#' @export
has_floating_parts <- function(x) {
  if (inherits(x, "igraph")) {
    return(has_floating_parts_graph(x))
  }
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

  allowed <- c("root", "nodes", "linkage", "parents")
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
  root <- as.integer(root)

  nodes <- part$nodes
  if (!is.null(nodes)) {
    if (
      !is.numeric(nodes) ||
        length(nodes) == 0 ||
        anyNA(nodes) ||
        any(!is.finite(nodes)) ||
        any(nodes > .Machine$integer.max) ||
        any(nodes != as.integer(nodes)) ||
        any(nodes < 1) ||
        any(nodes > graph_size)
    ) {
      cli::cli_abort(
        "Floating part {.field nodes} must contain valid vertex indices."
      )
    }
    nodes <- as.integer(nodes)
    if (anyDuplicated(nodes) > 0) {
      cli::cli_abort("Floating part node indices must be unique.")
    }
    if (!root %in% nodes) {
      cli::cli_abort(
        "Floating part {.field nodes} must contain its {.field root}."
      )
    }
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
    root = root,
    nodes = nodes,
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

  parts <- purrr::map(
    parts,
    normalize_floating_part,
    graph_size = igraph::vcount(graph)
  )
  missing_nodes <- purrr::map_lgl(parts, ~ is.null(.x$nodes))
  if (any(missing_nodes)) {
    membership <- igraph::components(graph, mode = "weak")$membership
    parts[missing_nodes] <- purrr::map(
      parts[missing_nodes],
      function(part) {
        component <- membership[[part$root]]
        part$nodes <- as.integer(which(membership == component))
        part
      }
    )
  }

  parts
}

floating_main_vertices <- function(
  graph,
  parts = normalize_floating_parts(graph)
) {
  floating_nodes <- unlist(
    purrr::map(parts, "nodes"),
    use.names = FALSE
  )
  as.integer(setdiff(
    seq_len(igraph::vcount(graph)),
    floating_nodes
  ))
}

floating_part_candidate_parents <- function(graph, part) {
  if (length(part$parents) > 0) {
    return(as.integer(part$parents))
  }

  as.integer(setdiff(seq_len(igraph::vcount(graph)), part$nodes))
}

floating_part_membership <- function(graph, parts) {
  membership <- rep(NA_integer_, igraph::vcount(graph))
  for (part_id in seq_along(parts)) {
    membership[parts[[part_id]]$nodes] <- as.integer(part_id)
  }
  membership
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
  raw_parts <- igraph::graph_attr(graph, "floating_parts")
  if (is.null(raw_parts)) {
    if (!is_out_tree(graph)) {
      cli::cli_abort("Glycan structure must be an out tree.")
    }
    return(invisible(NULL))
  }

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

  for (part_id in seq_along(parts)) {
    component_nodes <- which(
      info$membership == info$floating_components[[part_id]]
    )
    if (!setequal(parts[[part_id]]$nodes, component_nodes)) {
      cli::cli_abort(c(
        "Floating part node metadata does not match its graph component.",
        "x" = "Floating part {part_id} must contain exactly node{?s} {.val {component_nodes}}."
      ))
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
        any(part$parents %in% part$nodes)
    ) {
      cli::cli_abort(
        "Floating part parent indices cannot refer to its own component."
      )
    }
  }

  invisible(NULL)
}

floating_linkage_acceptor_positions <- function(linkage) {
  acceptor <- sub(".*-", "", linkage)
  if (identical(acceptor, "?")) {
    return(integer())
  }

  unique(as.integer(strsplit(acceptor, "/", fixed = TRUE)[[1]]))
}

main_attachment_domains <- function(
  graph,
  vertices = seq_len(igraph::vcount(graph))
) {
  if (igraph::ecount(graph) == 0) {
    return(list())
  }

  endpoints <- igraph::as_edgelist(graph, names = FALSE)
  linkages <- igraph::edge_attr(graph, "linkage")
  included <- endpoints[, 1] %in% vertices & endpoints[, 2] %in% vertices

  lapply(which(included), function(edge_id) {
    positions <- floating_linkage_acceptor_positions(linkages[[edge_id]])
    list(
      known = length(positions) > 0,
      slots = paste(endpoints[edge_id, 1], positions, sep = "\r")
    )
  })
}

definitely_occupied_main_slots <- function(domains) {
  slot_sets <- lapply(domains, `[[`, "slots")
  unique(unlist(slot_sets[lengths(slot_sets) == 1], use.names = FALSE))
}

floating_attachment_domain <- function(
  part,
  part_id,
  candidate_parents,
  occupied_slots
) {
  positions <- floating_linkage_acceptor_positions(part$linkage)
  if (length(positions) == 0) {
    return(list(known = FALSE, slots = character()))
  }

  explicit <- length(part$parents) > 0
  parents <- candidate_parents
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
  slot_sets <- lapply(domains, `[[`, "slots")
  has_conflict_free_assignment(slot_sets)
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

  main_domains <- main_attachment_domains(graph)
  occupied_slots <- definitely_occupied_main_slots(main_domains)
  floating_domains <- lapply(seq_along(info$parts), function(part_id) {
    part <- info$parts[[part_id]]
    floating_attachment_domain(
      part,
      part_id,
      floating_part_candidate_parents(graph, part),
      occupied_slots
    )
  })
  domains <- c(main_domains, floating_domains)

  if (!has_conflict_free_attachment_assignment(domains)) {
    cli::cli_abort(c(
      "Floating parts cannot be attached simultaneously.",
      "x" = "No conflict-free assignment exists for the declared parent and acceptor positions."
    ))
  }

  invisible(NULL)
}

selected_floating_attachment_domain <- function(
  part,
  parent,
  occupied_slots
) {
  positions <- floating_linkage_acceptor_positions(part$linkage)
  if (length(positions) == 0) {
    return(list(known = FALSE, slots = character()))
  }

  list(
    known = TRUE,
    slots = setdiff(paste(parent, positions, sep = "\r"), occupied_slots)
  )
}

selected_floating_substituent_domain <- function(
  substituent,
  parent,
  occupied_slots
) {
  position <- substituent_position_tokens(substituent$substituent)
  if (identical(position, "?")) {
    return(list(known = FALSE, slots = character()))
  }
  positions <- stringr::str_split(position, stringr::fixed("/"))[[1]]

  list(
    known = TRUE,
    slots = setdiff(paste(parent, positions, sep = "\r"), occupied_slots)
  )
}

floating_part_assignment_is_acyclic <- function(
  parent_nodes,
  membership
) {
  parent_parts <- membership[parent_nodes]
  parent_parts[is.na(parent_parts)] <- 0L

  for (part_id in seq_along(parent_parts)) {
    current <- as.integer(part_id)
    visited <- integer()
    while (current != 0L) {
      if (current %in% visited) {
        return(FALSE)
      }
      visited <- c(visited, current)
      current <- parent_parts[[current]]
    }
  }

  TRUE
}

floating_parent_component_domains <- function(
  part_domains,
  membership
) {
  purrr::map(part_domains, function(parents) {
    components <- membership[parents]
    components[is.na(components)] <- 0L
    as.integer(components)
  })
}

floating_components_reach_main <- function(parent_component_domains) {
  if (length(parent_component_domains) == 0) {
    return(TRUE)
  }

  reaches_main <- rep(FALSE, length(parent_component_domains))
  repeat {
    newly_reachable <- vapply(
      seq_along(parent_component_domains),
      function(part_id) {
        if (reaches_main[[part_id]]) {
          return(FALSE)
        }
        parents <- parent_component_domains[[part_id]]
        floating_parents <- parents[parents > 0L]
        any(parents == 0L) || any(reaches_main[floating_parents])
      },
      logical(1)
    )
    if (!any(newly_reachable)) {
      break
    }
    reaches_main[newly_reachable] <- TRUE
  }

  all(reaches_main)
}

floating_part_assignment_closes_cycle <- function(
  part_id,
  selected_parent_components
) {
  visited <- rep(FALSE, length(selected_parent_components))
  current <- as.integer(part_id)
  while (current != 0L && !is.na(selected_parent_components[[current]])) {
    if (visited[[current]]) {
      return(TRUE)
    }
    visited[[current]] <- TRUE
    current <- selected_parent_components[[current]]
  }

  FALSE
}

has_valid_floating_assignment <- function(
  graph,
  parts,
  substituents
) {
  edge_domains <- main_attachment_domains(graph)
  substituent_domains <- main_substituent_domains(graph)
  fixed_domains <- c(edge_domains, substituent_domains)
  occupied_slots <- definitely_occupied_main_slots(fixed_domains)

  part_domains <- purrr::map(
    parts,
    ~ floating_part_candidate_parents(graph, .x)
  )
  floating_sub_domains <- purrr::map(
    substituents,
    ~ floating_substituent_candidate_parents(graph, .x)
  )

  for (part_id in seq_along(parts)) {
    floating_attachment_domain(
      parts[[part_id]],
      part_id,
      part_domains[[part_id]],
      occupied_slots
    )
  }
  for (substituent_id in seq_along(substituents)) {
    floating_substituent_domain(
      substituents[[substituent_id]],
      substituent_id,
      floating_sub_domains[[substituent_id]],
      occupied_slots
    )
  }

  membership <- floating_part_membership(graph, parts)
  parent_component_domains <- floating_parent_component_domains(
    part_domains,
    membership
  )
  if (!floating_components_reach_main(parent_component_domains)) {
    return(FALSE)
  }

  part_slot_constraints <- purrr::map_lgl(
    parts,
    ~ length(floating_linkage_acceptor_positions(.x$linkage)) > 0
  )
  substituent_slot_constraints <- purrr::map_lgl(
    substituents,
    ~ substituent_position_tokens(.x$substituent) != "?"
  )
  if (!any(part_slot_constraints) && !any(substituent_slot_constraints)) {
    return(has_conflict_free_attachment_assignment(fixed_domains))
  }

  selected_parts <- integer(length(parts))
  selected_parent_components <- rep(NA_integer_, length(parts))
  selected_substituents <- integer(length(substituents))

  assignment_has_free_slots <- function() {
    part_slot_domains <- purrr::map2(
      parts,
      selected_parts,
      ~ selected_floating_attachment_domain(.x, .y, occupied_slots)
    )
    substituent_slot_domains <- purrr::map2(
      substituents,
      selected_substituents,
      ~ selected_floating_substituent_domain(.x, .y, occupied_slots)
    )

    has_conflict_free_attachment_assignment(c(
      fixed_domains,
      part_slot_domains,
      substituent_slot_domains
    ))
  }

  assign_substituent <- function(substituent_id) {
    if (substituent_id > length(floating_sub_domains)) {
      return(assignment_has_free_slots())
    }
    for (parent in floating_sub_domains[[substituent_id]]) {
      selected_substituents[[substituent_id]] <<- parent
      if (assign_substituent(substituent_id + 1L)) {
        return(TRUE)
      }
    }
    FALSE
  }

  assign_part <- function(part_id) {
    if (part_id > length(part_domains)) {
      if (!floating_part_assignment_is_acyclic(selected_parts, membership)) {
        return(FALSE)
      }
      return(assign_substituent(1L))
    }
    for (candidate_id in seq_along(part_domains[[part_id]])) {
      parent <- part_domains[[part_id]][[candidate_id]]
      selected_parts[[part_id]] <<- parent
      selected_parent_components[[part_id]] <<-
        parent_component_domains[[part_id]][[candidate_id]]
      if (
        !floating_part_assignment_closes_cycle(
          part_id,
          selected_parent_components
        ) &&
          assign_part(part_id + 1L)
      ) {
        return(TRUE)
      }
    }
    selected_parts[[part_id]] <<- 0L
    selected_parent_components[[part_id]] <<- NA_integer_
    FALSE
  }

  assign_part(1L)
}

validate_floating_metadata_assignments <- function(graph) {
  if (
    !any(
      c("floating_parts", "floating_substituents") %in%
        igraph::graph_attr_names(graph)
    )
  ) {
    return(invisible(NULL))
  }

  parts <- normalize_floating_parts(graph)
  substituents <- normalize_floating_substituents(graph)
  if (length(parts) == 0 && length(substituents) == 0) {
    return(invisible(NULL))
  }

  if (!has_valid_floating_assignment(graph, parts, substituents)) {
    if (length(parts) > 0 && length(substituents) == 0) {
      cli::cli_abort(c(
        "Floating parts cannot be attached simultaneously.",
        "x" = "No conflict-free acyclic assignment connects every floating component to the main tree."
      ))
    }
    if (length(parts) == 0) {
      cli::cli_abort(c(
        "Floating substituents cannot be localized simultaneously.",
        "x" = "No conflict-free assignment exists for the declared parent residues and carbon positions."
      ))
    }
    cli::cli_abort(c(
      "Floating parts and substituents cannot be localized simultaneously.",
      "x" = "No conflict-free acyclic assignment exists for the declared parent domains."
    ))
  }

  invisible(NULL)
}

component_sequence_order <- function(
  graph,
  vertices,
  symmetry_labels = NULL
) {
  subgraph <- igraph::induced_subgraph(graph, vertices)
  subgraph <- delete_floating_parts_attr(subgraph)
  subgraph <- delete_floating_substituents_attr(subgraph)
  cache <- build_seq_cache(subgraph)
  order <- if (is.null(symmetry_labels)) {
    seq_glycan_order(cache$root, cache)
  } else {
    graph_names <- igraph::V(graph)$name
    subgraph_vertices <- match(igraph::V(subgraph)$name, graph_names)
    seq_glycan_order_with_floating_symmetry(
      cache$root,
      cache,
      symmetry_labels[subgraph_vertices]
    )
  }

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

attach_floating_parts <- function(
  graph,
  parts,
  part_ids,
  parent_nodes
) {
  for (selection_id in seq_along(part_ids)) {
    part_id <- part_ids[[selection_id]]
    part <- parts[[part_id]]
    graph <- igraph::add_edges(
      graph,
      c(parent_nodes[[selection_id]], part$root),
      linkage = part$linkage
    )
  }

  remaining <- parts[-part_ids]
  if (length(remaining) == 0) {
    return(delete_floating_parts_attr(graph))
  }

  membership <- igraph::components(graph, mode = "weak")$membership
  remaining_components <- membership[purrr::map_int(remaining, "root")]
  if (anyDuplicated(remaining_components) > 0) {
    cli::cli_abort(
      "Selected floating-part assignments merge unresolved component roots."
    )
  }

  remaining <- purrr::map2(
    remaining,
    remaining_components,
    function(part, component) {
      part$nodes <- as.integer(which(membership == component))
      if (length(part$parents) > 0) {
        part$parents <- setdiff(part$parents, part$nodes)
        if (length(part$parents) == 0) {
          cli::cli_abort(
            "Selected floating-part assignments leave an explicit parent domain empty."
          )
        }
      }
      part
    }
  )
  set_floating_parts_attr(graph, remaining)
}

resolve_single_parent_floating_parts <- function(graph) {
  repeat {
    parts <- normalize_floating_parts(graph)
    if (length(parts) == 0) {
      return(graph)
    }

    candidate_parents <- purrr::map(
      parts,
      ~ floating_part_candidate_parents(graph, .x)
    )
    resolved <- which(lengths(candidate_parents) == 1)
    if (length(resolved) == 0) {
      return(graph)
    }

    graph <- attach_floating_parts(
      graph,
      parts,
      resolved,
      purrr::map_int(candidate_parents[resolved], 1L)
    )
  }
}

canonicalize_floating_graph <- function(graph) {
  graph <- resolve_single_parent_floating_substituents(graph)
  graph <- resolve_single_parent_floating_parts(graph)
  if (!has_floating_metadata(graph)) {
    seq_cache <- build_seq_cache(graph)
    order <- seq_glycan_order(seq_cache$root, seq_cache)
    return(.reorder_by_sequence_order(graph, order))
  }

  parts <- normalize_floating_parts(graph)
  substituents <- normalize_floating_substituents(graph)
  main_vertices <- floating_metadata_main_vertices(graph, parts)
  labels <- floating_augmented_structure_labels(
    graph,
    parts,
    substituents,
    main_vertices
  )
  original_names <- igraph::V(graph)$name

  main_order <- component_sequence_order(
    graph,
    main_vertices,
    labels$vertices
  )

  floating_orders <- purrr::map2(
    parts,
    seq_along(parts),
    function(part, part_id) {
      order <- component_sequence_order(
        graph,
        part$nodes,
        labels$vertices
      )

      list(
        part_id = as.integer(part_id),
        order = order,
        linkage = part$linkage,
        parents = part$parents,
        key = paste(
          order$iupac,
          part$linkage,
          if (length(part$parents) == 0) "all" else "explicit",
          paste(
            sort(labels$vertices[
              floating_part_candidate_parents(graph, part)
            ]),
            collapse = ","
          ),
          sep = "\r"
        ),
        label = labels$parts[[part_id]]
      )
    }
  )

  floating_order <- order(
    purrr::map_chr(floating_orders, "key"),
    purrr::map_int(floating_orders, "label")
  )
  floating_orders <- floating_orders[floating_order]

  floating_substituent_orders <- purrr::map2(
    substituents,
    seq_along(substituents),
    function(substituent, substituent_id) {
      list(
        substituent = substituent$substituent,
        parents = substituent$parents,
        key = paste(
          substituent$substituent,
          if (length(substituent$parents) == 0) "all" else "explicit",
          paste(
            sort(labels$vertices[
              floating_substituent_candidate_parents(graph, substituent)
            ]),
            collapse = ","
          ),
          sep = "\r"
        ),
        label = labels$substituents[[substituent_id]]
      )
    }
  )
  floating_substituent_order <- order(
    purrr::map_chr(floating_substituent_orders, "key"),
    purrr::map_int(floating_substituent_orders, "label")
  )
  floating_substituent_orders <- floating_substituent_orders[
    floating_substituent_order
  ]

  floating_vertices <- unlist(
    purrr::map(floating_orders, c("order", "vertices")),
    use.names = FALSE
  )
  floating_edges <- unlist(
    purrr::map(floating_orders, c("order", "edges")),
    use.names = FALSE
  )
  vertex_order <- c(floating_vertices, main_order$vertices)
  edge_order <- c(floating_edges, main_order$edges)

  ordered_names <- original_names[vertex_order]
  canonical_index <- stats::setNames(seq_along(ordered_names), ordered_names)

  graph <- delete_floating_parts_attr(graph)
  graph <- delete_floating_substituents_attr(graph)
  graph <- .reorder_by_sequence_order(
    graph,
    list(vertices = vertex_order, edges = edge_order)
  )

  canonical_parts <- vector("list", length(floating_orders))
  for (i in seq_along(floating_orders)) {
    floating <- floating_orders[[i]]
    node_names <- original_names[floating$order$vertices]
    parent_names <- original_names[floating$parents]
    canonical_parts[[i]] <- list(
      root = as.integer(canonical_index[[floating$order$root_name]]),
      nodes = as.integer(unname(canonical_index[node_names])),
      linkage = floating$linkage,
      parents = sort(as.integer(unname(canonical_index[parent_names])))
    )
  }

  graph <- set_floating_parts_attr(graph, canonical_parts)
  canonical_substituents <- purrr::map(
    floating_substituent_orders,
    function(substituent) {
      list(
        substituent = substituent$substituent,
        parents = sort(as.integer(unname(canonical_index[
          original_names[substituent$parents]
        ])))
      )
    }
  )
  set_floating_substituents_attr(graph, canonical_substituents)
}

floating_part_iupac <- function(
  graph,
  part
) {
  subgraph <- igraph::induced_subgraph(graph, part$nodes)
  subgraph <- delete_floating_parts_attr(subgraph)
  subgraph <- delete_floating_substituents_attr(subgraph)
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
