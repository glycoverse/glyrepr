# Internal helpers for floating glycan substructures.

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

  required <- c("root", "linkage")
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
  if (is.null(parents)) {
    parents <- integer()
  }
  if (
    !is.numeric(parents) ||
      anyNA(parents) ||
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

canonicalize_floating_graph <- function(graph) {
  parts <- normalize_floating_parts(graph)
  info <- floating_graph_info(graph, parts)
  original_names <- igraph::V(graph)$name

  main_order <- component_sequence_order(graph, info$main_vertices)
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
