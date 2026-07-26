# Canonicalize main-tree ties using explicit floating-part parent sets.

#' Order a Main Glycan Tree with Floating-Parent Symmetry
#'
#' Resolve otherwise indistinguishable main-tree branch orderings using all
#' explicit floating-part parent sets jointly. This keeps parent indices stable
#' across vertex and edge permutations without collapsing relationships between
#' multiple floating parts.
#'
#' @param graph A validated glycan forest with vertex names.
#' @param info Floating-graph component information from
#'   [floating_graph_info()].
#'
#' @returns A sequence-order list in the same form as
#'   [component_sequence_order()].
#' @noRd
floating_symmetry_main_order <- function(graph, info) {
  has_explicit_parents <- purrr::some(
    info$parts,
    ~ length(.x$parents) > 0
  )
  if (!has_explicit_parents) {
    return(component_sequence_order(graph, info$main_vertices))
  }

  main <- igraph::induced_subgraph(graph, info$main_vertices)
  main <- delete_floating_parts_attr(main)
  symmetry_labels <- floating_augmented_main_labels(graph, main, info)

  cache <- build_seq_cache(main)
  order <- seq_glycan_order_with_floating_symmetry(
    cache$root,
    cache,
    symmetry_labels
  )

  vertex_names <- igraph::V(main)$name[order$vertices]
  edge_table <- igraph::as_data_frame(main, what = "edges")
  edge_keys <- paste(edge_table$from, edge_table$to, sep = "\r")
  ordered_edge_keys <- edge_keys[order$edges]

  graph_edges <- igraph::as_data_frame(graph, what = "edges")
  graph_edge_keys <- paste(graph_edges$from, graph_edges$to, sep = "\r")

  list(
    vertices = match(vertex_names, igraph::V(graph)$name),
    edges = match(ordered_edge_keys, graph_edge_keys),
    iupac = seq_glycan_iupac(cache$root, cache),
    root_name = igraph::V(main)$name[cache$root]
  )
}


#' Label Main Vertices in an Augmented Constraint Graph
#'
#' Build an undirected colored graph containing the rooted, edge-labeled main
#' tree and one constraint vertex per explicit floating parent set. BLISS
#' canonical labels then provide an isomorphism-invariant final tie-breaker for
#' otherwise identical glycan branches.
#'
#' @param graph The complete glycan forest.
#' @param main The induced main-tree graph.
#' @param info Floating-graph component information.
#'
#' @returns An integer canonical label for each vertex of `main`.
#' @noRd
floating_augmented_main_labels <- function(graph, main, info) {
  main_count <- igraph::vcount(main)
  main_edges <- igraph::as_edgelist(main, names = FALSE)
  edge_count <- nrow(main_edges)
  explicit <- which(
    purrr::map_lgl(info$parts, ~ length(.x$parents) > 0)
  )
  part_count <- length(explicit)

  edge_nodes <- main_count + seq_len(edge_count)
  part_nodes <- main_count + edge_count + seq_len(part_count)
  augmented <- igraph::make_empty_graph(
    main_count + edge_count + part_count,
    directed = FALSE
  )

  augmented_edges <- integer()
  if (edge_count > 0) {
    augmented_edges <- c(
      augmented_edges,
      as.integer(rbind(
        main_edges[, 1],
        edge_nodes,
        edge_nodes,
        main_edges[, 2]
      ))
    )
  }

  main_names <- igraph::V(main)$name
  graph_names <- igraph::V(graph)$name
  for (i in seq_along(explicit)) {
    part <- info$parts[[explicit[[i]]]]
    parent_names <- graph_names[part$parents]
    parents <- match(parent_names, main_names)
    augmented_edges <- c(
      augmented_edges,
      as.integer(rbind(rep(part_nodes[[i]], length(parents)), parents))
    )
  }
  if (length(augmented_edges) > 0) {
    augmented <- igraph::add_edges(augmented, augmented_edges)
  }

  cache <- build_seq_cache(main)
  main_keys <- paste(
    "main",
    igraph::V(main)$mono,
    igraph::V(main)$sub,
    ifelse(seq_len(main_count) == cache$root, "root", "node"),
    sep = "\r"
  )
  edge_keys <- if (edge_count == 0) {
    character()
  } else {
    paste(
      "main-edge",
      igraph::E(main)$linkage,
      sep = "\r"
    )
  }
  part_keys <- purrr::map_chr(
    explicit,
    function(i) {
      part <- info$parts[[i]]
      part$parents <- integer()
      paste(
        "floating-parent-set",
        floating_part_iupac(graph, part, info$membership),
        sep = "\r"
      )
    }
  )
  color_keys <- c(main_keys, edge_keys, part_keys)
  colors <- match(
    color_keys,
    sort(unique(color_keys), method = "radix")
  )

  labels <- igraph::canonical_permutation(
    augmented,
    colors = colors
  )$labeling
  as.integer(labels[seq_len(main_count)])
}


#' Generate IUPAC Vertex and Edge Order with a Symmetry Tie-Breaker
#'
#' @param node Current main-tree node.
#' @param cache Precomputed sequence cache.
#' @param symmetry_labels Canonical augmented-graph labels.
#'
#' @returns A list with integer `vertices` and `edges` vectors.
#' @noRd
seq_glycan_order_with_floating_symmetry <- function(
  node,
  cache,
  symmetry_labels
) {
  children <- cache$children[[node]]
  if (length(children) == 0) {
    return(list(vertices = node, edges = integer()))
  }

  if (length(children) > 1) {
    children_order <- order_branches_with_floating_symmetry(
      node,
      cache,
      symmetry_labels
    )
    backbone_child <- children[[children_order$backbone]]
    backbone_edge <- cache$edge_ids[[node]][[children_order$backbone]]
    backbone_order <- seq_glycan_order_with_floating_symmetry(
      backbone_child,
      cache,
      symmetry_labels
    )
    branch_orders <- lapply(
      children_order$branches,
      function(branch_index) {
        branch_child <- children[[branch_index]]
        branch_order <- seq_glycan_order_with_floating_symmetry(
          branch_child,
          cache,
          symmetry_labels
        )
        list(
          vertices = branch_order$vertices,
          edges = c(
            branch_order$edges,
            cache$edge_ids[[node]][[branch_index]]
          )
        )
      }
    )
  } else {
    backbone_child <- children[[1]]
    backbone_edge <- cache$edge_ids[[node]][[1]]
    backbone_order <- seq_glycan_order_with_floating_symmetry(
      backbone_child,
      cache,
      symmetry_labels
    )
    branch_orders <- list()
  }

  list(
    vertices = c(
      backbone_order$vertices,
      unlist(lapply(branch_orders, `[[`, "vertices"), use.names = FALSE),
      node
    ),
    edges = c(
      backbone_order$edges,
      backbone_edge,
      unlist(lapply(branch_orders, `[[`, "edges"), use.names = FALSE)
    )
  )
}


#' Order Structurally Tied Branches with Floating-Parent Constraints
#'
#' @param node Current main-tree node.
#' @param cache Precomputed sequence cache.
#' @param symmetry_labels Canonical augmented-graph labels.
#'
#' @returns Backbone and branch indices into `cache$children[[node]]`.
#' @noRd
order_branches_with_floating_symmetry <- function(
  node,
  cache,
  symmetry_labels
) {
  children <- cache$children[[node]]
  child_depths <- cache$depths[children]
  linkage_ranks <- cache$linkage_ranks[[node]]
  child_signatures <- cache$signatures[children]
  child_symmetry <- symmetry_labels[children]

  backbone <- order(
    child_depths,
    linkage_ranks,
    child_signatures,
    child_symmetry,
    decreasing = TRUE
  )[[1]]

  branches <- order(
    linkage_ranks,
    child_signatures,
    child_symmetry,
    decreasing = TRUE
  )
  branches <- branches[branches != backbone]

  list(backbone = backbone, branches = branches)
}
