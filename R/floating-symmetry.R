# Canonicalize complete floating forests and their parent constraints.

#' Order a Main Glycan Tree with Floating-Parent Symmetry
#'
#' Resolve otherwise indistinguishable main-tree branch orderings using every
#' component, floating-part parent relation, and floating substituent jointly.
#' This keeps complete-sequence parent indices stable across vertex, edge,
#' component, and metadata permutations.
#'
#' @param graph A validated glycan forest with vertex names.
#' @param info Complete floating metadata and main-tree vertex indices.
#'
#' @returns A sequence-order list in the same form as
#'   `component_sequence_order()`.
#' @noRd
floating_symmetry_main_order <- function(graph, info) {
  substituents <- info$substituents
  if (is.null(substituents)) {
    substituents <- list()
  }
  labels <- floating_augmented_structure_labels(
    graph,
    info$parts,
    substituents,
    info$main_vertices
  )
  component_sequence_order(
    graph,
    info$main_vertices,
    labels$vertices
  )
}

floating_augmented_structure_labels <- function(
  graph,
  parts,
  substituents,
  main_vertices
) {
  vertex_count <- igraph::vcount(graph)
  edge_count <- igraph::ecount(graph)
  part_count <- length(parts)
  substituent_count <- length(substituents)

  membership <- floating_part_membership(graph, parts)
  indegree <- igraph::degree(graph, mode = "in")
  part_roots <- purrr::map_int(parts, "root")
  roles <- rep("floating-node", vertex_count)
  roles[main_vertices] <- "main-node"
  roles[intersect(which(indegree == 0), main_vertices)] <- "main-root"
  roles[part_roots] <- "floating-root"
  color_keys <- paste(
    "residue",
    roles,
    igraph::V(graph)$mono,
    igraph::V(graph)$sub,
    sep = "\r"
  )

  edge_nodes <- vertex_count + seq_len(edge_count)
  if (edge_count > 0) {
    color_keys <- c(
      color_keys,
      paste("glycan-edge", igraph::E(graph)$linkage, sep = "\r")
    )
  }
  part_nodes <- vertex_count + edge_count + seq_len(part_count)
  if (part_count > 0) {
    color_keys <- c(
      color_keys,
      paste(
        "floating-part",
        purrr::map_chr(parts, "linkage"),
        sep = "\r"
      )
    )
  }
  substituent_nodes <-
    vertex_count + edge_count + part_count + seq_len(substituent_count)
  if (substituent_count > 0) {
    color_keys <- c(
      color_keys,
      paste(
        "floating-substituent",
        purrr::map_chr(substituents, "substituent"),
        sep = "\r"
      )
    )
  }

  augmented_edges <- integer()
  if (edge_count > 0) {
    endpoints <- igraph::as_edgelist(graph, names = FALSE)
    augmented_edges <- c(
      augmented_edges,
      as.integer(rbind(
        endpoints[, 1],
        edge_nodes,
        edge_nodes,
        endpoints[, 2]
      ))
    )
  }

  add_relation <- function(constraint, target, relation) {
    relation_node <- length(color_keys) + 1L
    color_keys <<- c(color_keys, relation)
    augmented_edges <<- c(
      augmented_edges,
      constraint,
      relation_node,
      relation_node,
      target
    )
  }

  for (part_id in seq_along(parts)) {
    part <- parts[[part_id]]
    add_relation(
      part_nodes[[part_id]],
      part$root,
      "floating-part-root"
    )
    candidates <- floating_part_candidate_parents(graph, part)
    for (parent in candidates) {
      add_relation(
        part_nodes[[part_id]],
        parent,
        "floating-part-parent"
      )
    }
  }
  for (substituent_id in seq_along(substituents)) {
    substituent <- substituents[[substituent_id]]
    candidates <- floating_substituent_candidate_parents(graph, substituent)
    for (parent in candidates) {
      add_relation(
        substituent_nodes[[substituent_id]],
        parent,
        "floating-substituent-parent"
      )
    }
  }

  augmented <- igraph::make_empty_graph(
    length(color_keys),
    directed = FALSE
  )
  if (length(augmented_edges) > 0) {
    augmented <- igraph::add_edges(augmented, augmented_edges)
  }
  colors <- match(
    color_keys,
    sort(unique(color_keys), method = "radix")
  )
  labels <- igraph::canonical_permutation(
    augmented,
    colors = colors
  )$labeling

  list(
    vertices = as.integer(labels[seq_len(vertex_count)]),
    parts = as.integer(labels[part_nodes]),
    substituents = as.integer(labels[substituent_nodes]),
    membership = membership
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
  substituents <- info$substituents
  if (is.null(substituents)) {
    substituents <- list()
  }
  explicit_substituents <- which(
    purrr::map_lgl(substituents, ~ length(.x$parents) > 0)
  )
  constraint_count <- length(explicit) + length(explicit_substituents)

  edge_nodes <- main_count + seq_len(edge_count)
  part_nodes <- main_count + edge_count + seq_len(constraint_count)
  augmented <- igraph::make_empty_graph(
    main_count + edge_count + constraint_count,
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
  offset <- length(explicit)
  for (i in seq_along(explicit_substituents)) {
    substituent <- substituents[[explicit_substituents[[i]]]]
    parent_names <- graph_names[substituent$parents]
    parents <- match(parent_names, main_names)
    constraint_node <- part_nodes[[offset + i]]
    augmented_edges <- c(
      augmented_edges,
      as.integer(rbind(rep(constraint_node, length(parents)), parents))
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
        floating_part_iupac(
          graph,
          part
        ),
        sep = "\r"
      )
    }
  )
  substituent_keys <- purrr::map_chr(
    explicit_substituents,
    ~ paste(
      "floating-substituent-parent-set",
      substituents[[.x]]$substituent,
      sep = "\r"
    )
  )
  color_keys <- c(main_keys, edge_keys, part_keys, substituent_keys)
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
