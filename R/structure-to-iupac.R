#' Convert Glycan Structure to IUPAC-like Sequence
#'
#' @description
#' Convert a glycan structure to a sequence representation in the form of 
#' mono(linkage)mono, with branches represented by square brackets [].
#' The backbone is chosen as the longest path, and for branches, linkages are 
#' ordered lexicographically with smaller linkages on the backbone.
#'
#' @details
#' # Sequence Format
#' 
#' The sequence follows the format mono(linkage)mono, where:
#' - mono: monosaccharide name with optional substituents (e.g., Glc, GlcNAc, Glc3Me)
#' - linkage: glycosidic linkage (e.g., b1-4, a1-3)
#' - Branches are enclosed in square brackets []
#' - Substituents are appended directly to monosaccharide names (e.g., Glc3Me for Glc with 3Me substituent)
#' 
#' # Backbone Selection
#' 
#' The backbone is selected as the longest path in the tree. For branches,
#' the same rule applies recursively.
#' 
#' # Linkage Comparison
#' 
#' Linkages are compared lexicographically:
#' 1. First by anomeric configuration: ? > b > a
#' 2. Then by first position: ? > numbers (numerically)
#' 3. Finally by second position: ? > numbers (numerically)
#' 
#' Smaller linkages are placed on the backbone, larger ones in branches.
#'
#' @param glycan A glyrepr_structure vector.
#'
#' @returns A character vector representing the IUPAC sequences.
#'
#' @examples
#' # Simple linear structure
#' structure_to_iupac(o_glycan_core_1())
#' 
#' # Branched structure  
#' structure_to_iupac(n_glycan_core())
#' 
#' # Structure with substituents
#' graph <- igraph::make_graph(~ 1-+2)
#' igraph::V(graph)$mono <- c("Glc", "GlcNAc")
#' igraph::V(graph)$sub <- c("3Me", "6Ac")
#' igraph::E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' glycan <- glycan_structure(graph)
#' structure_to_iupac(glycan)  # Returns "GlcNAc6Ac(b1-4)Glc3Me(a1-"
#' 
#' # Vectorized structures
#' structs <- glycan_structure(o_glycan_core_1(), n_glycan_core())
#' structure_to_iupac(structs)
#'
#' @export
structure_to_iupac <- function(glycan) {
  if (!is_glycan_structure(glycan)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }

  data <- vctrs::vec_data(glycan)
  unname(vctrs::field(data, "iupac"))
}

# Internal function to convert a single igraph to IUPAC
.structure_to_iupac_single <- function(glycan) {
  root <- which(igraph::degree(glycan, mode = "in") == 0)
  seq_cache <- build_seq_cache(glycan, root)

  # Step 1: Generate pseudo-IUPAC sequence starting from root
  # `pseudo_seq` is a string like "V1E1[V2E2]V3E3".
  # V1, V2, V3 are the vertex indices, E1, E2, E3 are the edge indices.
  pseudo_seq <- seq_glycan(root, seq_cache)

  # Step 2: Replace vertex and edge indices with actual monosaccharides and linkages
  real_seq <- replace_mono_and_link(pseudo_seq, glycan)

  anomer <- glycan$anomer
  paste0(real_seq, "(", anomer, "-")
}

#' Calculate linkage rank used for ordering linkages
#'
#' This is the reciprocal of the second position of the linkage.
#' "?" always rank the highest, so assigned 1.
#' For example, "b1-4" has rank 1/4, "a2-3" has rank 1/3, "a?-?" has rank 1.
#'
#' @param linkages Character vector in format "xy-z" (e.g., "b1-4", "a2-3")
#' @returns Numeric rank of the linkage.
#' @noRd
calculate_linkage_rank <- function(linkages) {
  pos2 <- stringr::str_split_i(linkages, "-", 2)  # second position
  suppressWarnings(
    dplyr::if_else(
      pos2 == "?" | stringr::str_detect(pos2, "/"),
      1,
      1 / as.numeric(pos2)
    )
  )
}

#' Calculate depth (longest path to leaf) for each node
#'
#' @param glycan An igraph object representing a glycan structure
#' @param root Root vertex
#' @returns Numeric vector of depths for each node (indexed by vertex id)
#' @noRd
calculate_depths <- function(glycan, root) {
  vcount <- igraph::vcount(glycan)
  depths <- rep(NA_real_, vcount)

  # Use a local function to calculate depths with memoization
  calculate_single_depth <- function(node) {
    node_index <- as.integer(node)

    if (!is.na(depths[node_index])) {
      return(depths[node_index])
    }

    # Get child nodes
    children <- as.integer(igraph::neighbors(glycan, node_index, mode = "out"))

    if (length(children) == 0) {
      # Leaf node
      depths[node_index] <<- 0
    } else {
      # Calculate depth for all children first
      child_depths <- vapply(children, calculate_single_depth, numeric(1))

      # This node's depth is 1 + max child depth
      depths[node_index] <<- 1 + max(child_depths)
    }

    depths[node_index]
  }

  # Calculate depths for all nodes starting from root
  calculate_single_depth(root)

  # Fill in any remaining nodes (shouldn't be needed for connected graph)
  for (v in seq_len(vcount)) {
    if (is.na(depths[v])) {
      calculate_single_depth(v)
    }
  }

  depths
}

#' Build adjacency cache for sequence generation
#'
#' It caches the following information:
#' - children: children[[i]] is all the children ids of vertex i.
#' - edge_ids: edge_ids[[i]] is all the edge ids directed from vertex i.
#' - linkages: linkages[[i]] is all the linkage attributes (e.g. "b1-4") of the edges directed from vertex i.
#' - depths: depths[i] is the maximum depth of the subtree rooted at vertex i.
#' - signatures: signatures[i] is the signature of the subtree rooted at vertex i.
#'   The signature of a node is a string used for branch ordering in ties breaking.
#'
#' @param glycan An igraph object representing a glycan structure
#' @param root Root vertex index used for depth calculation
#' @returns List containing child vertices, edge ids, linkages per parent, and node depths
#' @noRd
build_seq_cache <- function(glycan, root) {
  vcount <- igraph::vcount(glycan)
  edge_ids <- seq_len(igraph::ecount(glycan))
  edge_vertices <- igraph::ends(glycan, igraph::E(glycan), names = FALSE)
  edge_linkages <- igraph::edge_attr(glycan, "linkage")

  children <- vector("list", vcount)
  parent_edge_ids <- vector("list", vcount)
  parent_linkages <- vector("list", vcount)

  for (edge_id in edge_ids) {
    parent <- edge_vertices[edge_id, 1]
    child <- edge_vertices[edge_id, 2]

    children[[parent]] <- c(children[[parent]], child)
    parent_edge_ids[[parent]] <- c(parent_edge_ids[[parent]], edge_id)
    parent_linkages[[parent]] <- c(parent_linkages[[parent]], edge_linkages[edge_id])
  }

  # ===== Calculate node signatures =====
  mono_vec <- igraph::vertex_attr(glycan, "mono")
  sub_vec  <- igraph::vertex_attr(glycan, "sub")
  mono_sub <- ifelse(
    is.na(sub_vec) | sub_vec == "",
    mono_vec,
    paste0(mono_vec, stringr::str_remove_all(sub_vec, ","))
  )

  signature <- rep(NA_character_, vcount)

  compute_sig <- function(node) {
    if (!is.na(signature[node])) {
      return(signature[node])
    }

    kids <- children[[node]]

    if (is.null(kids) || length(kids) == 0) {
      signature[node] <<- mono_sub[node]
    } else {
      toks <- character(length(kids))
      for (i in seq_along(kids)) {
        kid <- kids[i]
        link_i <- parent_linkages[[node]][i]
        toks[i] <- paste0(link_i, "->", compute_sig(kid))
      }

      toks <- sort(toks)

      signature[node] <<- paste0(
        mono_sub[node],
        "{",
        paste0(toks, collapse = ","),
        "}"
      )
    }

    signature[node]
  }

  # trigger calculation of all node signatures
  for (n in seq_len(vcount)) {
    compute_sig(n)
  }
  # ===== End of calculating node signatures =====

  list(
    children = children,
    edge_ids = parent_edge_ids,
    linkages = parent_linkages,
    depths = calculate_depths(glycan, root),
    signatures = signature
  )
}

#' Generate glycan sequence recursively using pseudo-IUPAC format
#'
#' @param node Current node
#' @param cache Precomputed adjacency and edge metadata
#' @returns Character string representing the pseudo-IUPAC sequence
#' @noRd
seq_glycan <- function(node, cache) {
  children <- cache$children[[node]]
  edge_ids <- cache$edge_ids[[node]]

  # Base case: leaf node
  if (is.null(children) || length(children) == 0) {
    # Return vertex index with V prefix
    return(paste0("V", as.character(node)))
  }

  if (length(children) > 1) {
    children_order <- order_branches(node, cache)
    backbone_child <- children[[children_order$backbone]]
    backbone_edge_id <- edge_ids[[children_order$backbone]]
    backbone_seq <- seq_glycan(backbone_child, cache)
    branch_children <- children[children_order$branches]
    branch_edge_ids <- edge_ids[children_order$branches]
    branch_seqs <- purrr::map_chr(branch_children, ~ seq_glycan(.x, cache))
    branch_seqs <- paste0("[", branch_seqs, "E", branch_edge_ids, "]")
  } else {
    backbone_child <- children[[1]]
    backbone_edge_id <- edge_ids[[1]]
    backbone_seq <- seq_glycan(backbone_child, cache)
    branch_seqs <- ""
  }

  # Combine: backbone + E<edge_index> + branches + V<node_index>
  branches_str <- paste0(branch_seqs, collapse = "")
  paste0(backbone_seq, "E", backbone_edge_id, branches_str, "V", as.character(node))
}

#' Order branches
#'
#' This function orders all branches of a node by depth, then by linkages.
#'
#' @param node A node index.
#' @param cache Precomputed adjacency and edge metadata
#' @returns A list of two elements:
#'   - "backbone": the index of the backbone child.
#'   - "branches": the indices of the branches in order.
#'   All indices can only be used on `cache$children[[node]]`.
#'   For example, when `backbone` is 2, it means `cache$children[[node]][[2]]` is the backbone child.
#' @noRd
order_branches <- function(node, cache) {
  # Find the backbone child
  children <- cache$children[[node]]
  child_depths <- cache$depths[children]
  child_linkages <- cache$linkages[[node]]
  linkage_ranks <- calculate_linkage_rank(child_linkages)
  child_sigs <- cache$signatures[children]
  backbone_child_index <- order(child_depths, linkage_ranks, child_sigs, decreasing = TRUE)[[1]]

  # Order the rest of the children
  if (length(children) > 1) {
    branch_order <- order(linkage_ranks, child_sigs, decreasing = TRUE)
    branch_order <- branch_order[branch_order != backbone_child_index]
    list(backbone = backbone_child_index, branches = branch_order)
  } else {
    list(backbone = backbone_child_index, branches = integer(0))
  }
}

#' Replace vertex and edge indices with actual monosaccharides and linkages
#'
#' @param pseudo_seq Character string containing pseudo-IUPAC sequence with V and E prefixes
#' @param glycan An igraph object representing a glycan structure
#' @returns Character string with V<index> replaced by monosaccharides and E<index> replaced by linkages
#' @noRd
replace_mono_and_link <- function(pseudo_seq, glycan) {
  # Replace vertex indices (V<index>) with monosaccharides and substituents
  vertex_pattern <- "V(\\d+)"
  vertex_matches <- stringr::str_match_all(pseudo_seq, vertex_pattern)[[1]]

  if (nrow(vertex_matches) > 0) {
    for (i in seq_len(nrow(vertex_matches))) {
      vertex_index <- as.numeric(vertex_matches[i, 2])
      mono <- igraph::vertex_attr(glycan, "mono", vertex_index)
      sub <- igraph::vertex_attr(glycan, "sub", vertex_index)

      # Combine monosaccharide with substituent if present
      mono_with_sub <- if (sub == "") {
        mono
      } else {
        # Remove commas to get IUPAC format
        iupac_sub <- stringr::str_remove_all(sub, ",")
        paste0(mono, iupac_sub)
      }

      # Replace V<index> with actual monosaccharide
      pseudo_seq <- stringr::str_replace(
        pseudo_seq, paste0("V", vertex_index), mono_with_sub
      )
    }
  }

  # Replace edge indices (E<index>) with linkages
  edge_pattern <- "E(\\d+)"
  edge_matches <- stringr::str_match_all(pseudo_seq, edge_pattern)[[1]]

  if (nrow(edge_matches) > 0) {
    for (i in seq_len(nrow(edge_matches))) {
      edge_index <- as.numeric(edge_matches[i, 2])
      linkage <- igraph::edge_attr(glycan, "linkage", edge_index)
      # Replace E<index> with actual linkage wrapped in parentheses
      pseudo_seq <- stringr::str_replace(
        pseudo_seq, paste0("E", edge_index), paste0("(", linkage, ")")
      )
    }
  }

  pseudo_seq
}
