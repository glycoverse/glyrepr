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
  depths <- calculate_depths(glycan, root)
  seq_cache <- build_seq_cache(glycan)

  # Step 1: Generate pseudo-IUPAC sequence starting from root
  # `pseudo_seq` is a string like "V1E1[V2E2]V3E3".
  # V1, V2, V3 are the vertex indices, E1, E2, E3 are the edge indices.
  pseudo_seq <- seq_glycan(root, depths, seq_cache)

  # Step 2: Replace vertex and edge indices with actual monosaccharides and linkages
  real_seq <- replace_mono_and_link(pseudo_seq, glycan)

  anomer <- glycan$anomer
  paste0(real_seq, "(", anomer, "-")
}

#' Parse linkage string into comparable components
#'
#' @param linkage Character string in format "xy-z" (e.g., "b1-4", "a2-3")
#' @returns Named list with x, y, z components and their numeric ranks
#' @noRd
parse_linkage <- function(linkage) {
  # Parse linkage format: xy-z where x is a/b/?, y is digit/?, z is digit/?

  x <- stringr::str_sub(linkage, 1, 1)  # anomeric configuration
  link_part <- stringr::str_sub(linkage, 2, -1)  # first position  
  y <- stringr::str_split_i(link_part, "-", 1)  # first position  
  z <- stringr::str_split_i(link_part, "-", 2)  # second position

  # Convert to numeric ranks for comparison
  # For x: ? > b > a, so assign ? = 3, b = 2, a = 1
  x_rank <- switch(x, "a" = 1, "b" = 2, "?" = 3)

  # For y and z: ? is greater than any number
  y_rank <- if (y == "?" || stringr::str_detect(y, "/")) 0 else as.numeric(y)
  z_rank <- if (z == "?" || stringr::str_detect(z, "/")) 0 else as.numeric(z)

  list(x = x, y = y, z = z, x_rank = x_rank, y_rank = y_rank, z_rank = z_rank)
}

#' Order linkages
#'
#' @param linkages A character vector of linkages.
#' @param decreasing Logical. If TRUE, the linkages are ordered in decreasing order.
#' @returns A character vector of ordered linkages.
#' @noRd
order_linkages <- function(linkages, decreasing = FALSE) {
  parsed <- purrr::map(linkages, parse_linkage)
  rank <- purrr::map_dbl(parsed, ~ .x$z_rank)
  order(rank, decreasing = decreasing)
}

#' Calculate depth (longest path to leaf) for each node
#' 
#' @param glycan An igraph object representing a glycan structure
#' @param root Root vertex
#' @returns Named vector of depths for each node
#' @noRd
calculate_depths <- function(glycan, root) {
  all_vertices <- igraph::V(glycan)
  depths <- rep(NA_real_, length(all_vertices))
  names(depths) <- names(all_vertices)

  # Use a local function to calculate depths with memoization
  calculate_single_depth <- function(node) {
    node_name <- as.character(node)

    if (!is.na(depths[node_name])) {
      return(depths[node_name])
    }

    # Get child nodes
    children <- igraph::neighbors(glycan, node, mode = "out")

    if (length(children) == 0) {
      # Leaf node
      depths[node_name] <<- 0
    } else {
      # Calculate depth for all children first
      child_depths <- sapply(children, calculate_single_depth)

      # This node's depth is 1 + max child depth
      depths[node_name] <<- 1 + max(child_depths)
    }

    depths[node_name]
  }

  # Calculate depths for all nodes starting from root
  calculate_single_depth(root)

  # Fill in any remaining nodes (shouldn't be needed for connected graph)
  for (v in all_vertices) {
    if (is.na(depths[as.character(v)])) {
      calculate_single_depth(v)
    }
  }

  depths
}

#' Build adjacency cache for sequence generation
#'
#' @param glycan An igraph object representing a glycan structure
#' @returns List containing child vertices, edge ids, and linkages per parent
#' @noRd
build_seq_cache <- function(glycan) {
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

  list(
    children = children,
    edge_ids = parent_edge_ids,
    linkages = parent_linkages
  )
}

#' Generate glycan sequence recursively using pseudo-IUPAC format
#'
#' @param node Current node
#' @param depths Vector of depths for all nodes
#' @param cache Precomputed adjacency and edge metadata
#' @returns Character string representing the pseudo-IUPAC sequence
#' @noRd
seq_glycan <- function(node, depths, cache) {
  children <- cache$children[[node]]

  # Base case: leaf node
  if (is.null(children) || length(children) == 0) {
    # Return vertex index with V prefix
    return(paste0("V", as.character(node)))
  }

  # Find backbone child (deepest, break ties by linkage)
  child_names <- as.character(children)
  child_depths <- depths[child_names]
  max_depth <- max(child_depths)
  backbone_candidate_idx <- which(child_depths == max_depth)
  child_linkages <- cache$linkages[[node]]
  child_edge_ids <- cache$edge_ids[[node]]

  # Choose backbone child by linkage if multiple candidates
  if (length(backbone_candidate_idx) > 1) {
    candidate_linkages <- child_linkages[backbone_candidate_idx]
    # Sort by linkage (ascending order) - smaller linkage wins
    linkage_order <- order_linkages(candidate_linkages, decreasing = FALSE)
    backbone_idx <- backbone_candidate_idx[linkage_order[1]]
  } else {
    backbone_idx <- backbone_candidate_idx[1]
  }

  backbone_child <- children[backbone_idx]
  # Other children are branches
  branch_idx <- setdiff(seq_along(children), backbone_idx)

  # Generate backbone sequence
  backbone_seq <- seq_glycan(backbone_child, depths, cache)

  # Get backbone edge index
  backbone_edge_id <- child_edge_ids[backbone_idx]

  # Generate branch sequences
  branch_parts <- character()
  if (length(branch_idx) > 0) {
    # Sort branches by linkage
    branch_linkages <- child_linkages[branch_idx]
    branch_order <- order_linkages(branch_linkages, decreasing = FALSE)
    branch_idx <- branch_idx[branch_order]

    for (i in seq_along(branch_idx)) {
      idx <- branch_idx[i]
      branch_child <- children[idx]
      edge_id <- child_edge_ids[idx]

      # Generate branch sequence with edge index
      branch_seq <- seq_glycan(branch_child, depths, cache)
      branch_parts[i] <- paste0("[", branch_seq, "E", edge_id, "]")
    }
  }

  # Combine: backbone + E<edge_index> + branches + V<node_index>
  branches_str <- paste0(branch_parts, collapse = "")
  paste0(backbone_seq, "E", backbone_edge_id, branches_str, "V", as.character(node))
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
