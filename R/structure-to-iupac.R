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
#' @return A character vector representing the IUPAC sequences.
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
#' graph$alditol <- FALSE
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
    rlang::abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }

  unname(vctrs::field(glycan, "codes"))
}

# Internal function to convert a single igraph to IUPAC
.structure_to_iupac_single <- function(glycan) {
  # Find root (node with in-degree 0)
  root <- which(igraph::degree(glycan, mode = "in") == 0)
  if (length(root) != 1) {
    rlang::abort("Glycan structure must have exactly one root node.")
  }
  
  # Calculate depth for each node
  depths <- calculate_depths(glycan, root)
  
  # Generate sequence starting from root
  seq_result <- seq_glycan(glycan, root, depths)
  
  # Add root anomer at the end
  anomer <- glycan$anomer
  paste0(seq_result, "(", anomer, "-")
}

#' Parse linkage string into comparable components
#' 
#' @param linkage Character string in format "xy-z" (e.g., "b1-4", "a2-3")
#' @return Named list with x, y, z components and their numeric ranks
#' @keywords internal
parse_linkage <- function(linkage) {
  # Parse linkage format: xy-z where x is a/b/?, y is digit/?, z is digit/?
  pattern <- "^([ab\\?])(\\d|\\?)[-](\\d|\\?)$"
  match <- stringr::str_match(linkage, pattern)
  
  if (is.na(match[1])) {
    rlang::abort(glue::glue("Invalid linkage format: {linkage}"))
  }
  
  x <- match[2]  # anomeric configuration
  y <- match[3]  # first position  
  z <- match[4]  # second position
  
  # Convert to numeric ranks for comparison
  # For x: ? > b > a, so assign ? = 3, b = 2, a = 1
  x_rank <- switch(x, "a" = 1, "b" = 2, "?" = 3)
  
  # For y and z: ? is greater than any number
  y_rank <- if (y == "?") Inf else as.numeric(y)
  z_rank <- if (z == "?") Inf else as.numeric(z)
  
  list(x = x, y = y, z = z, x_rank = x_rank, y_rank = y_rank, z_rank = z_rank)
}

#' Compare two linkages
#' 
#' @param linkage1,linkage2 Character strings representing linkages
#' @return Integer: -1 if linkage1 < linkage2, 0 if equal, 1 if linkage1 > linkage2
#' @keywords internal
compare_linkages <- function(linkage1, linkage2) {
  p1 <- parse_linkage(linkage1)
  p2 <- parse_linkage(linkage2)
  
  # Compare x_rank first
  if (p1$x_rank != p2$x_rank) {
    return(sign(p1$x_rank - p2$x_rank))
  }
  
  # Then y_rank
  if (p1$y_rank != p2$y_rank) {
    return(sign(p1$y_rank - p2$y_rank))
  }
  
  # Finally z_rank
  sign(p1$z_rank - p2$z_rank)
}

#' Calculate depth (longest path to leaf) for each node
#' 
#' @param glycan An igraph object representing a glycan structure
#' @param root Root vertex
#' @return Named vector of depths for each node
#' @keywords internal
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

#' Generate glycan sequence recursively
#' 
#' @param glycan An igraph object representing a glycan structure
#' @param node Current node
#' @param depths Vector of depths for all nodes
#' @return Character string representing the sequence
#' @keywords internal
seq_glycan <- function(glycan, node, depths) {
  children <- igraph::neighbors(glycan, node, mode = "out")
  
  # Base case: leaf node
  if (length(children) == 0) {
    mono <- igraph::vertex_attr(glycan, "mono", node)
    sub <- igraph::vertex_attr(glycan, "sub", node)
    # Combine monosaccharide with substituent if present
    if (sub == "") {
      return(mono)
    } else {
      return(paste0(mono, sub))
    }
  }
  
  # Find backbone child (deepest, break ties by linkage)
  child_names <- as.character(children)
  child_depths <- depths[child_names]
  max_depth <- max(child_depths)
  backbone_candidates <- children[child_depths == max_depth]
  
  # Choose backbone child by linkage if multiple candidates
  if (length(backbone_candidates) > 1) {
    candidate_linkages <- character(length(backbone_candidates))
    for (i in seq_along(backbone_candidates)) {
      edge_id <- igraph::get_edge_ids(glycan, vp = c(as.character(node), as.character(backbone_candidates[i])))
      candidate_linkages[i] <- igraph::edge_attr(glycan, "linkage", edge_id)
    }
    
    # Sort by linkage (ascending order) - smaller linkage wins
    linkage_order <- order(candidate_linkages, decreasing = FALSE)
    backbone_child <- backbone_candidates[linkage_order[1]]
  } else {
    backbone_child <- backbone_candidates[1]
  }
  
  # Other children are branches
  branch_children <- children[!children %in% backbone_child]
  
  # Generate backbone sequence
  backbone_seq <- seq_glycan(glycan, backbone_child, depths)
  
  # Get backbone linkage
  backbone_edge_id <- igraph::get_edge_ids(glycan, vp = c(as.character(node), as.character(backbone_child)))
  backbone_linkage <- igraph::edge_attr(glycan, "linkage", backbone_edge_id)
  
  # Generate branch sequences
  branch_parts <- character()
  if (length(branch_children) > 0) {
    # Sort branches by linkage
    branch_linkages <- character(length(branch_children))
    for (i in seq_along(branch_children)) {
      edge_id <- igraph::get_edge_ids(glycan, vp = c(as.character(node), as.character(branch_children[i])))
      branch_linkages[i] <- igraph::edge_attr(glycan, "linkage", edge_id)
    }
    
    branch_order <- order(branch_linkages, decreasing = FALSE)
    branch_children <- branch_children[branch_order]
    
    for (i in seq_along(branch_children)) {
      branch_child <- branch_children[i]
      edge_id <- igraph::get_edge_ids(glycan, vp = c(as.character(node), as.character(branch_child)))
      branch_linkage <- igraph::edge_attr(glycan, "linkage", edge_id)
      
      # Generate branch sequence
      branch_seq <- seq_glycan(glycan, branch_child, depths)
      branch_parts[i] <- paste0("[", branch_seq, "(", branch_linkage, ")]")
    }
  }
  
  # Get current node's monosaccharide and substituent
  current_mono <- igraph::vertex_attr(glycan, "mono", node)
  current_sub <- igraph::vertex_attr(glycan, "sub", node)
  
  # Combine monosaccharide with substituent if present
  # Format: mono + substituent (e.g., "Glc" + "3Me" = "Glc3Me")
  current_mono_with_sub <- if (current_sub == "") {
    current_mono
  } else {
    paste0(current_mono, current_sub)
  }
  
  # Combine: backbone + (backbone_linkage) + branches + current_mono_with_sub
  branches_str <- paste0(branch_parts, collapse = "")
  paste0(backbone_seq, "(", backbone_linkage, ")", branches_str, current_mono_with_sub)
}

