#' Create a Glycan Structure Vector
#'
#' A vector of glycan structures with efficient storage using hash-based deduplication.
#'
#' @details
#' The underlying implementation uses hash values of IUPAC codes of the glycan structures.
#' This prevents redundant storage and computation, which is very useful for glycan structures.
#' 
#' Each glycan structure is a directed graph modeling of a glycan structure,
#' where nodes are monosaccharides and edges are linkages.
#'
#' # S3 Classes
#'
#' A glycan structure vector is a vctrs record with additional S3 class `glyrepr_structure`.
#' Therefore, `sloop::s3_class()` of a glycan structure vector is 
#' `c("glyrepr_structure", "vctrs_rcrd")`.
#'
#' Constraints for individual structures:
#' - The graph must be directed and an outward tree (reducing end as root).
#' - The graph must have a graph attribute `anomer`, in the form of "a1".
#' - The graph must have a graph attribute `alditol`, a logical value.
#' - The graph must have a vertex attribute `mono` for monosaccharide names.
#' - The graph must have a vertex attribute `sub` for substituents.
#' - The graph must have an edge attribute `linkage` for linkages.
#' - Monosaccharide name must be known, either generic (Hex, HexNAc, etc.)
#'   or concrete (Glc, Gal, etc.), but not a mixture of both.
#'   NA is not allowed.
#' - Substituent must be valid, in the form of "xY", where x is position
#'   and Y is substituent name, e.g. "2Ac", "3S", etc.
#' - Linkages must be valid, in the form of "a/bX-Y", where X and Y are integers,
#'   e.g. "b1-4", "a2-3", etc.
#'   NA is not allowed.
#'
#' @param ... igraph graph objects to be converted to glycan structures, or existing glycan structure vectors.
#' @param x An object to check or convert.
#'
#' @return A glycan structure vector object.
#'
#' @examples
#' library(igraph)
#'
#' # A simple glycan structure: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2)  # 1 and 2 are monosaccharides
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")
#' V(graph)$sub <- ""
#' E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' 
#' # Create glycan structure vector
#' glycan_structure(graph)
#' 
#' # Create vector with multiple structures
#' glycan_structure(graph, o_glycan_core_1())
#'
#' @importFrom magrittr %>%
#' @export
glycan_structure <- function(...) {
  args <- list(...)
  
  # Handle different input types
  graphs <- list()
  for (arg in args) {
    if (is_glycan_structure(arg)) {
      # Extract individual structures from existing glycan_structure vector
      graphs <- c(graphs, get_structures_from_vector(arg))
    } else if (inherits(arg, "igraph")) {
      graphs <- c(graphs, list(arg))
    } else {
      rlang::abort("All arguments must be igraph objects or glycan_structure vectors.")
    }
  }
  
  if (length(graphs) == 0) {
    # Return empty vector
    return(new_glycan_structure(character(), character()))
  }
  
  # Validate and process each graph
  processed_graphs <- purrr::map(graphs, function(graph) {
    checkmate::assert_class(graph, "igraph")
    graph %>%
      validate_single_glycan_structure() %>%
      ensure_name_vertex_attr()
  })
  
  # Use IUPAC codes directly as data for the rcrd structure
  iupacs <- purrr::map_chr(processed_graphs, .structure_to_iupac_single)
  
  # Get mono types for each structure
  mono_types <- purrr::map_chr(processed_graphs, get_graph_mono_type)
  
  # Create a unique list based on uniqueness of IUPAC codes for structures storage
  unique_indices <- which(!duplicated(iupacs))
  unique_graphs <- processed_graphs[unique_indices]
  unique_iupacs <- iupacs[unique_indices]
  names(unique_graphs) <- unique_iupacs
  
  new_glycan_structure(iupacs, mono_types, unique_graphs)
}

# Helper function to validate a single glycan structure
validate_single_glycan_structure <- function(glycan) {
  checkmate::assert_class(glycan, "igraph")

  # Check if it is a directed graph
  if (!is_directed_graph(glycan)) {
    rlang::abort("Glycan structure must be directed.")
  }
  
  # Check if it is an out tree
  if (!is_out_tree(glycan)) {
    rlang::abort("Glycan structure must be an out tree.")
  }
  
  # Check if graph has a vertex attribute "mono"
  # This is the monosaccharide name, e.g. "GlcNAc", "Man", etc.
  if (!has_vertex_attrs(glycan, "mono")) {
    rlang::abort("Glycan structure must have a vertex attribute 'mono'.")
  }
  
  # Check if no NA in vertex attribute "mono"
  mono_names <- igraph::vertex_attr(glycan, "mono")
  if (any(is.na(mono_names))) {
    rlang::abort("Glycan structure must have no NA in vertex attribute 'mono'.")
  }
  
  # Check if all monosaccharides are known
  if (!all(is_known_mono(mono_names))) {
    unknown_monos <- unique(igraph::V(glycan)$mono[!is_known_mono(igraph::V(glycan)$mono)])
    msg <- glue::glue("Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}")
    rlang::abort(msg, monos = unknown_monos)
  }
  
  # Check if mixed use of generic and concrete monosaccharides
  if (mix_generic_concrete(mono_names)) {
    rlang::abort("Monosaccharides must be either all generic or all concrete.")
  }
  
  # Check if graph has a vertex attribute "sub"
  # This is the substituent name, e.g. "Ac", "S", "P", or "" (no).
  if (!has_vertex_attrs(glycan, "sub")) {
    rlang::abort("Glycan structure must have a vertex attribute 'sub'.")
  }
  
  # Check if no NA in vertex attribute "sub"
  subs <- igraph::vertex_attr(glycan, "sub")
  if (any(is.na(subs))) {
    rlang::abort("Glycan structure must have no NA in vertex attribute 'sub'.")
  }
  
  # Check if all substituents are valid
  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue("Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}")
    rlang::abort(msg, subs = invalid_subs)
  }
  
  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attrs(glycan, "linkage")) {
    rlang::abort("Glycan structure must have an edge attribute 'linkage'.")
  }
  
  # Check if no NA in edge attribute "linkage"
  linkages <- igraph::edge_attr(glycan, "linkage")
  if (any(is.na(linkages))) {
    rlang::abort("Glycan structure must have no NA in edge attribute 'linkage'.")
  }
  
  # Check if all linkages are valid
  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
  
  # Check if "anomer" attribute exists
  if (is.null(glycan$anomer)) {
    rlang::abort("Glycan structure must have a graph attribute 'anomer'.")
  }
  
  # Check if "anomer" attribute is valid
  if (!valid_anomer(glycan$anomer)) {
    rlang::abort(glue::glue("Invalid anomer: {glycan$anomer}"))
  }
  
  # Check if "alditol" attribute exists
  if (is.null(glycan$alditol)) {
    rlang::abort("Glycan structure must have a graph attribute 'alditol'.")
  }
  
  # Check if "alditol" attribute is valid
  if (!is.logical(glycan$alditol)) {
    rlang::abort("Glycan structure attribute 'alditol' must be logical.")
  }

  glycan
}

# Helper function to create a new glycan structure vector
new_glycan_structure <- function(iupac = character(), mono_type = character(), structures = list()) {
  vctrs::new_rcrd(
    list(iupac = iupac, mono_type = mono_type),
    structures = structures,
    class = "glyrepr_structure"
  )
}

# Helper function to extract structures from existing vector
get_structures_from_vector <- function(x) {
  if (!is_glycan_structure(x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(x, "structures")
  
  # Return list of individual structures corresponding to each element
  purrr::map(codes, ~ structures[[.x]])
}

ensure_name_vertex_attr <- function(glycan) {
  if (!("name" %in% igraph::vertex_attr_names(glycan))) {
    names <- as.character(seq_len(igraph::vcount(glycan)))
    glycan <- igraph::set_vertex_attr(glycan, "name", value = names)
  }
  glycan
}

#' @export 
#' @rdname glycan_structure
is_glycan_structure <- function(x) {
  inherits(x, "glyrepr_structure")
}

#' Convert to Glycan Structure Vector
#'
#' Convert an object to a glycan structure vector.
#'
#' @param x An object to convert to a glycan structure vector.
#'   Can be an igraph object, a list of igraph objects,
#'   a character vector of IUPAC-condensed strings,
#'   or an existing glyrepr_structure object.
#'
#' @return A glyrepr_structure object.
#'
#' @examples
#' library(igraph)
#' 
#' # Convert a single igraph
#' graph <- make_graph(~ 1-+2)
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")
#' V(graph)$sub <- ""
#' E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' as_glycan_structure(graph)
#' 
#' # Convert a list of igraphs
#' o_glycan_vec <- o_glycan_core_1()
#' o_glycan_graph <- get_structure_graphs(o_glycan_vec, 1)
#' as_glycan_structure(list(graph, o_glycan_graph))
#' 
#' # Convert a character vector of IUPAC-condensed strings
#' as_glycan_structure(c("GlcNAc(b1-4)GlcNAc", "Man(a1-2)GlcNAc"))
#'
#' @export
as_glycan_structure <- function(x) {
  UseMethod("as_glycan_structure")
}

#' @export
as_glycan_structure.glyrepr_structure <- function(x) {
  x
}

#' @export
as_glycan_structure.igraph <- function(x) {
  glycan_structure(x)
}

#' @export
as_glycan_structure.list <- function(x) {
  # Validate that all elements are igraph objects
  if (!all(purrr::map_lgl(x, ~ inherits(.x, "igraph")))) {
    rlang::abort(c(
      "All elements in the list must be igraph objects.",
      "i" = "Each graph in the list should be a valid glycan structure."
    ))
  }
  do.call(glycan_structure, x)
}

#' @export
as_glycan_structure.character <- function(x) {
  if (length(x) == 1) {
    graph <- .parse_iupac_condensed_single(x)
    return(glycan_structure(graph))
  } else {
    # Multiple characters - return list of structures
    graphs <- purrr::map(x, .parse_iupac_condensed_single)
    return(do.call(glycan_structure, graphs))
  }
}

#' @export
as_glycan_structure.default <- function(x) {
  rlang::abort(c(
    "Cannot convert object of class {.cls {class(x)}} to glyrepr_structure.",
    "i" = "Supported types: igraph object, list of igraph objects, character vector (IUPAC-condensed), or existing glyrepr_structure."
  ))
}

#' @export
vec_ptype_full.glyrepr_structure <- function(x, ...) "glycan_structure"

#' @export
vec_ptype_abbr.glyrepr_structure <- function(x, ...) "structure"

#' @export
format.glyrepr_structure <- function(x, ...) {
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "iupac")
  unname(codes)
}

#' Format a Subset of Glycan Structures with Optional Colors
#'
#' @param x A glyrepr_structure object
#' @param indices Indices of structures to format
#' @param colored A logical value indicating whether to add colors
#' @return A character vector of formatted structures for the specified indices
#' @keywords internal
format_glycan_structure_subset <- function(x, indices, colored = TRUE) {
  if (!colored) {
    return(format(x)[indices])
  }
  
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "iupac")[indices]
  mono_types <- vctrs::field(data, "mono_type")[indices]
  structures <- attr(x, "structures")
  
  # For each structure, add colors if concrete type
  purrr::map2_chr(codes, mono_types, function(code, mono_type) {
    structure <- structures[[code]]
    mono_names <- igraph::V(structure)$mono
    
    # Add colors to monosaccharides and gray linkages
    if (colored) {
      colorize_iupac_string(code, mono_names)
    } else {
      code
    }
  })
}


#' @export
obj_print_footer.glyrepr_structure <- function(x, ...) {
  cat("# Unique structures: ", format(length(attr(x, "structures"))), "\n", sep = "")
}

#' @export
obj_print_data.glyrepr_structure <- function(x, ..., max_n = 10, colored = TRUE) {
  if (length(x) == 0) {
    return()
  }

  n <- length(x)
  n_show <- min(n, max_n)
  
  # Only format the structures that need to be shown to improve performance
  indices_to_show <- seq_len(n_show)
  formatted <- format_glycan_structure_subset(x, indices_to_show, colored = colored)
  
  # Print each IUPAC structure on its own line with indexing, up to max_n
  for (i in seq_len(n_show)) {
    cat("[", i, "] ", formatted[i], "\n", sep = "")
  }
  if (n > max_n) {
    cat("... (", n - max_n, " more not shown)\n", sep = "")
  }
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.glyrepr_structure <- function(x, ...) {
  pillar::pillar_shaft(format(x))
}

#' @export
vec_ptype2.glyrepr_structure.glyrepr_structure <- function(x, y, ...) {
  x_structures <- attr(x, "structures")
  y_structures <- attr(y, "structures")
  
  # Combine all structures from both x and y
  all_structures <- c(x_structures, y_structures)
  
  # Remove duplicates (keep first occurrence)
  unique_names <- unique(names(all_structures))
  combined_structures <- all_structures[unique_names]
  names(combined_structures) <- unique_names
  
  # Create prototype with all structures
  new_glycan_structure(character(), character(), structures = combined_structures)
}

#' @export
vec_cast.glyrepr_structure.glyrepr_structure <- function(x, to, ...) {
  x
}

#' @export
vec_restore.glyrepr_structure <- function(x, to, ...) {
  # x is the proxy data (data.frame) after vctrs operations
  # to is the original glyrepr_structure object for reference

  # Extract data from the proxy
  data <- vctrs::vec_data(x)
  iupacs <- vctrs::field(data, "iupac")
  mono_types <- vctrs::field(data, "mono_type")

  # Get available structures
  available_structures <- attr(to, "structures")

  # Find which unique structures are still needed
  needed_codes <- unique(iupacs)

  # For slicing operations, optimize by removing unused structures
  # For combination operations, we typically get all the needed structures already
  if (length(needed_codes) > 0 && length(available_structures) > 0) {
    # Only keep structures that are actually used
    retrenched_structures <- available_structures[needed_codes[needed_codes %in% names(available_structures)]]
  } else {
    retrenched_structures <- available_structures
  }

  # Create new vector with structures
  new_glycan_structure(iupacs, mono_types, structures = retrenched_structures)
}


#' @export
as.character.glyrepr_structure <- function(x, ...) {
  data <- vctrs::vec_data(x)
  vctrs::field(data, "iupac")
}

#' Access Individual Glycan Structures
#'
#' Extract individual glycan structure graphs from a glycan structure vector.
#'
#' @param x A glycan structure vector.
#' @param i Index or indices of structures to extract.
#'
#' @return A list of glycan_structure (igraph) objects.
#'
#' @examples
#' structures <- glycan_structure(o_glycan_core_1(), n_glycan_core())
#' get_structure_graphs(structures, 1)
#' get_structure_graphs(structures, c(1, 2))
#'
#' @export
get_structure_graphs <- function(x, i = NULL) {
  if (!is_glycan_structure(x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(x, "structures")
  
  if (is.null(i)) {
    i <- seq_along(codes)
  }
  
  selected_codes <- codes[i]
  result <- purrr::map(selected_codes, ~ structures[[.x]])
  
  if (length(result) == 1) {
    return(result[[1]])
  } else {
    return(result)
  }
}
