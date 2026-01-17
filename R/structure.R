#' Create a Glycan Structure Vector
#'
#' @description
#' `glycan_structure()` creates an efficient glycan structure vector for storing and
#' processing glycan molecular structures. The function employs hash-based deduplication
#' mechanisms, making it suitable for glycoproteomics, glycomics analysis, and glycan
#' structure comparison studies.
#'
#' @details
#' # Core Features
#'
#' - **Efficient Storage**: Uses hash values of IUPAC codes for deduplication, 
#'   avoiding redundant storage of identical glycan structures
#' - **Graph Model Representation**: Each glycan structure is represented as a directed 
#'   graph where nodes are monosaccharides and edges are glycosidic linkages
#' - **Vectorized Operations**: Supports R's vectorized operations for batch 
#'   processing of glycan data
#' - **Type Safety**: Built on the vctrs package, providing type-safe operations
#'
#' # Data Structure Overview
#'
#' A glycan structure vector is a vctrs record with an additional S3 class 
#' `glyrepr_structure`. Therefore, `sloop::s3_class()` returns the class hierarchy 
#' `c("glyrepr_structure", "vctrs_rcrd")`.
#'
#' Each glycan structure must satisfy the following constraints:
#'
#' ## Graph Structure Requirements
#' - Must be a directed graph with an outward tree structure (reducing end as root)
#' - Must have a graph attribute `anomer` in the format "a1" or "b1"
#'   - Unknown parts can be represented with "?", e.g., "?1", "a?", "??"
#'
#' ## Node Attributes
#' - `mono`: Monosaccharide names, must be known monosaccharide types
#'   - Generic names: Hex, HexNAc, dHex, NeuAc, etc.
#'   - Concrete names: Glc, Gal, Man, GlcNAc, etc.
#'   - Cannot mix generic and concrete names
#'   - NA values are not allowed
#' - `sub`: Substituent information
#'   - Single substituent format: "xY" (x = position, Y = substituent name), 
#'     e.g., "2Ac", "3S"
#'   - Multiple substituents separated by commas and ordered by position, 
#'     e.g., "3Me,4Ac", "2S,6P"
#'   - No substituents represented by empty string ""
#'
#' ## Edge Attributes
#' - `linkage`: Glycosidic linkage information in format "a/bX-Y"
#'   - Standard format: e.g., "b1-4", "a2-3"
#'   - Unknown positions allowed: "a1-?", "b?-3", "??-?"
#'   - Partially unknown positions: "a1-3/6", "a1-3/6/9"
#'   - NA values are not allowed
#'
#' # Node and Edge Order
#'
#' The indices of vertices and linkages in a glycan correspond directly to their
#' order in the IUPAC-condensed string, which is printed when you print a
#' [glyrepr::glycan_structure()].
#' For example, for the glycan `Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-`,
#' the vertices are "Man", "Man", "Man", "GlcNAc", "GlcNAc",
#' and the linkages are "a1-3", "a1-6", "b1-4", "b1-4".
#'
#' # Use Cases
#'
#' - **Glycoproteomics Analysis**: Processing glycan structure information from 
#'   mass spectrometry data
#' - **Glycomics Research**: Comparing glycan expression profiles across different 
#'   samples or conditions
#' - **Structure-Function Analysis**: Studying relationships between glycan 
#'   structures and biological functions
#' - **Database Queries**: Performing structure matching and searches in glycan 
#'   databases
#'
#' @param ... igraph graph objects to be converted to glycan structures, or existing 
#'   glycan structure vectors. Supports mixed input of multiple objects.
#' @param x An object to check or convert.
#'
#' @returns A `glyrepr_structure` class glycan structure vector object.
#'
#' @examples
#' library(igraph)
#'
#' # Example 1: Create a simple glycan structure GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2)  # Create graph with two monosaccharides
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")  # Set monosaccharide types
#' V(graph)$sub <- ""  # No substituents
#' E(graph)$linkage <- "b1-4"  # b1-4 glycosidic linkage
#' graph$anomer <- "a1"  # a anomeric carbon
#'
#' # Create glycan structure vector
#' simple_struct <- glycan_structure(graph)
#' print(simple_struct)
#'
#' # Example 2: Use predefined glycan core structures
#' n_core <- n_glycan_core()  # N-glycan core structure
#' o_core1 <- o_glycan_core_1()  # O-glycan Core 1 structure
#'
#' # Create vector with multiple structures
#' multi_struct <- glycan_structure(n_core, o_core1)
#' print(multi_struct)
#'
#' # Example 3: Create complex structure with substituents
#' complex_graph <- make_graph(~ 1-+2-+3)
#' V(complex_graph)$mono <- c("GlcNAc", "Gal", "Neu5Ac")
#' V(complex_graph)$sub <- c("", "", "")  # Add substituents as needed
#' E(complex_graph)$linkage <- c("b1-4", "a2-3")
#' complex_graph$anomer <- "b1"
#'
#' complex_struct <- glycan_structure(complex_graph)
#' print(complex_struct)
#'
#' # Example 4: Check if object is a glycan structure
#' is_glycan_structure(simple_struct)  # TRUE
#' is_glycan_structure(graph)          # FALSE
#'
#' # Example 5: Mix different input types
#' mixed_struct <- glycan_structure(graph, o_glycan_core_2(), simple_struct)
#' print(mixed_struct)
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
      cli::cli_abort("All arguments must be igraph objects or glycan_structure vectors.")
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

  res <- new_glycan_structure(iupacs, mono_types, unique_graphs)
  reorder_graphs(res)
}

# Helper function to validate a single glycan structure
validate_single_glycan_structure <- function(glycan) {
  checkmate::assert_class(glycan, "igraph")

  # Check if it is a directed graph
  if (!is_directed_graph(glycan)) {
    cli::cli_abort("Glycan structure must be directed.")
  }

  # Check if it is an out tree
  if (!is_out_tree(glycan)) {
    cli::cli_abort("Glycan structure must be an out tree.")
  }

  # Check if graph has a vertex attribute "mono"
  # This is the monosaccharide name, e.g. "GlcNAc", "Man", etc.
  if (!has_vertex_attrs(glycan, "mono")) {
    cli::cli_abort("Glycan structure must have a vertex attribute 'mono'.")
  }

  # Check if no NA in vertex attribute "mono"
  mono_names <- igraph::vertex_attr(glycan, "mono")
  if (any(is.na(mono_names))) {
    cli::cli_abort("Glycan structure must have no NA in vertex attribute 'mono'.")
  }

  # Check if all monosaccharides are known
  if (!all(is_known_mono(mono_names))) {
    unknown_monos <- unique(igraph::V(glycan)$mono[!is_known_mono(igraph::V(glycan)$mono)])
    msg <- glue::glue("Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}")
    cli::cli_abort(msg, monos = unknown_monos)
  }

  # Check if mixed use of generic and concrete monosaccharides
  if (mix_generic_concrete(mono_names)) {
    cli::cli_abort("Monosaccharides must be either all generic or all concrete.")
  }

  # Check if graph has a vertex attribute "sub"
  # This is the substituent name, e.g. "Ac", "S", "P", or "" (no).
  if (!has_vertex_attrs(glycan, "sub")) {
    cli::cli_abort("Glycan structure must have a vertex attribute 'sub'.")
  }

  # Check if no NA in vertex attribute "sub"
  subs <- igraph::vertex_attr(glycan, "sub")
  if (any(is.na(subs))) {
    cli::cli_abort("Glycan structure must have no NA in vertex attribute 'sub'.")
  }

  # Check if all substituents are valid
  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue("Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}")
    cli::cli_abort(msg, subs = invalid_subs)
  }

  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attrs(glycan, "linkage")) {
    cli::cli_abort("Glycan structure must have an edge attribute 'linkage'.")
  }

  # Check if no NA in edge attribute "linkage"
  linkages <- igraph::edge_attr(glycan, "linkage")
  if (any(is.na(linkages))) {
    cli::cli_abort("Glycan structure must have no NA in edge attribute 'linkage'.")
  }

  # Check if all linkages are valid
  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    cli::cli_abort(msg, linkages = invalid_linkages)
  }

  # Check if any duplicated linkage positions exist
  if (any_dup_linkage_pos(glycan)) {
    cli::cli_abort("Duplicated linkage positions.")
  }

  # Check if "anomer" attribute exists
  if (is.null(glycan$anomer)) {
    cli::cli_abort("Glycan structure must have a graph attribute 'anomer'.")
  }

  # Check if "anomer" attribute is valid
  if (!valid_anomer(glycan$anomer)) {
    cli::cli_abort(glue::glue("Invalid anomer: {glycan$anomer}"))
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
    cli::cli_abort("Input must be a glycan_structure vector.")
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
#' @returns A glyrepr_structure object.
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
#' as_glycan_structure(graph)
#'
#' # Convert a list of igraphs
#' o_glycan_vec <- o_glycan_core_1()
#' o_glycan_graph <- get_structure_graphs(o_glycan_vec)
#' as_glycan_structure(list(graph, o_glycan_graph))
#'
#' # Convert a character vector of IUPAC-condensed strings
#' as_glycan_structure(c("GlcNAc(b1-4)GlcNAc(b1-", "Man(a1-2)GlcNAc(b1-"))
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
    cli::cli_abort(c(
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
  cli::cli_abort(c(
    "Cannot convert object of class {.cls {class(x)}} to glyrepr_structure.",
    "i" = "Supported types: igraph object, list of igraph objects, character vector (IUPAC-condensed), or existing glyrepr_structure."
  ))
}

#' @export
vec_ptype_full.glyrepr_structure <- function(x, ...) "glycan_structure"

#' @export
vec_ptype_abbr.glyrepr_structure <- function(x, ...) "struct"

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
#' @returns A character vector of formatted structures for the specified indices
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
  if (length(x) == 0) {
    return(pillar::pillar_shaft(character()))
  }

  # Get formatted strings with colors
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "iupac")
  mono_types <- vctrs::field(data, "mono_type")
  structures <- attr(x, "structures")

  # For each structure, add colors if concrete type
  formatted <- purrr::map2_chr(codes, mono_types, function(code, mono_type) {
    structure <- structures[[code]]
    mono_names <- igraph::V(structure)$mono

    # Add colors to monosaccharides and gray linkages
    colorize_iupac_string(code, mono_names)
  })

  pillar::new_pillar_shaft_simple(formatted, align = "left", min_width = 10)
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
#' @param return_list If `TRUE`, always returns a list.
#'   If `FALSE` and `x` has a length of 1, return the igraph object directly.
#'   If not provided (default), `FALSE` when `x` has a length of 1 and `TRUE` otherwise.
#'
#' @returns A list of igraph objects or an igraph object directly (see `return_list` parameter).
#'
#' @examples
#' structures <- glycan_structure(o_glycan_core_1(), n_glycan_core())
#' get_structure_graphs(structures)
#' get_structure_graphs(structures)
#'
#' @export
get_structure_graphs <- function(x, return_list = NULL) {
  checkmate::assert_class(x, "glyrepr_structure")
  checkmate::assert_flag(return_list, null.ok = TRUE)

  if (is.null(return_list)) {
    return_list <- length(x) > 1
  } else {
    if (!return_list && length(x) > 1) {
      cli::cli_abort(c(
        "{.arg return_list} must be `TRUE` or `NULL` if {.arg x} has a length greater than 1.",
        "i" = "Length of {.arg x}: {.val {length(x)}}"
      ))
    }
  }

  data <- vctrs::vec_data(x)
  iupacs <- vctrs::field(data, "iupac")
  structures <- attr(x, "structures")
  res <- purrr::map(iupacs, ~ structures[[.x]])
  if (!return_list) {
    res <- res[[1]]
  }
  res
}
