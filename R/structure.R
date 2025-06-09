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
#' `c("glyrepr_structure", "vctrs_rcrd", "vctrs_vctr")`.
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
    return(new_glycan_structure(character(), character(), list()))
  }
  
  # Validate and process each graph
  processed_graphs <- purrr::map(graphs, function(graph) {
    checkmate::assert_class(graph, "igraph")
    graph %>%
      validate_single_glycan_structure() %>%
      ensure_name_vertex_attr()
  })
  
  # Use hash values of IUPAC codes for deduplication
  iupacs <- purrr::map_chr(processed_graphs, .structure_to_iupac_single)
  codes <- purrr::map_chr(iupacs, rlang::hash)
  
  # Create a unique list based on uniqueness of codes
  unique_indices <- which(!duplicated(codes))
  unique_graphs <- processed_graphs[unique_indices]
  names(unique_graphs) <- codes[unique_indices]
  unique_iupacs <- iupacs[unique_indices]
  names(unique_iupacs) <- codes[unique_indices]
  
  new_glycan_structure(codes, unique_iupacs, unique_graphs)
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
new_glycan_structure <- function(codes = character(), iupacs = character(), structures = list()) {
  vctrs::new_rcrd(
    list(codes = codes),
    iupacs = iupacs,
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
  codes <- vctrs::field(data, "codes")
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
as_glycan_structure.default <- function(x) {
  rlang::abort(c(
    "Cannot convert object of class {.cls {class(x)}} to glyrepr_structure.",
    "i" = "Supported types: igraph object, list of igraph objects, or existing glyrepr_structure."
  ))
}

#' @export
vec_ptype_full.glyrepr_structure <- function(x, ...) "glycan_structure"

#' @export
vec_ptype_abbr.glyrepr_structure <- function(x, ...) "structure"

#' @export
format.glyrepr_structure <- function(x, ...) {
  data <- vctrs::vec_data(x)
  codes <- vctrs::field(data, "codes")
  iupacs <- attr(x, "iupacs")
  unname(iupacs[codes])
}

#' @export
obj_print_footer.glyrepr_structure <- function(x, ...) {
  cat("# Unique structures: ", format(length(attr(x, "iupacs"))), "\n", sep = "")
}

#' @export
obj_print_data.glyrepr_structure <- function(x, ..., max_n = 10) {
  if (length(x) == 0) {
    return()
  }

  formatted <- format(x)
  n <- length(formatted)
  n_show <- min(n, max_n)
  # Print each IUPAC structure on its own line with indexing, up to max_n
  for (i in seq_len(n_show)) {
    cat("[", i, "] ", formatted[i], "\n", sep = "")
  }
  if (n > max_n) {
    cat("... (", n - max_n, " more not shown)\n", sep = "")
  }
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
  codes <- vctrs::field(data, "codes")
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

#' Validate Glycan Structure (Legacy Function)
#'
#' @description
#' This is a legacy function for backward compatibility. 
#' It validates individual glycan structure igraph objects.
#'
#' @param glycan An igraph object or glyrepr_structure vector
#'
#' @return The validated glycan structure(s)
#'
#' @keywords internal
#' @export
validate_glycan_structure <- function(glycan) {
  UseMethod("validate_glycan_structure")
}

#' @export
validate_glycan_structure.igraph <- function(glycan) {
  validate_single_glycan_structure(glycan)
}

#' @export
validate_glycan_structure.glyrepr_structure <- function(glycan) {
  # For vectorized structures, validate each individual structure
  data <- vctrs::vec_data(glycan)
  codes <- vctrs::field(data, "codes")
  structures <- attr(glycan, "structures")
  
  # Validate each unique structure
  for (code in names(structures)) {
    validate_single_glycan_structure(structures[[code]])
  }
  
  glycan
}

#' @export
validate_glycan_structure.default <- function(glycan) {
  # Try to handle legacy glycan_structure objects (single igraphs with glycan_structure class)
  if (inherits(glycan, "glycan_structure") && inherits(glycan, "igraph")) {
    return(validate_single_glycan_structure(glycan))
  }
  
  rlang::abort(c(
    "Cannot validate object of class {.cls {class(glycan)}}.",
    "i" = "Supported types: igraph object or glyrepr_structure vector."
  ))
}

#' Map Functions Over Glycan Structure Vectors
#'
#' @description
#' These functions apply a function to each unique structure in a glycan structure vector,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr mapping functions, but optimized for glycan structure vectors.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .f A function that takes an igraph object and returns a result.
#' @param ... Additional arguments passed to `.f`.
#' @param .ptype A prototype for the return type (for `structure_map_vec`).
#'
#' @details
#' These functions only compute `.f` once for each unique structure, then map
#' the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element individually when there are duplicate structures.
#'
#' - `structure_map()`: Returns a list with the same length as `.x`
#' - `structure_map_vec()`: Returns an atomic vector with the same length as `.x`
#' - `structure_map_lgl()`: Returns a logical vector
#' - `structure_map_int()`: Returns an integer vector  
#' - `structure_map_dbl()`: Returns a double vector
#' - `structure_map_chr()`: Returns a character vector
#' - `structure_map_structure()`: Returns a new glycan structure vector (`.f` must return igraph objects)
#'
#' @return 
#' - `structure_map()`: A list
#' - `structure_map_vec()`: An atomic vector of type specified by `.ptype`
#' - `structure_map_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
#' - `structure_map_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' 
#' # Map a function that counts vertices - only computed twice, not three times
#' structure_map_int(structures, igraph::vcount)
#' 
#' # Map a function that returns logical
#' structure_map_lgl(structures, function(g) igraph::vcount(g) > 5)
#' 
#' # Map a function that modifies structure (must return igraph)
#' add_vertex_names <- function(g) {
#'   if (!("name" %in% igraph::vertex_attr_names(g))) {
#'     igraph::set_vertex_attr(g, "name", value = paste0("v", seq_len(igraph::vcount(g))))
#'   } else {
#'     g
#'   }
#' }
#' structure_map_structure(structures, add_vertex_names)
#'
#' @name structure_map
NULL

#' @rdname structure_map
#' @export
structure_map <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  if (!is.function(.f)) {
    rlang::abort("`.f` must be a function.")
  }
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "codes")
  structures <- attr(.x, "structures")
  
  # Apply function only to unique structures
  unique_codes <- names(structures)
  unique_results <- purrr::map(unique_codes, function(code) {
    .f(structures[[code]], ...)
  })
  names(unique_results) <- unique_codes
  
  # Map results back to original vector positions
  purrr::map(codes, ~ unique_results[[.x]])
}

#' @rdname structure_map
#' @export
structure_map_vec <- function(.x, .f, ..., .ptype = NULL) {
  results <- structure_map(.x, .f, ...)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname structure_map
#' @export
structure_map_lgl <- function(.x, .f, ...) {
  structure_map_vec(.x, .f, ..., .ptype = logical())
}

#' @rdname structure_map
#' @export
structure_map_int <- function(.x, .f, ...) {
  structure_map_vec(.x, .f, ..., .ptype = integer())
}

#' @rdname structure_map
#' @export
structure_map_dbl <- function(.x, .f, ...) {
  structure_map_vec(.x, .f, ..., .ptype = double())
}

#' @rdname structure_map
#' @export
structure_map_chr <- function(.x, .f, ...) {
  structure_map_vec(.x, .f, ..., .ptype = character())
}

#' @rdname structure_map
#' @export
structure_map_structure <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  if (!is.function(.f)) {
    rlang::abort("`.f` must be a function.")
  }
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "codes")
  structures <- attr(.x, "structures")
  
  # Apply function only to unique structures
  unique_codes <- names(structures)
  transformed_structures <- purrr::map(unique_codes, function(code) {
    result <- .f(structures[[code]], ...)
    if (!inherits(result, "igraph")) {
      rlang::abort("Function `.f` must return an igraph object when using `structure_map_structure()`.")
    }
    result
  })
  
  # Create new glycan structure vector from transformed structures
  # Map back to original positions
  individual_structures <- purrr::map(codes, ~ transformed_structures[[which(unique_codes == .x)]])
  
  do.call(glycan_structure, individual_structures)
}

#' Apply Function to Unique Structures Only
#'
#' @description
#' Apply a function only to the unique structures in a glycan structure vector,
#' returning results in the same order as the unique structures appear.
#' This is useful when you need to perform expensive computations but only
#' care about unique results.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .f A function that takes an igraph object and returns a result.
#' @param ... Additional arguments passed to `.f`.
#'
#' @return A list with results for each unique structure, named by their hash codes.
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' structures <- glycan_structure(core1, core1, core1)  # same structure 3 times
#' 
#' # Only compute once for the unique structure
#' unique_results <- structure_map_unique(structures, igraph::vcount)
#' length(unique_results)  # 1, not 3
#'
#' @export
structure_map_unique <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  if (!is.function(.f)) {
    rlang::abort("`.f` must be a function.")
  }
  
  structures <- attr(.x, "structures")
  
  # Apply function only to unique structures
  results <- purrr::map(structures, .f, ...)
  results
}
