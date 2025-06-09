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
  
  codes <- vctrs::vec_data(.x)
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
  
  codes <- vctrs::vec_data(.x)
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
