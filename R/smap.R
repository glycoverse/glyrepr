#' Map Functions Over Glycan Structure Vectors
#'
#' @description
#' These functions apply a function to each unique structure in a glycan structure vector,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr mapping functions, but optimized for glycan structure vectors.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .f A function that takes an igraph object and returns a result. 
#'   Can be a function, purrr-style lambda (`~ .x$attr`), or a character string naming a function.
#' @param ... Additional arguments passed to `.f`.
#' @param .ptype A prototype for the return type (for `smap_vec`).
#' @param .parallel Logical; whether to use parallel processing. If `FALSE` (default), 
#'   parallel processing is disabled. Set to `TRUE` to enable parallel processing.
#'
#' @details
#' These functions only compute `.f` once for each unique structure, then map
#' the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element individually when there are duplicate structures.
#'
#'
#' **Return Types:**
#' - `smap()`: Returns a list with the same length as `.x`
#' - `smap_vec()`: Returns an atomic vector with the same length as `.x`
#' - `smap_lgl()`: Returns a logical vector
#' - `smap_int()`: Returns an integer vector  
#' - `smap_dbl()`: Returns a double vector
#' - `smap_chr()`: Returns a character vector
#' - `smap_structure()`: Returns a new glycan structure vector (`.f` must return igraph objects)
#'
#' @return 
#' - `smap()`: A list
#' - `smap_vec()`: An atomic vector of type specified by `.ptype`
#' - `smap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
#' - `smap_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' 
#' # Map a function that counts vertices - only computed twice, not three times
#' smap_int(structures, igraph::vcount)
#' 
#' # Map a function that returns logical
#' smap_lgl(structures, function(g) igraph::vcount(g) > 5)
#' 
#' # Use purrr-style lambda functions  
#' smap_int(structures, ~ igraph::vcount(.x))
#' smap_lgl(structures, ~ igraph::vcount(.x) > 5)
#' 
#' # Map a function that modifies structure (must return igraph)
#' add_vertex_names <- function(g) {
#'   if (!("name" %in% igraph::vertex_attr_names(g))) {
#'     igraph::set_vertex_attr(g, "name", value = paste0("v", seq_len(igraph::vcount(g))))
#'   } else {
#'     g
#'   }
#' }
#' smap_structure(structures, add_vertex_names)
#' 
#' \dontrun{
#' # Parallel processing examples
#' 
#' # Example 1: Manual control of parallel processing
#' # For small datasets, you can still force parallel processing
#' 
#' # Force sequential processing (default for small datasets)
#' result_seq <- smap_int(structures, igraph::vcount, .parallel = FALSE)
#' 
#' # Force parallel processing (requires parallel backend setup)
#' # First set up a parallel backend:
#' # future::plan(future::multisession, workers = 4)
#' result_par <- smap_int(structures, igraph::vcount, .parallel = TRUE)
#' 
#' # Example 2: Manual parallel processing for large datasets
#' # Parallel processing is disabled by default - you need to explicitly enable it
#' 
#' # Simulate a large dataset (in practice, this might come from parsing 
#' # a large glycomics dataset with many different structures)
#' # Note: This is just for demonstration - in real use, you'd have 
#' # genuinely different structures from your data
#' 
#' # For large datasets, you can enable parallel processing manually:
#' # large_structures <- your_large_glycan_dataset  # Many unique structures
#' # result_parallel <- smap_int(large_structures, igraph::vcount, .parallel = TRUE)
#' 
#' # Example 3: Complex computation that benefits from parallelization
#' complex_analysis <- function(g) {
#'   # Simulate expensive computation (e.g., complex graph analysis)
#'   Sys.sleep(0.01)  # 10ms delay per unique structure
#'   list(
#'     vcount = igraph::vcount(g),
#'     ecount = igraph::ecount(g),
#'     diameter = igraph::diameter(g),
#'     transitivity = igraph::transitivity(g)
#'   )
#' }
#' 
#' # Even with few unique structures, complex computations benefit from parallel
#' # Set up parallel backend first: future::plan(future::multisession, workers = 4)
#' system.time({
#'   results_par <- smap(structures, complex_analysis, .parallel = TRUE)
#' })
#' 
#' # Compare with sequential processing
#' system.time({
#'   results_seq <- smap(structures, complex_analysis, .parallel = FALSE)
#' })
#' 
#' # Example 4: Practical use case with real data
#' # When processing large glycomics datasets:
#' # glycan_data <- parse_glycan_file("large_dataset.txt")  # Many structures
#' # 
#' # # Compute molecular properties - manually enable parallel if needed
#' # molecular_weights <- smap_dbl(glycan_data, calculate_molecular_weight, .parallel = TRUE)
#' # 
#' # # Force parallel for expensive computations even with fewer structures  
#' # complexity_scores <- smap_dbl(glycan_data, expensive_complexity_analysis, 
#' #                               .parallel = TRUE)
#' }
#'
#' @name smap
NULL

# Helper function to rebuild glycan_structure with proper deduplication
# after modifications that may create identical structures
.rebuild_structure_with_dedup <- function(modified_structures, idx_mapping) {
  # Get new IUPACs for all modified structures
  new_unique_iupacs <- purrr::map_chr(modified_structures, .structure_to_iupac_single)
  new_iupacs <- new_unique_iupacs[idx_mapping]

  # Get new mono types for all modified structures
  new_unique_mono_types <- unname(purrr::map_chr(modified_structures, get_graph_mono_type))
  new_mono_types <- new_unique_mono_types[idx_mapping]

  # Re-deduplicate structures based on new IUPACs to handle cases where
  # modifications create identical structures
  unique_new_indices <- which(!duplicated(new_unique_iupacs))
  final_unique_structures <- modified_structures[unique_new_indices]
  final_unique_iupacs <- new_unique_iupacs[unique_new_indices]
  names(final_unique_structures) <- final_unique_iupacs

  # Create result using the standard glycan_structure constructor to ensure
  # proper deduplication and validation
  new_glycan_structure(new_iupacs, new_mono_types, final_unique_structures)
}

# Helper function for parallel processing
.smap_apply <- function(data_list, func, use_parallel = FALSE, auto_threshold = 100, ...) {
  n_tasks <- length(data_list)
  
  # Handle NULL values for backward compatibility
  if (is.null(use_parallel)) {
    use_parallel <- FALSE
  }
  
  if (use_parallel && n_tasks > 1) {
    # Check if future backend is set up
    if (!future::nbrOfWorkers() > 1) {
      cli::cli_inform("No parallel backend detected. Using sequential processing.")
      use_parallel <- FALSE
    }
  } else if (use_parallel && n_tasks <= 1) {
    use_parallel <- FALSE
  }
  
  if (use_parallel) {
    furrr::future_map(data_list, func, ...)
  } else {
    purrr::map(data_list, func, ...)
  }
}

# Helper function for common smap logic
.smap_base <- function(.x, .f, ..., .parallel = FALSE, .convert_fn = NULL) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Extract dots without .parallel
  dots <- list(...)
  
  # Apply function only to unique structures
  unique_codes <- names(structures)
  unique_results <- .smap_apply(unique_codes, function(code) {
    do.call(.f, c(list(structures[[code]]), dots))
  }, use_parallel = .parallel)
  names(unique_results) <- unique_codes
  
  # If no conversion function provided, return list (for smap)
  if (is.null(.convert_fn)) {
    # Optimized mapping: use match() instead of individual lookups
    idx <- match(codes, unique_codes)
    return(unname(unique_results[idx]))
  }
  
  # Convert to target type and map back to original positions
  unique_converted <- .convert_fn(unique_results)
  names(unique_converted) <- unique_codes
  
  # Optimized mapping: use match() instead of individual lookups
  idx <- match(codes, unique_codes)
  unname(unique_converted[idx])
}

#' @rdname smap
#' @export
smap <- function(.x, .f, ..., .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel, .convert_fn = NULL)
}

#' @rdname smap
#' @export
smap_vec <- function(.x, .f, ..., .ptype = NULL, .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel,
             .convert_fn = function(results) vctrs::vec_c(!!!results, .ptype = .ptype))
}

#' @rdname smap
#' @export
smap_lgl <- function(.x, .f, ..., .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel, 
             .convert_fn = function(results) as.logical(unlist(results, use.names = FALSE)))
}

#' @rdname smap
#' @export
smap_int <- function(.x, .f, ..., .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel,
             .convert_fn = function(results) as.integer(unlist(results, use.names = FALSE)))
}

#' @rdname smap
#' @export
smap_dbl <- function(.x, .f, ..., .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel,
             .convert_fn = function(results) as.double(unlist(results, use.names = FALSE)))
}

#' @rdname smap
#' @export
smap_chr <- function(.x, .f, ..., .parallel = FALSE) {
  .smap_base(.x, .f, ..., .parallel = .parallel,
             .convert_fn = function(results) as.character(unlist(results, use.names = FALSE)))
}

#' @rdname smap
#' @export
smap_structure <- function(.x, .f, ..., .parallel = FALSE) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .f <- rlang::as_function(.f)

  data <- vctrs::vec_data(.x)
  iupacs <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Extract dots without .parallel
  dots <- list(...)
  
  # Apply function only to unique structures
  unique_iupacs <- names(structures)
  new_structures <- .smap_apply(unique_iupacs, function(iupac) {
    result <- do.call(.f, c(list(structures[[iupac]]), dots))
    if (!inherits(result, "igraph")) {
      cli::cli_abort("Function `.f` must return an igraph object when using `smap_structure()`.")
    }
    result
  }, use_parallel = .parallel)
  
  # Rebuild glycan_structure with proper deduplication
  idx <- match(iupacs, unique_iupacs)
  .rebuild_structure_with_dedup(new_structures, idx)
}

#' Apply Function to Unique Structures Only
#'
#' @description
#' Apply a function only to the unique structures in a glycan structure vector,
#' returning results in the same order as the unique structures appear.
#' This is useful when you need to perform expensive computations but only
#' care about unique results.
#'
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .f A function that takes an igraph object and returns a result. 
#'   Can be a function, purrr-style lambda (`~ .x$attr`), or a character string naming a function.
#' @param ... Additional arguments passed to `.f`.
#' @param .parallel Logical; whether to use parallel processing. If `FALSE` (default), 
#'   parallel processing is disabled. Set to `TRUE` to enable parallel processing.
#'   See examples in \code{\link{smap}} for how to set up and use parallel processing.
#'
#' @return A list with results for each unique structure, named by their hash codes.
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' structures <- glycan_structure(core1, core1, core1)  # same structure 3 times
#' 
#' # Only compute once for the unique structure
#' unique_results <- smap_unique(structures, igraph::vcount)
#' length(unique_results)  # 1, not 3
#' 
#' # Use purrr-style lambda
#' unique_results2 <- smap_unique(structures, ~ igraph::vcount(.x))
#' length(unique_results2)  # 1, not 3
#'
#' @export
smap_unique <- function(.x, .f, ..., .parallel = FALSE) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }
  
  .f <- rlang::as_function(.f)
  
  structures <- attr(.x, "structures")
  
  # Extract dots without .parallel
  dots <- list(...)
  
  # Apply function only to unique structures
  results <- .smap_apply(structures, function(g) {
    do.call(.f, c(list(g), dots))
  }, use_parallel = .parallel)
  results
}

#' Test Predicates on Glycan Structure Vectors
#'
#' @description
#' These functions test predicates on unique structures in a glycan structure vector,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr predicate functions, but optimized for glycan structure vectors.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .p A predicate function that takes an igraph object and returns a logical value. 
#'   Can be a function, purrr-style lambda (`~ .x$attr`), or a character string naming a function.
#' @param ... Additional arguments passed to `.p`.
#'
#' @details
#' These functions only evaluate `.p` once for each unique structure, making them
#' much more efficient than applying `.p` to each element individually when there
#' are duplicate structures.
#' 
#' 
#' **Return Values:**
#' - `ssome()`: Returns `TRUE` if at least one unique structure satisfies the predicate
#' - `severy()`: Returns `TRUE` if all unique structures satisfy the predicate  
#' - `snone()`: Returns `TRUE` if no unique structures satisfy the predicate
#'
#' @return A single logical value.
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' 
#' # Test if some structures have more than 5 vertices
#' ssome(structures, function(g) igraph::vcount(g) > 5)
#' 
#' # Test if all structures have at least 3 vertices
#' severy(structures, function(g) igraph::vcount(g) >= 3)
#' 
#' # Test if no structures have more than 20 vertices
#' snone(structures, function(g) igraph::vcount(g) > 20)
#' 
#' # Use purrr-style lambda functions
#' ssome(structures, ~ igraph::vcount(.x) > 5)
#' severy(structures, ~ igraph::vcount(.x) >= 3)
#' snone(structures, ~ igraph::vcount(.x) > 20)
#'
#' @name smap_predicates
NULL

#' @rdname smap_predicates
#' @export
ssome <- function(.x, .p, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }
  
  .p <- rlang::as_function(.p)
  
  structures <- attr(.x, "structures")
  
  # Apply predicate only to unique structures using purrr::some
  purrr::some(structures, .p, ...)
}

#' @rdname smap_predicates
#' @export
severy <- function(.x, .p, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }
  
  .p <- rlang::as_function(.p)
  
  structures <- attr(.x, "structures")
  
  # Apply predicate only to unique structures using purrr::every
  purrr::every(structures, .p, ...)
}

#' @rdname smap_predicates
#' @export
snone <- function(.x, .p, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }
  
  .p <- rlang::as_function(.p)
  
  structures <- attr(.x, "structures")
  
  # Apply predicate only to unique structures using purrr::none
  purrr::none(structures, .p, ...)
}

#' Map Functions Over Two Glycan Structure Vectors
#'
#' @description
#' These functions apply a function to each unique structure combination in two glycan structure vectors,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr map2 functions, but optimized for glycan structure vectors.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .y A vector of the same length as `.x`, or length 1 (will be recycled).
#' @param .f A function that takes an igraph object (from `.x`) and a value (from `.y`) and returns a result. 
#'   Can be a function, purrr-style lambda (`~ .x + .y`), or a character string naming a function.
#' @param ... Additional arguments passed to `.f`.
#' @param .ptype A prototype for the return type (for `smap2_vec`).
#' @param .parallel Logical; whether to use parallel processing. If `FALSE` (default), 
#'   parallel processing is disabled. Set to `TRUE` to enable parallel processing.
#'   See examples in \code{\link{smap}} for how to set up and use parallel processing.
#'
#' @details
#' These functions only compute `.f` once for each unique combination of structure and corresponding
#' `.y` value, then map the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element pair individually when there are duplicate structure-value combinations.
#'
#'
#' **Return Types:**
#' - `smap2()`: Returns a list with the same length as `.x`
#' - `smap2_vec()`: Returns an atomic vector with the same length as `.x`
#' - `smap2_lgl()`: Returns a logical vector
#' - `smap2_int()`: Returns an integer vector  
#' - `smap2_dbl()`: Returns a double vector
#' - `smap2_chr()`: Returns a character vector
#' - `smap2_structure()`: Returns a new glycan structure vector (`.f` must return igraph objects)
#'
#' @return 
#' - `smap2()`: A list
#' - `smap2_vec()`: An atomic vector of type specified by `.ptype`
#' - `smap2_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
#' - `smap2_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create structure vectors with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' weights <- c(1.0, 2.0, 1.0)  # corresponding weights
#' 
#' # Map a function that uses both structure and weight
#' smap2_dbl(structures, weights, function(g, w) igraph::vcount(g) * w)
#' 
#' # Use purrr-style lambda functions  
#' smap2_dbl(structures, weights, ~ igraph::vcount(.x) * .y)
#' 
#' # Test with recycling (single weight for all structures)
#' smap2_dbl(structures, 2.5, ~ igraph::vcount(.x) * .y)
#' 
#' # Map a function that modifies structure based on second argument
#' # This example adds a graph attribute instead of modifying topology
#' add_weight_attr <- function(g, weight) {
#'   igraph::set_graph_attr(g, "weight", weight)
#' }
#' weights_to_add <- c(1.5, 2.5, 1.5)
#' smap2_structure(structures, weights_to_add, add_weight_attr)
#'
#' @name smap2
NULL

#' @rdname smap2
#' @export
smap2 <- function(.x, .y, .f, ..., .parallel = FALSE) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(list())
  }
  
  # Recycle .y to match length of .x
  .y <- vctrs::vec_recycle(.y, length(.x))
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    y_val = .y,
    stringsAsFactors = FALSE
  )
  combinations_df$combo_key <- paste0(combinations_df$code, "|||", combinations_df$y_val)
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  unique_results <- .smap_apply(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    .f(structures[[row$code]], row$y_val, ...)
  }, use_parallel = .parallel)
  names(unique_results) <- unique_combinations_df$combo_key
  
  # Map results back to original vector positions
  purrr::map(combinations_df$combo_key, ~ unique_results[[.x]])
}

#' @rdname smap2
#' @export
smap2_vec <- function(.x, .y, .f, ..., .ptype = NULL, .parallel = FALSE) {
  results <- smap2(.x, .y, .f, ..., .parallel = .parallel)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname smap2
#' @export
smap2_lgl <- function(.x, .y, .f, ..., .parallel = FALSE) {
  smap2_vec(.x, .y, .f, ..., .ptype = logical(), .parallel = .parallel)
}

#' @rdname smap2
#' @export
smap2_int <- function(.x, .y, .f, ..., .parallel = FALSE) {
  smap2_vec(.x, .y, .f, ..., .ptype = integer(), .parallel = .parallel)
}

#' @rdname smap2
#' @export
smap2_dbl <- function(.x, .y, .f, ..., .parallel = FALSE) {
  smap2_vec(.x, .y, .f, ..., .ptype = double(), .parallel = .parallel)
}

#' @rdname smap2
#' @export
smap2_chr <- function(.x, .y, .f, ..., .parallel = FALSE) {
  smap2_vec(.x, .y, .f, ..., .ptype = character(), .parallel = .parallel)
}

#' @rdname smap2
#' @export
smap2_structure <- function(.x, .y, .f, ..., .parallel = FALSE) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(glycan_structure())
  }
  
  # Recycle .y to match length of .x
  .y <- vctrs::vec_recycle(.y, length(.x))
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    y_val = .y,
    stringsAsFactors = FALSE
  )
  combinations_df$combo_key <- paste0(combinations_df$code, "|||", combinations_df$y_val)
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  new_structures <- .smap_apply(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    result <- .f(structures[[row$code]], row$y_val, ...)
    if (!inherits(result, "igraph")) {
      cli::cli_abort("Function `.f` must return an igraph object when using `smap2_structure()`.")
    }
    result
  }, use_parallel = .parallel)
  
  # Rebuild glycan_structure with proper deduplication
  idx <- match(combinations_df$combo_key, unique_combinations_df$combo_key)
  .rebuild_structure_with_dedup(new_structures, idx)
}

#' Map Functions Over Glycan Structure Vectors and Multiple Arguments
#'
#' @description
#' These functions apply a function to each unique structure in a glycan structure vector
#' along with corresponding elements from multiple other vectors,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr pmap functions, but optimized for glycan structure vectors.
#'
#' @param .l A list where the first element is a glycan structure vector (glyrepr_structure)
#'   and the remaining elements are vectors of the same length or length 1 (will be recycled).
#' @param .f A function that takes an igraph object (from first element of `.l`) and 
#'   values from other elements, returning a result. 
#'   Can be a function, purrr-style lambda (`~ .x + .y + .z`), or a character string naming a function.
#' @param ... Additional arguments passed to `.f`.
#' @param .ptype A prototype for the return type (for `spmap_vec`).
#' @param .parallel Logical; whether to use parallel processing. If `FALSE` (default), 
#'   parallel processing is disabled. Set to `TRUE` to enable parallel processing.
#'   See examples in \code{\link{smap}} for how to set up and use parallel processing.
#'
#' @details
#' These functions only compute `.f` once for each unique combination of structure and corresponding
#' values from other vectors, then map the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element combination individually when there are duplicate combinations.
#'
#' **Time Complexity Performance:**
#' 
#' Performance scales with unique combinations of all arguments rather than total vector length.
#' When argument vectors are highly redundant, performance approaches O(unique_structures).
#' Scaling factor shows time increase when vector size increases 20x.
#'
#' **Return Types:**
#' - `spmap()`: Returns a list with the same length as the input vectors
#' - `spmap_vec()`: Returns an atomic vector with the same length as the input vectors
#' - `spmap_lgl()`: Returns a logical vector
#' - `spmap_int()`: Returns an integer vector  
#' - `spmap_dbl()`: Returns a double vector
#' - `spmap_chr()`: Returns a character vector
#' - `spmap_structure()`: Returns a new glycan structure vector (`.f` must return igraph objects)
#'
#' @return 
#' - `spmap()`: A list
#' - `spmap_vec()`: An atomic vector of type specified by `.ptype`
#' - `spmap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
#' - `spmap_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create structure vectors with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' weights <- c(1.0, 2.0, 1.0)  # corresponding weights
#' factors <- c(2, 3, 2)  # corresponding factors
#' 
#' # Map a function that uses structure, weight, and factor
#' spmap_dbl(list(structures, weights, factors), 
#'           function(g, w, f) igraph::vcount(g) * w * f)
#' 
#' # Use purrr-style lambda functions  
#' spmap_dbl(list(structures, weights, factors), ~ igraph::vcount(..1) * ..2 * ..3)
#' 
#' # Test with recycling
#' spmap_dbl(list(structures, 2.0, 3), ~ igraph::vcount(..1) * ..2 * ..3)
#'
#' @name spmap
NULL

#' @rdname spmap
#' @export
spmap <- function(.l, .f, ..., .parallel = FALSE) {
  # Check if it's actually a plain list, not a vctrs object that looks like a list
  if (!inherits(.l, "list") || inherits(.l, "vctrs_vctr") || length(.l) == 0) {
    cli::cli_abort("Input `.l` must be a non-empty list.")
  }
  
  if (!is_glycan_structure(.l[[1]])) {
    cli::cli_abort("First element of `.l` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.l[[1]]) == 0) {
    return(list())
  }
  
  # Recycle all vectors to match length of first vector
  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.l[[1]])
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.l[[1]], "structures")
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    stringsAsFactors = FALSE
  )
  
  # Add other arguments as columns
  for (i in seq_along(.l)[-1]) {
    combinations_df[[paste0("arg", i)]] <- .l[[i]]
  }
  
  # Create combination key
  combinations_df$combo_key <- do.call(paste, c(combinations_df[setdiff(names(combinations_df), "combo_key")], sep = "|||"))
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  unique_results <- .smap_apply(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    # Build full argument list: first is structure, then other args
    args <- list(structures[[row$code]])
    if (length(.l) > 1) {
      for (j in 2:length(.l)) {
        args[[j]] <- row[[paste0("arg", j)]]
      }
    }
    do.call(.f, c(args, list(...)))
  }, use_parallel = .parallel)
  names(unique_results) <- unique_combinations_df$combo_key
  
  # Map results back to original vector positions
  purrr::map(combinations_df$combo_key, ~ unique_results[[.x]])
}

#' @rdname spmap
#' @export
spmap_vec <- function(.l, .f, ..., .ptype = NULL, .parallel = FALSE) {
  results <- spmap(.l, .f, ..., .parallel = .parallel)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname spmap
#' @export
spmap_lgl <- function(.l, .f, ..., .parallel = FALSE) {
  spmap_vec(.l, .f, ..., .ptype = logical(), .parallel = .parallel)
}

#' @rdname spmap
#' @export
spmap_int <- function(.l, .f, ..., .parallel = FALSE) {
  spmap_vec(.l, .f, ..., .ptype = integer(), .parallel = .parallel)
}

#' @rdname spmap
#' @export
spmap_dbl <- function(.l, .f, ..., .parallel = FALSE) {
  spmap_vec(.l, .f, ..., .ptype = double(), .parallel = .parallel)
}

#' @rdname spmap
#' @export
spmap_chr <- function(.l, .f, ..., .parallel = FALSE) {
  spmap_vec(.l, .f, ..., .ptype = character(), .parallel = .parallel)
}

#' @rdname spmap
#' @export
spmap_structure <- function(.l, .f, ..., .parallel = FALSE) {
  if (!is.list(.l) || length(.l) == 0) {
    cli::cli_abort("Input `.l` must be a non-empty list.")
  }
  
  if (!is_glycan_structure(.l[[1]])) {
    cli::cli_abort("First element of `.l` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.l[[1]]) == 0) {
    return(glycan_structure())
  }
  
  # Recycle all vectors to match length of first vector
  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.l[[1]])
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.l[[1]], "structures")
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    stringsAsFactors = FALSE
  )
  
  # Add other arguments as columns
  for (i in seq_along(.l)[-1]) {
    combinations_df[[paste0("arg", i)]] <- .l[[i]]
  }
  
  # Create combination key
  combinations_df$combo_key <- do.call(paste, c(combinations_df[setdiff(names(combinations_df), "combo_key")], sep = "|||"))
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  new_structures <- .smap_apply(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    # Build full argument list: first is structure, then other args
    args <- list(structures[[row$code]])
    if (length(.l) > 1) {
      for (j in 2:length(.l)) {
        args[[j]] <- row[[paste0("arg", j)]]
      }
    }
    result <- do.call(.f, c(args, list(...)))
    if (!inherits(result, "igraph")) {
      cli::cli_abort("Function `.f` must return an igraph object when using `spmap_structure()`.")
    }
    result
  }, use_parallel = .parallel)
  
  # Rebuild glycan_structure with proper deduplication
  idx <- match(combinations_df$combo_key, unique_combinations_df$combo_key)
  .rebuild_structure_with_dedup(new_structures, idx)
}

#' Map Functions Over Glycan Structure Vectors with Indices
#'
#' @description
#' These functions apply a function to each unique structure in a glycan structure vector
#' along with their corresponding indices,
#' taking advantage of hash-based deduplication to avoid redundant computation.
#' Similar to purrr imap functions, but optimized for glycan structure vectors.
#'
#' @param .x A glycan structure vector (glyrepr_structure).
#' @param .f A function that takes an igraph object (from `.x`) and an index/name,
#'   returning a result. 
#'   Can be a function, purrr-style lambda (`~ paste(.x, .y)`), or a character string naming a function.
#' @param ... Additional arguments passed to `.f`.
#' @param .ptype A prototype for the return type (for `simap_vec`).
#'
#' @details
#' These functions only compute `.f` once for each unique combination of structure and corresponding
#' index/name, then map the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element individually when there are duplicate structures.
#'
#' **IMPORTANT PERFORMANCE NOTE:**
#' Due to the inclusion of position indices, `simap` functions have **O(total_structures)** 
#' time complexity because each position creates a unique combination, even with identical structures.
#'
#' **Alternative:** Consider `smap()` functions if position information is not required.
#'
#' The index passed to `.f` is the position in the original vector (1-based).
#' If the vector has names, the names are passed instead of indices.
#'
#' **Return Types:**
#' - `simap()`: Returns a list with the same length as `.x`
#' - `simap_vec()`: Returns an atomic vector with the same length as `.x`
#' - `simap_lgl()`: Returns a logical vector
#' - `simap_int()`: Returns an integer vector  
#' - `simap_dbl()`: Returns a double vector
#' - `simap_chr()`: Returns a character vector
#' - `simap_structure()`: Returns a new glycan structure vector (`.f` must return igraph objects)
#'
#' @return 
#' - `simap()`: A list
#' - `simap_vec()`: An atomic vector of type specified by `.ptype`
#' - `simap_lgl()`: Returns a logical vector
#' - `simap_int()`: Returns an integer vector  
#' - `simap_dbl()`: Returns a double vector
#' - `simap_chr()`: Returns a character vector
#' - `simap_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create structure vectors with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
#' 
#' # Map a function that uses both structure and index
#' simap_chr(structures, function(g, i) paste0("Structure_", i, "_vcount_", igraph::vcount(g)))
#' 
#' # Use purrr-style lambda functions  
#' simap_chr(structures, ~ paste0("Pos", .y, "_vertices", igraph::vcount(.x)))
#' 
#' @name simap
NULL

#' @rdname simap
#' @export
simap <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(list())
  }
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Get indices or names
  if (!is.null(names(.x))) {
    indices <- names(.x)
  } else {
    indices <- seq_along(.x)
  }
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    index = indices,
    stringsAsFactors = FALSE
  )
  combinations_df$combo_key <- paste0(combinations_df$code, "|||", combinations_df$index)
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  unique_results <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    .f(structures[[row$code]], row$index, ...)
  })
  names(unique_results) <- unique_combinations_df$combo_key
  
  # Map results back to original vector positions
  purrr::map(combinations_df$combo_key, ~ unique_results[[.x]])
}

#' @rdname simap
#' @export
simap_vec <- function(.x, .f, ..., .ptype = NULL) {
  results <- simap(.x, .f, ...)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname simap
#' @export
simap_lgl <- function(.x, .f, ...) {
  simap_vec(.x, .f, ..., .ptype = logical())
}

#' @rdname simap
#' @export
simap_int <- function(.x, .f, ...) {
  simap_vec(.x, .f, ..., .ptype = integer())
}

#' @rdname simap
#' @export
simap_dbl <- function(.x, .f, ...) {
  simap_vec(.x, .f, ..., .ptype = double())
}

#' @rdname simap
#' @export
simap_chr <- function(.x, .f, ...) {
  simap_vec(.x, .f, ..., .ptype = character())
}

#' @rdname simap
#' @export
simap_structure <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(glycan_structure())
  }
  
  .f <- rlang::as_function(.f)
  
  data <- vctrs::vec_data(.x)
  codes <- vctrs::field(data, "iupac")
  structures <- attr(.x, "structures")
  
  # Get indices or names
  if (!is.null(names(.x))) {
    indices <- names(.x)
  } else {
    indices <- seq_along(.x)
  }
  
  # Create unique combinations data frame for proper handling
  combinations_df <- data.frame(
    code = codes,
    index = indices,
    stringsAsFactors = FALSE
  )
  combinations_df$combo_key <- paste0(combinations_df$code, "|||", combinations_df$index)
  
  unique_combinations_df <- combinations_df[!duplicated(combinations_df$combo_key), ]
  
  # Apply function only to unique combinations
  new_structures <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    result <- .f(structures[[row$code]], row$index, ...)
    if (!inherits(result, "igraph")) {
      cli::cli_abort("Function `.f` must return an igraph object when using `simap_structure()`.")
    }
    result
  })
  
  # Rebuild glycan_structure with proper deduplication
  idx <- match(combinations_df$combo_key, unique_combinations_df$combo_key)
  .rebuild_structure_with_dedup(new_structures, idx)
}
