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
#'
#' @details
#' These functions only compute `.f` once for each unique structure, then map
#' the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element individually when there are duplicate structures.
#'
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
#' @name smap
NULL

#' @rdname smap
#' @export
smap <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  .f <- rlang::as_function(.f)
  
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

#' @rdname smap
#' @export
smap_vec <- function(.x, .f, ..., .ptype = NULL) {
  results <- smap(.x, .f, ...)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname smap
#' @export
smap_lgl <- function(.x, .f, ...) {
  smap_vec(.x, .f, ..., .ptype = logical())
}

#' @rdname smap
#' @export
smap_int <- function(.x, .f, ...) {
  smap_vec(.x, .f, ..., .ptype = integer())
}

#' @rdname smap
#' @export
smap_dbl <- function(.x, .f, ...) {
  smap_vec(.x, .f, ..., .ptype = double())
}

#' @rdname smap
#' @export
smap_chr <- function(.x, .f, ...) {
  smap_vec(.x, .f, ..., .ptype = character())
}

#' @rdname smap
#' @export
smap_structure <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.x)
  structures <- attr(.x, "structures")
  
  # Apply function only to unique structures
  unique_codes <- names(structures)
  transformed_structures <- purrr::map(unique_codes, function(code) {
    result <- .f(structures[[code]], ...)
    if (!inherits(result, "igraph")) {
      rlang::abort("Function `.f` must return an igraph object when using `smap_structure()`.")
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
#'   Can be a function, purrr-style lambda (`~ .x$attr`), or a character string naming a function.
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
#' unique_results <- smap_unique(structures, igraph::vcount)
#' length(unique_results)  # 1, not 3
#' 
#' # Use purrr-style lambda
#' unique_results2 <- smap_unique(structures, ~ igraph::vcount(.x))
#' length(unique_results2)  # 1, not 3
#'
#' @export
smap_unique <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input must be a glycan_structure vector.")
  }
  
  .f <- rlang::as_function(.f)
  
  structures <- attr(.x, "structures")
  
  # Apply function only to unique structures
  results <- purrr::map(structures, .f, ...)
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
    rlang::abort("Input must be a glycan_structure vector.")
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
    rlang::abort("Input must be a glycan_structure vector.")
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
    rlang::abort("Input must be a glycan_structure vector.")
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
#'
#' @details
#' These functions only compute `.f` once for each unique combination of structure and corresponding
#' `.y` value, then map the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element pair individually when there are duplicate structure-value combinations.
#'
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
smap2 <- function(.x, .y, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(list())
  }
  
  # Recycle .y to match length of .x
  .y <- vctrs::vec_recycle(.y, length(.x))
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.x)
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
  unique_results <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    .f(structures[[row$code]], row$y_val, ...)
  })
  names(unique_results) <- unique_combinations_df$combo_key
  
  # Map results back to original vector positions
  purrr::map(combinations_df$combo_key, ~ unique_results[[.x]])
}

#' @rdname smap2
#' @export
smap2_vec <- function(.x, .y, .f, ..., .ptype = NULL) {
  results <- smap2(.x, .y, .f, ...)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname smap2
#' @export
smap2_lgl <- function(.x, .y, .f, ...) {
  smap2_vec(.x, .y, .f, ..., .ptype = logical())
}

#' @rdname smap2
#' @export
smap2_int <- function(.x, .y, .f, ...) {
  smap2_vec(.x, .y, .f, ..., .ptype = integer())
}

#' @rdname smap2
#' @export
smap2_dbl <- function(.x, .y, .f, ...) {
  smap2_vec(.x, .y, .f, ..., .ptype = double())
}

#' @rdname smap2
#' @export
smap2_chr <- function(.x, .y, .f, ...) {
  smap2_vec(.x, .y, .f, ..., .ptype = character())
}

#' @rdname smap2
#' @export
smap2_structure <- function(.x, .y, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(glycan_structure())
  }
  
  # Recycle .y to match length of .x
  .y <- vctrs::vec_recycle(.y, length(.x))
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.x)
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
  transformed_structures <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    result <- .f(structures[[row$code]], row$y_val, ...)
    if (!inherits(result, "igraph")) {
      rlang::abort("Function `.f` must return an igraph object when using `smap2_structure()`.")
    }
    result
  })
  names(transformed_structures) <- unique_combinations_df$combo_key
  
  # Create new glycan structure vector from transformed structures
  # Map back to original positions
  individual_structures <- purrr::map(combinations_df$combo_key, ~ transformed_structures[[.x]])
  
  do.call(glycan_structure, individual_structures)
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
#'
#' @details
#' These functions only compute `.f` once for each unique combination of structure and corresponding
#' values from other vectors, then map the results back to the original vector positions. This is much more efficient
#' than applying `.f` to each element combination individually when there are duplicate combinations.
#'
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
spmap <- function(.l, .f, ...) {
  if (!is.list(.l) || length(.l) == 0) {
    rlang::abort("Input `.l` must be a non-empty list.")
  }
  
  if (!is_glycan_structure(.l[[1]])) {
    rlang::abort("First element of `.l` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.l[[1]]) == 0) {
    return(list())
  }
  
  # Recycle all vectors to match length of first vector
  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.l[[1]])
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
  unique_results <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    # Build full argument list: first is structure, then other args
    args <- list(structures[[row$code]])
    if (length(.l) > 1) {
      for (j in 2:length(.l)) {
        args[[j]] <- row[[paste0("arg", j)]]
      }
    }
    do.call(.f, c(args, list(...)))
  })
  names(unique_results) <- unique_combinations_df$combo_key
  
  # Map results back to original vector positions
  purrr::map(combinations_df$combo_key, ~ unique_results[[.x]])
}

#' @rdname spmap
#' @export
spmap_vec <- function(.l, .f, ..., .ptype = NULL) {
  results <- spmap(.l, .f, ...)
  vctrs::vec_c(!!!results, .ptype = .ptype)
}

#' @rdname spmap
#' @export
spmap_lgl <- function(.l, .f, ...) {
  spmap_vec(.l, .f, ..., .ptype = logical())
}

#' @rdname spmap
#' @export
spmap_int <- function(.l, .f, ...) {
  spmap_vec(.l, .f, ..., .ptype = integer())
}

#' @rdname spmap
#' @export
spmap_dbl <- function(.l, .f, ...) {
  spmap_vec(.l, .f, ..., .ptype = double())
}

#' @rdname spmap
#' @export
spmap_chr <- function(.l, .f, ...) {
  spmap_vec(.l, .f, ..., .ptype = character())
}

#' @rdname spmap
#' @export
spmap_structure <- function(.l, .f, ...) {
  if (!is.list(.l) || length(.l) == 0) {
    rlang::abort("Input `.l` must be a non-empty list.")
  }
  
  if (!is_glycan_structure(.l[[1]])) {
    rlang::abort("First element of `.l` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.l[[1]]) == 0) {
    return(glycan_structure())
  }
  
  # Recycle all vectors to match length of first vector
  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.l[[1]])
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
  transformed_structures <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
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
      rlang::abort("Function `.f` must return an igraph object when using `spmap_structure()`.")
    }
    result
  })
  names(transformed_structures) <- unique_combinations_df$combo_key
  
  # Create new glycan structure vector from transformed structures
  # Map back to original positions
  individual_structures <- purrr::map(combinations_df$combo_key, ~ transformed_structures[[.x]])
  
  do.call(glycan_structure, individual_structures)
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
#' The index passed to `.f` is the position in the original vector (1-based).
#' If the vector has names, the names are passed instead of indices.
#'
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
#' - `simap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
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
#' # With named vectors
#' names(structures) <- c("first", "second", "third")
#' simap_chr(structures, ~ paste0(.y, "_has_", igraph::vcount(.x), "_vertices"))
#'
#' @name simap
NULL

#' @rdname simap
#' @export
simap <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    rlang::abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(list())
  }
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.x)
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
    rlang::abort("Input `.x` must be a glycan_structure vector.")
  }
  
  # Handle empty input
  if (length(.x) == 0) {
    return(glycan_structure())
  }
  
  .f <- rlang::as_function(.f)
  
  codes <- vctrs::vec_data(.x)
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
  transformed_structures <- purrr::map(seq_len(nrow(unique_combinations_df)), function(i) {
    row <- unique_combinations_df[i, ]
    result <- .f(structures[[row$code]], row$index, ...)
    if (!inherits(result, "igraph")) {
      rlang::abort("Function `.f` must return an igraph object when using `simap_structure()`.")
    }
    result
  })
  names(transformed_structures) <- unique_combinations_df$combo_key
  
  # Create new glycan structure vector from transformed structures
  # Map back to original positions
  individual_structures <- purrr::map(combinations_df$combo_key, ~ transformed_structures[[.x]])
  
  do.call(glycan_structure, individual_structures)
}
