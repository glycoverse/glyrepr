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
#' @returns
#' - `smap()`: A list
#' - `smap_vec()`: An atomic vector of type specified by `.ptype`
#' - `smap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type
#' - `smap_structure()`: A new glyrepr_structure object
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' core2 <- n_glycan_core()
#' structures <- c(core1, core2, core1)  # core1 appears twice
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

# Helper function to rebuild glycan_structure with proper deduplication
# after modifications that may create identical graphs
.rebuild_structure_with_dedup <- function(modified_graphs, idx_mapping) {
  # Get new IUPACs for all modified graphs
  new_unique_iupacs <- purrr::map_chr(
    modified_graphs,
    graph_to_iupac
  )
  new_iupacs <- new_unique_iupacs[idx_mapping]

  # Re-deduplicate graphs based on new IUPACs to handle cases where
  # modifications create identical graphs
  unique_new_indices <- which(!duplicated(new_unique_iupacs))
  final_unique_graphs <- modified_graphs[unique_new_indices]
  final_unique_iupacs <- new_unique_iupacs[unique_new_indices]
  names(final_unique_graphs) <- final_unique_iupacs

  # Create result glycan_structure
  new_glycan_structure(new_iupacs, final_unique_graphs)
}

#' Split a glycan structure vector into valid and missing positions
#'
#' @param x A glycan structure vector.
#' @returns A list with vector data, graph data, names, and NA masks.
#' @noRd
.structure_map_input <- function(x) {
  codes <- vctrs::vec_data(x)
  na_mask <- is.na(codes)

  list(
    codes = codes,
    graphs = attr(x, "graphs"),
    input_names = names(x),
    na_mask = na_mask,
    na_count = sum(na_mask),
    has_na = any(na_mask),
    all_na = length(codes) > 0 && all(na_mask),
    valid_codes = codes[!na_mask],
    valid_names = names(x)[!na_mask]
  )
}

#' Build a deduplicated-combination table for structure mapping
#'
#' @param valid_codes Non-missing structure codes.
#' @param args Additional recycled argument vectors aligned with `valid_codes`.
#' @returns A tibble with a `code`, one `arg*` column per argument, and `combo_key`.
#' @noRd
.structure_combo_table <- function(valid_codes, args = list()) {
  combinations_df <- tibble::tibble(code = valid_codes)
  key_components <- list(valid_codes)

  for (i in seq_along(args)) {
    arg_name <- paste0("arg", i)
    combinations_df[[arg_name]] <- args[[i]]
    key_components <- append(
      key_components,
      list(purrr::map_chr(args[[i]], .generate_value_key))
    )
  }

  combinations_df$combo_key <- do.call(paste, c(key_components, sep = "|||"))
  combinations_df
}

#' Return unique rows from a structure combination table
#'
#' @param combinations_df A table returned by `.structure_combo_table()`.
#' @returns A tibble containing unique combination keys.
#' @noRd
.unique_structure_combos <- function(combinations_df) {
  combinations_df[!duplicated(combinations_df$combo_key), ]
}

#' Extract mapped arguments from one combination-table row
#'
#' @param row A one-row tibble from `.structure_combo_table()`.
#' @param n_args Number of argument columns to extract.
#' @returns A list of mapped argument values.
#' @noRd
.extract_combo_args <- function(row, n_args) {
  purrr::map(seq_len(n_args), function(i) {
    .extract_from_list_column(row[[paste0("arg", i)]])
  })
}

#' Restore a list result to the original structure vector shape
#'
#' @param valid_results Results for non-missing positions.
#' @param map_input Metadata from `.structure_map_input()`.
#' @param na_value Value to place at missing positions.
#' @returns A list with original length and names.
#' @noRd
.restore_list_with_na <- function(valid_results, map_input, na_value = NA) {
  if (!map_input$has_na) {
    names(valid_results) <- map_input$input_names
    return(valid_results)
  }

  result <- vector("list", length(map_input$codes))
  result[!map_input$na_mask] <- valid_results
  result[map_input$na_mask] <- list(na_value)
  names(result) <- map_input$input_names
  result
}

#' Restore a structure result to the original structure vector shape
#'
#' @param valid_result Result for non-missing positions.
#' @param map_input Metadata from `.structure_map_input()`.
#' @returns A glycan structure vector with original length and names.
#' @noRd
.restore_structure_with_na <- function(valid_result, map_input) {
  if (map_input$all_na) {
    result <- new_na_glycan_structure(length(map_input$codes))
    names(result) <- map_input$input_names
    return(result)
  }

  if (map_input$has_na) {
    result_iupacs <- rep(NA_character_, length(map_input$codes))
    result_iupacs[!map_input$na_mask] <- vctrs::vec_data(valid_result)

    result <- new_glycan_structure(result_iupacs, attr(valid_result, "graphs"))
    names(result) <- map_input$input_names
    return(result)
  }

  names(valid_result) <- map_input$input_names
  valid_result
}

#' Map unique non-missing structure and argument combinations
#'
#' @param map_input Metadata from `.structure_map_input()`.
#' @param args Additional argument vectors aligned to non-missing structures.
#' @param .f A function called with one graph followed by values from `args`.
#' @param dots Additional arguments captured from `...`.
#' @param .structure Logical; whether `.f` must return igraph objects.
#' @param .caller Public function name used in graph-return validation messages.
#' @returns A list or glycan structure vector restored to the input shape.
#' @noRd
.map_structure_combinations <- function(
  map_input,
  args = list(),
  .f,
  dots = list(),
  .structure = FALSE,
  .caller = "map function"
) {
  if (map_input$all_na) {
    if (.structure) {
      return(.restore_structure_with_na(NULL, map_input))
    }

    result <- vector("list", length(map_input$codes))
    result[map_input$na_mask] <- list(NA)
    names(result) <- map_input$input_names
    return(result)
  }

  combinations_df <- .structure_combo_table(map_input$valid_codes, args)
  unique_combinations_df <- .unique_structure_combos(combinations_df)
  n_args <- length(args)

  unique_results <- purrr::map(
    seq_len(nrow(unique_combinations_df)),
    function(i) {
      row <- unique_combinations_df[i, ]
      mapped_args <- .extract_combo_args(row, n_args)
      result <- do.call(
        .f,
        c(list(map_input$graphs[[row$code]]), mapped_args, dots)
      )

      if (.structure && !inherits(result, "igraph")) {
        cli::cli_abort(paste0(
          "Function `.f` must return an igraph object when using `",
          .caller,
          "`."
        ))
      }

      result
    }
  )
  names(unique_results) <- unique_combinations_df$combo_key

  if (.structure) {
    idx <- match(combinations_df$combo_key, unique_combinations_df$combo_key)
    valid_result <- .rebuild_structure_with_dedup(unique_results, idx)
    names(valid_result) <- map_input$valid_names
    return(.restore_structure_with_na(valid_result, map_input))
  }

  valid_results <- purrr::map(
    combinations_df$combo_key,
    ~ unique_results[[.x]]
  )
  names(valid_results) <- map_input$valid_names

  .restore_list_with_na(valid_results, map_input, na_value = NA)
}

# Helper function for common smap logic
.smap_base <- function(.x, .f, ..., .convert_fn = NULL) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .f <- rlang::as_function(.f)

  map_input <- .structure_map_input(.x)

  # If all elements are NA, return NA results
  if (map_input$all_na) {
    if (is.null(.convert_fn)) {
      result <- vector("list", length(map_input$codes))
    } else {
      # Get typed NA by running .convert_fn on a properly-typed NA
      # Use NA_character_ for character output (as.character(NA) gives "NA" string)
      # For other types, NA coerces correctly
      test_na <- NA_character_ # Default to character type
      typed_na <- .convert_fn(list(test_na))
      result <- rep(typed_na, map_input$na_count)
    }
    names(result) <- map_input$input_names
    return(result)
  }

  dots <- list(...)

  # Apply function only to unique graphs (NA already filtered out)
  unique_codes <- names(map_input$graphs)
  unique_results <- purrr::map(
    unique_codes,
    function(code) {
      do.call(.f, c(list(map_input$graphs[[code]]), dots))
    }
  )
  names(unique_results) <- unique_codes

  # If no conversion function provided, return list (for smap)
  if (is.null(.convert_fn)) {
    # Optimized mapping: use match() instead of individual lookups
    idx <- match(map_input$valid_codes, unique_codes)
    valid_results <- unique_results[idx]
    names(valid_results) <- map_input$valid_names

    return(.restore_list_with_na(valid_results, map_input, na_value = NULL))
  }

  # Convert to target type and map back to original positions
  unique_converted <- .convert_fn(unique_results)
  names(unique_converted) <- unique_codes

  # Optimized mapping: use match() instead of individual lookups
  idx <- match(map_input$valid_codes, unique_codes)
  valid_result <- unique_converted[idx]

  # Initialize result with NA and assign valid results
  result <- rep(NA, length(map_input$codes))
  result[!map_input$na_mask] <- valid_result
  names(result) <- map_input$input_names
  return(result)
}

# Helper function to generate a hash-based key for complex objects
# Used for deduplication when dealing with nested lists and other complex structures
.generate_object_hash <- function(obj) {
  serialized_data <- serialize(obj, connection = NULL)
  hash_input <- sum(as.integer(serialized_data))
  hex_hash <- as.hexmode(hash_input)
  format(hex_hash, width = 8)
}

# Helper function to extract values from tibble list-columns
# Tibble stores single-element lists as list(data) in list-columns, requiring unwrapping
.extract_from_list_column <- function(value) {
  if (is.list(value) && length(value) == 1) {
    value[[1]] # Extract from list-column wrapper
  } else {
    value # Use as-is for non-list or multi-element cases
  }
}

# Helper function to generate deduplication keys for arbitrary values
# Handles both simple values (convert to character) and complex objects (hash-based keys)
.generate_value_key <- function(value, prefix = "list_") {
  if (is.list(value)) {
    paste0(prefix, .generate_object_hash(value))
  } else {
    as.character(value)
  }
}

#' @rdname smap
#' @export
smap <- function(.x, .f, ...) {
  .smap_base(.x, .f, ..., .convert_fn = NULL)
}

#' @rdname smap
#' @export
smap_vec <- function(.x, .f, ..., .ptype = NULL) {
  .smap_base(
    .x,
    .f,
    ...,
    .convert_fn = function(results) vctrs::vec_c(!!!results, .ptype = .ptype)
  )
}

#' @rdname smap
#' @export
smap_lgl <- function(.x, .f, ...) {
  .smap_base(
    .x,
    .f,
    ...,
    .convert_fn = function(results) {
      as.logical(unlist(results, use.names = FALSE))
    }
  )
}

#' @rdname smap
#' @export
smap_int <- function(.x, .f, ...) {
  .smap_base(
    .x,
    .f,
    ...,
    .convert_fn = function(results) {
      as.integer(unlist(results, use.names = FALSE))
    }
  )
}

#' @rdname smap
#' @export
smap_dbl <- function(.x, .f, ...) {
  .smap_base(
    .x,
    .f,
    ...,
    .convert_fn = function(results) {
      as.double(unlist(results, use.names = FALSE))
    }
  )
}

#' @rdname smap
#' @export
smap_chr <- function(.x, .f, ...) {
  .smap_base(
    .x,
    .f,
    ...,
    .convert_fn = function(results) {
      as.character(unlist(results, use.names = FALSE))
    }
  )
}

#' @rdname smap
#' @export
smap_structure <- function(.x, .f, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .f <- rlang::as_function(.f)
  map_input <- .structure_map_input(.x)

  dots <- list(...)

  if (map_input$all_na) {
    return(.restore_structure_with_na(NULL, map_input))
  }

  # Apply function only to unique graphs
  unique_iupacs <- names(map_input$graphs)
  new_graphs <- purrr::map(
    unique_iupacs,
    function(iupac) {
      result <- do.call(.f, c(list(map_input$graphs[[iupac]]), dots))
      if (!inherits(result, "igraph")) {
        cli::cli_abort(
          "Function `.f` must return an igraph object when using `smap_structure()`."
        )
      }
      result
    }
  )

  # Rebuild glycan_structure with proper deduplication
  idx <- match(map_input$valid_codes, unique_iupacs)
  valid_result <- .rebuild_structure_with_dedup(new_graphs, idx)
  names(valid_result) <- map_input$valid_names

  .restore_structure_with_na(valid_result, map_input)
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
#' @return A list with results for each unique structure, named by their hash codes.
#'
#' @examples
#' # Create a structure vector with duplicates
#' core1 <- o_glycan_core_1()
#' structures <- c(core1, core1, core1)  # same structure 3 times
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
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .f <- rlang::as_function(.f)

  graphs <- attr(.x, "graphs")

  dots <- list(...)

  # Apply function only to unique graphs
  results <- purrr::map(
    graphs,
    function(g) {
      do.call(.f, c(list(g), dots))
    }
  )
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
#' structures <- c(core1, core2, core1)  # core1 appears twice
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

  graphs <- attr(.x, "graphs")

  # Apply predicate only to unique graphs using purrr::some
  purrr::some(graphs, .p, ...)
}

#' @rdname smap_predicates
#' @export
severy <- function(.x, .p, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .p <- rlang::as_function(.p)

  graphs <- attr(.x, "graphs")

  # Apply predicate only to unique graphs using purrr::every
  purrr::every(graphs, .p, ...)
}

#' @rdname smap_predicates
#' @export
snone <- function(.x, .p, ...) {
  if (!is_glycan_structure(.x)) {
    cli::cli_abort("Input must be a glycan_structure vector.")
  }

  .p <- rlang::as_function(.p)

  graphs <- attr(.x, "graphs")

  # Apply predicate only to unique graphs using purrr::none
  purrr::none(graphs, .p, ...)
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
#' **NA Handling:**
#' NA elements in `.x` are preserved in the output - the function is not applied to NA positions,
#' and the corresponding results are set to NA.
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
#' structures <- c(core1, core2, core1)  # core1 appears twice
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
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }

  if (length(.x) == 0) {
    return(list())
  }

  .y <- vctrs::vec_recycle(.y, length(.x))
  .f <- rlang::as_function(.f)

  map_input <- .structure_map_input(.x)
  valid_y <- .y[!map_input$na_mask]

  .map_structure_combinations(
    map_input,
    list(valid_y),
    .f,
    dots = list(...)
  )
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
    cli::cli_abort("Input `.x` must be a glycan_structure vector.")
  }

  if (length(.x) == 0) {
    return(glycan_structure())
  }

  .y <- vctrs::vec_recycle(.y, length(.x))
  .f <- rlang::as_function(.f)

  map_input <- .structure_map_input(.x)
  valid_y <- .y[!map_input$na_mask]

  .map_structure_combinations(
    map_input,
    list(valid_y),
    .f,
    dots = list(...),
    .structure = TRUE,
    .caller = "smap2_structure()"
  )
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
#' values from other vectors, then map the results back to the original vector positions.
#'
#' **NA Handling:**
#' NA elements in the first argument (glycan structure vector) are preserved in the output.
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
#' structures <- c(core1, core2, core1)  # core1 appears twice
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
  if (!inherits(.l, "list") || inherits(.l, "vctrs_vctr") || length(.l) == 0) {
    cli::cli_abort("Input `.l` must be a non-empty list.")
  }

  if (!is_glycan_structure(.l[[1]])) {
    cli::cli_abort("First element of `.l` must be a glycan_structure vector.")
  }

  if (length(.l[[1]]) == 0) {
    return(list())
  }

  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  .f <- rlang::as_function(.f)

  map_input <- .structure_map_input(.l[[1]])
  valid_args <- purrr::map(.l[-1], ~ .x[!map_input$na_mask])

  .map_structure_combinations(
    map_input,
    valid_args,
    .f,
    dots = list(...)
  )
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
    cli::cli_abort("Input `.l` must be a non-empty list.")
  }

  if (!is_glycan_structure(.l[[1]])) {
    cli::cli_abort("First element of `.l` must be a glycan_structure vector.")
  }

  if (length(.l[[1]]) == 0) {
    return(glycan_structure())
  }

  target_length <- length(.l[[1]])
  .l <- purrr::map(.l, ~ vctrs::vec_recycle(.x, target_length))
  .f <- rlang::as_function(.f)

  map_input <- .structure_map_input(.l[[1]])
  valid_args <- purrr::map(.l[-1], ~ .x[!map_input$na_mask])

  .map_structure_combinations(
    map_input,
    valid_args,
    .f,
    dots = list(...),
    .structure = TRUE,
    .caller = "spmap_structure()"
  )
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
#' structures <- c(core1, core2, core1)  # core1 appears twice
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
  map_input <- .structure_map_input(.x)

  # Get indices or names
  if (!is.null(names(.x))) {
    indices <- names(.x)
  } else {
    indices <- seq_along(.x)
  }

  valid_indices <- indices[!map_input$na_mask]

  .map_structure_combinations(
    map_input,
    list(valid_indices),
    .f,
    dots = list(...)
  )
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
  map_input <- .structure_map_input(.x)

  # Get indices or names
  if (!is.null(names(.x))) {
    indices <- names(.x)
  } else {
    indices <- seq_along(.x)
  }

  valid_indices <- indices[!map_input$na_mask]

  .map_structure_combinations(
    map_input,
    list(valid_indices),
    .f,
    dots = list(...),
    .structure = TRUE,
    .caller = "simap_structure()"
  )
}
