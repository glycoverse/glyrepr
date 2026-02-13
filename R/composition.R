#' Create a Glycan Composition
#'
#' Create a glycan composition from a list of named integer vectors.
#' Compositions can contain both monosaccharides and substituents.
#'
#' @param ... Named integer vectors.
#'   Names are monosaccharides or substituents, values are numbers of residues.
#'   Monosaccharides and substituents can be mixed in the same composition.
#' @param x A list of named integer vectors.
#'
#' @returns A glyrepr_composition object.
#'
#' @details
#' Compositions can contain:
#'
#' - Monosaccharides: either generic (e.g., "Hex", "HexNAc") or concrete (e.g., "Glc", "Gal").
#'   All monosaccharides in a composition vector must be of the same type.
#' - Substituents: e.g., "Me", "Ac", "S". These can be mixed with either
#'   generic or concrete monosaccharides.
#'
#' Components are automatically sorted with monosaccharides first (according to
#' their order in the monosaccharides table), followed by substituents (according
#' to their order in `available_substituents()`).
#'
#' @examples
#' # A vector with one composition (generic monosaccharides)
#' glycan_composition(c(Hex = 5, HexNAc = 2))
#' # A vector with multiple compositions
#' glycan_composition(c(Hex = 5, HexNAc = 2), c(Hex = 5, HexNAc = 4, dHex = 2))
#' # Residues are reordered automatically
#' glycan_composition(c(HexNAc = 1, Hex = 2))
#' # An example for generic monosaccharides
#' glycan_composition(c(Hex = 2, HexNAc = 1))
#' # An example for concrete monosaccharides
#' glycan_composition(c(Glc = 2, Gal = 1))
#' # Compositions with substituents
#' glycan_composition(c(Glc = 1, S = 1))
#' glycan_composition(c(Hex = 3, HexNAc = 2, Me = 1, Ac = 1))
#' # Substituents are sorted after monosaccharides
#' glycan_composition(c(S = 1, Gal = 1, Ac = 1, Glc = 1))
#'
#' @seealso `available_monosaccharides()`, `available_substituents()`
#' @export
glycan_composition <- function(...) {
  args <- rlang::list2(...)

  # Handle NULL/NA inputs - convert to list elements that will be stored as NULL
  # Only unnamed NA/NULL are treated as missing; named NA (e.g., c(Hex = NA)) is invalid
  x <- purrr::map(args, ~ {
    if (is.null(.x)) {
      NULL  # Store as NULL to represent NA
    } else if (is.atomic(.x) && length(.x) == 1 && is.na(.x) && is.null(names(.x))) {
      NULL  # Unnamed scalar NA treated as missing
    } else {
      .valid_composition_element(.x)
    }
  })

  # Process non-NA elements
  x <- purrr::map(x, ~ {
    if (is.null(.x)) return(.x)
    result <- as.integer(.x)
    names(result) <- names(.x)
    .reorder_composition_components(result)
  })

  .valid_glycan_composition_input(x)
  new_glycan_composition(x)
}

#' @export
#' @rdname glycan_composition
is_glycan_composition <- function(x) {
  inherits(x, "glyrepr_composition")
}

#' Convert to Glycan Composition
#'
#' Converts an object to a glycan composition using vctrs casting framework.
#' This function provides a convenient way to convert various input types
#' to [glycan_composition()].
#'
#' @param x An object to convert to a glycan composition.
#'   Supported inputs include:
#'   - Named integer vectors or lists of named integer vectors
#'   - Character vectors with composition strings (e.g., "Hex(5)HexNAc(2)")
#'   - `glyrepr_structure` objects (counts both monosaccharides and substituents)
#'   - Existing `glyrepr_composition` objects (returned as-is)
#'
#' @returns A `glyrepr_composition` object.
#'
#' @details
#' This function uses the vctrs casting framework for type conversion.
#' When converting from glycan structures, both monosaccharides and substituents
#' are counted. Substituents are extracted from the `sub` attribute of each
#' vertex in the structure. For example, a vertex with `sub = "3Me"`
#' contributes one "Me" substituent to the composition.
#'
#' @examples
#' # From a single named vector
#' as_glycan_composition(c(Hex = 5, HexNAc = 2))
#'
#' # From a list of named vectors
#' as_glycan_composition(list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1)))
#'
#' # From a character vector of Byonic composition strings
#' as_glycan_composition(c("Hex(5)HexNAc(2)", "Hex(3)HexNAc(1)"))
#'
#' # From a character vector of simple composition strings
#' as_glycan_composition(c("H5N2", "H5N4S1F1"))
#'
#' # From an existing composition (returns as-is)
#' comp <- glycan_composition(c(Hex = 5, HexNAc = 2))
#' as_glycan_composition(comp)
#'
#' # From a glycan structure vector
#' strucs <- c(n_glycan_core(), o_glycan_core_1())
#' as_glycan_composition(strucs)
#'
#' @export
as_glycan_composition <- function(x) {
  vec_cast(x, new_glycan_composition())
}

#' @export
vec_ptype_full.glyrepr_composition <- function(x, ...) "glycan_composition"

#' @export
vec_ptype_abbr.glyrepr_composition <- function(x, ...) "comp"

#' @export
format.glyrepr_composition <- function(x, ...) {
  data <- vctrs::field(vctrs::vec_data(x), "data")

  formatted <- purrr::map_chr(data, ~ {
    if (.is_na_composition_elem(.x)) {
      "<NA>"
    } else {
      paste0(names(.x), "(", .x, ")", collapse = "")
    }
  })

  formatted
}

#' @export
vec_ptype2.glyrepr_composition.glyrepr_composition <- function(x, y, ...) {
  new_glycan_composition(list())
}

#' @export
vec_cast.glyrepr_composition.glyrepr_composition <- function(x, to, ...) {
  x
}

#' @export
vec_cast.character.glyrepr_composition <- function(x, to, ...) {
  convert_one <- function(comp) {
    if (.is_na_composition_elem(comp)) {
      "<NA>"
    } else {
      paste0(names(comp), "(", comp, ")", collapse = "")
    }
  }
  data <- vctrs::vec_data(x)
  purrr::map_chr(vctrs::field(data, "data"), convert_one)
}

#' @export
vec_cast.glyrepr_composition.character <- function(x, to, ...) {
  # Handle empty character vector
  if (length(x) == 0) {
    return(glycan_composition())
  }

  # Handle NA by creating NA composition elements
  na_mask <- is.na(x)
  non_na_positions <- which(!na_mask)
  if (any(na_mask)) {
    # Process non-NA characters
    non_na_x <- x[!na_mask]
    compositions <- list()

    if (length(non_na_x) > 0) {
      # Parse non-NA characters
      parse_result <- purrr::map(non_na_x, parse_single_composition)
      valid_flags <- purrr::map_lgl(parse_result, "valid")

      # Check for invalid characters
      invalid_indices <- which(!valid_flags)
      # Map back to original indices in x (not indices in non_na_x)
      original_invalid <- non_na_positions[invalid_indices]
      if (length(original_invalid) > 0) {
        cli::cli_abort(c(
          "Characters cannot be parsed as glycan compositions at index {original_invalid}",
          "i" = "Expected format: 'Hex(5)HexNAc(2)' with monosaccharide names followed by counts in parentheses."
        ))
      }

      compositions <- purrr::map(parse_result, "composition")
    }

    # Build result list preserving order - no global state needed
    result_list <- purrr::map(seq_along(x), ~ {
      if (na_mask[.x]) {
        NULL
      } else {
        # Find which composition this position corresponds to
        pos_in_non_na <- match(.x, non_na_positions)
        compositions[[pos_in_non_na]]
      }
    })

    return(do.call(glycan_composition, result_list))
  }

  # Original logic for non-NA case (unchanged)
  parse_result <- purrr::map(x, parse_single_composition)
  valid_flags <- purrr::map_lgl(parse_result, "valid")
  compositions <- purrr::map(parse_result, "composition")
  invalid_indices <- which(!valid_flags)
  if (length(invalid_indices) > 0) {
    cli::cli_abort(c(
      "Characters cannot be parsed as glycan compositions at index {invalid_indices}",
      "i" = "Expected format: 'Hex(5)HexNAc(2)' with monosaccharide names followed by counts in parentheses."
    ))
  }
  if (length(compositions) == 0) {
    return(glycan_composition())
  }
  do.call(glycan_composition, compositions)
}

#' @export
vec_cast.glyrepr_composition.glyrepr_structure <- function(x, to, ...) {
  # Use smap to convert each structure to composition
  compositions <- smap(x, function(graph) {
    # Count monosaccharides
    monos <- igraph::V(graph)$mono
    mono_tb <- table(monos)
    mono_result <- as.integer(mono_tb)
    names(mono_result) <- names(mono_tb)

    # Count substituents
    subs <- igraph::V(graph)$sub
    sub_types <- extract_substituent_types(subs)
    if (length(sub_types) > 0) {
      sub_tb <- table(sub_types)
      sub_result <- as.integer(sub_tb)
      names(sub_result) <- names(sub_tb)
    } else {
      sub_result <- integer(0)
    }

    # Combine monosaccharides and substituents
    result <- c(mono_result, sub_result)

    # Sort by composition component order (monosaccharides first, then substituents)
    result <- .reorder_composition_components(result)
    result
  })

  # Create composition object
  new_glycan_composition(compositions)
}

#' @export
vec_cast.glyrepr_composition.list <- function(x, to, ...) {
  if (length(x) == 0) {
    return(glycan_composition())
  }
  do.call(glycan_composition, x)
}

#' @export
vec_cast.glyrepr_composition.integer <- function(x, to, ...) {
  .cast_named_vector(x, to, ...)
}

#' @export
vec_cast.glyrepr_composition.double <- function(x, to, ...) {
  .cast_named_vector(x, to, ...)
}

#' @export
vec_ptype2.logical.glyrepr_composition <- function(x, y, ...) {
  new_glycan_composition(list())
}

#' @export
vec_ptype2.glyrepr_composition.logical <- function(x, y, ...) {
  new_glycan_composition(list())
}

#' @export
vec_cast.glyrepr_composition.logical <- function(x, to, ...) {
  # Cast each logical element to NA composition or error
  # TRUE and FALSE are invalid (not NA), only NA is valid
  result <- purrr::map(x, ~ {
    if (is.na(.x)) {
      NULL  # NA becomes NULL (NA composition)
    } else {
      cli::cli_abort("Cannot cast non-NA logical value to glyrepr_composition.")
    }
  })
  new_glycan_composition(result)
}

#' @export
vec_cast.logical.glyrepr_composition <- function(x, to, ...) {
  # Cast glyrepr_composition to logical:
  # - NA compositions become NA (preserving missingness)
  # - Valid compositions become FALSE
  is_na <- is.na(x)
  result <- rep(FALSE, length(x))
  result[is_na] <- NA
  result
}

#' @export
vec_restore.glyrepr_composition <- function(x, to, ...) {
  data <- vctrs::field(x, "data")

  # Skip NA elements (NULL) when checking types
  na_mask <- purrr::map_lgl(data, .is_na_composition_elem)
  non_na_data <- data[!na_mask]

  if (length(non_na_data) == 0) {
    return(new_glycan_composition(data))
  }

  monos_list <- purrr::map(non_na_data, ~ names(.x)[!names(.x) %in% available_substituents()])
  mono_types <- purrr::map_chr(monos_list, get_mono_type_impl)

  if (length(unique(mono_types)) > 1) {
    cli::cli_abort(c(
      "Can't combine `glyrepr_composition`s with different monosaccharide types.",
      "x" = "Found compositions with types: {.val {unique(mono_types)}}.",
      "i" = "Use {.fn convert_to_generic} to convert concrete to generic, or ensure both compositions use the same type."
    ))
  }

  new_glycan_composition(data)
}

#' @export
is.na.glyrepr_composition <- function(x, ...) {
  data <- vctrs::field(x, "data")
  purrr::map_lgl(data, .is_na_composition_elem)
}

.cast_named_vector <- function(x, to, ...) {
  if (length(x) == 0) {
    return(glycan_composition())
  }
  glycan_composition(x)
}

#' @export
obj_print_data.glyrepr_composition <- function(x, ..., max_n = 10, colored = TRUE) {
  if (length(x) == 0) {
    return()
  }

  n <- length(x)
  n_show <- min(n, max_n)

  # Only format the compositions that need to be shown to improve performance
  indices_to_show <- seq_len(n_show)
  formatted <- format_glycan_composition_subset(x, indices_to_show, colored = colored)

  # Print each composition on its own line with indexing, up to max_n
  for (i in seq_len(n_show)) {
    cat("[", i, "] ", formatted[i], "\n", sep = "")
  }
  if (n > max_n) {
    cat("... (", n - max_n, " more not shown)\n", sep = "")
  }
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.glyrepr_composition <- function(x, ...) {
  if (length(x) == 0) {
    return(pillar::pillar_shaft(character()))
  }

  indices <- seq_len(length(x))
  formatted <- format_glycan_composition_subset(x, indices, colored = TRUE)

  pillar::new_pillar_shaft_simple(formatted, align = "left", min_width = 10)
}

# --- Internal functions -------------------------------------------------------

#' Create a glycan composition object
#' @param x A list of named integer vectors.
#' @returns A glyrepr_composition object.
#' @noRd
new_glycan_composition <- function(x = list()) {
  # Use vctrs::list_of instead of new_list_of
  x <- vctrs::new_list_of(x, .ptype = integer(), class = "glyrepr_composition_list")
  vctrs::new_rcrd(list(data = x), class = "glyrepr_composition")
}

#' Validate glycan composition input
#'
#' Use by `glycan_composition()` to validate the input.
#'
#' @param x A list of named integer vectors.
#' @noRd
.valid_glycan_composition_input <- function(x) {
  # 0. Skip if empty
  if (length(x) == 0) {
    return()
  }

  # Filter out NA elements (NULL in list) for validation
  x_valid <- purrr::keep(x, ~ !is.null(.x))

  # 0b. If all elements are NA (NULL), skip further validation
  if (length(x_valid) == 0) {
    return()
  }

  # 1. Type check (skip NA elements)
  if (!purrr::every(x_valid, checkmate::test_integerish)) {
    cli::cli_abort(c(
      "Must be one or more named integer vectors.",
      "i" = "You might want to use {.fn as_glycan_composition} for more flexible input."
    ))
  }

  # 2. Empty check (skip NA elements)
  if (purrr::some(x_valid, ~ length(.x) == 0)) {
    cli::cli_abort("Each composition must have at least one residue.")
  }

  # 3. Name check (skip NA elements)
  if (!purrr::every(x_valid, ~ checkmate::test_named(.x))) {
    cli::cli_abort(c(
      "Must be one or more named integer vectors.",
      "x" = "The input doesn't have names."
    ))
  }

  # 4. Known monosaccharide check (skip NA elements)
  if (!purrr::every(x_valid, ~ all(is_known_composition_component(names(.x))))) {
    cli::cli_abort(c(
      "Must have only known monosaccharides",
      "i" = "Call {.fun available_monosaccharides} to see all known monosaccharides."
    ))
  }

  # 5. Mono type check (skip NA elements)
  mono_types <- .get_comp_mono_types(x_valid)
  if (any(mono_types == "mixed")) {
    cli::cli_abort(c(
      "Must have only one type of monosaccharide.",
      "x" = "Some compositions have mixed monosaccharide types (both generic and concrete)."
    ))
  }
  if (!all(mono_types == mono_types[[1]])) {
    cli::cli_abort(c(
      "Must have only one type of monosaccharide.",
      "x" = "Both generic and concrete compositions exist."
    ))
  }

  # 6. Positive number check (skip NA elements)
  if (!purrr::every(x_valid, ~ all(.x > 0))) {
    cli::cli_abort("Must have only positive numbers of residues.")
  }
}

#' Get monosaccharide types from a list of named integer vectors (composition components)
#' @param x A list of named integer vectors.
#' @returns A character vector of monosaccharide types.
#' @noRd
.get_comp_mono_types <- function(x) {
  x <- purrr::map(x, ~ .x[!names(.x) %in% available_substituents()])  # remove all substituents

  # After removing substituents, each composition must contain at least one monosaccharide.
  # Otherwise, get_mono_type_impl() would be called with an empty vector and fail.
  if (!purrr::every(x, ~ length(.x) > 0)) {
    cli::cli_abort("Each composition must contain at least one monosaccharide (non-substituent residue).")
  }
  purrr::map_chr(x, ~ get_mono_type_impl(names(.x)))
}

# Helper function to check if a name is a known composition component (monosaccharide or substituent)
is_known_composition_component <- function(names) {
  is_known_monosaccharide(names) | names %in% available_substituents()
}

# Helper function to extract substituent types from substituent strings
extract_substituent_types <- function(sub_strings) {
  # Extract substituent types from strings like "3Me", "6S", "4Ac,3Me"
  all_subs <- character(0)

  for (sub_str in sub_strings) {
    if (sub_str == "" || is.na(sub_str)) {
      next
    }

    # Split by commas for multiple substituents
    individual_subs <- stringr::str_split(sub_str, ",")[[1]]
    individual_subs <- individual_subs[individual_subs != ""]

    # Extract substituent names (remove position numbers)
    sub_names <- stringr::str_extract(individual_subs, "[A-Za-z]+$")
    all_subs <- c(all_subs, sub_names)
  }

  all_subs
}

# Helper function to parse a single composition string
parse_single_composition <- function(char) {
  # Reject empty strings
  if (char == "") {
    return(list(composition = NULL, valid = FALSE))
  }

  # Try each parser in sequence until one succeeds
  parsers <- list(.parse_byonic_comp, .parse_simple_comp)
  for (parser in parsers) {
    result <- tryCatch(
      {
        composition <- parser(char)
        list(composition = composition, valid = TRUE)
      },
      error = function(e) NULL
    )
    if (!is.null(result)) {
      return(result)
    }
  }

  # All parsers failed
  list(composition = NULL, valid = FALSE)
}

.parse_byonic_comp <- function(x) {
  # Use regex to find all patterns like "MonoName(number)"
  pattern <- "([A-Za-z0-9]+)\\((\\d+)\\)"
  matches <- stringr::str_extract_all(x, pattern, simplify = FALSE)[[1]]

  if (length(matches) == 0) {
    stop()
  }

  # Check if the entire string was matched (no remaining characters)
  total_matched_length <- sum(stringr::str_length(matches))
  if (total_matched_length != stringr::str_length(x)) {
    stop()
  }

  # Parse each match using stringr::str_match
  parsed_matches <- purrr::map(matches, function(match) {
    match_result <- stringr::str_match(match, pattern)
    list(
      name = match_result[1, 2],  # First capture group
      count = as.integer(match_result[1, 3])  # Second capture group
    )
  })

  # Extract names and counts
  mono_names <- purrr::map_chr(parsed_matches, "name")
  mono_counts <- purrr::map_int(parsed_matches, "count")

  # Create named vector
  comp <- mono_counts
  names(comp) <- mono_names
  comp
}

.parse_simple_comp <- function(x) {
  # "S" and "A" are both NeuAc, "G" is NeuGc
  mono_pattern <- "([HNFSAG])(\\d+)"
  matches <- stringr::str_extract_all(x, mono_pattern, simplify = FALSE)[[1]]

  # Check if no monos were matched
  if (length(matches) == 0) {
    stop()
  }

  # Check if the entire string was matched (no remaining characters)
  total_matched_length <- sum(stringr::str_length(matches))
  if (total_matched_length != stringr::str_length(x)) {
    stop()
  }

  parsed_matches <- purrr::map(matches, function(match) {
    match_result <- stringr::str_match(match, mono_pattern)
    list(
      name = match_result[1, 2],  # First capture group
      count = as.integer(match_result[1, 3])  # Second capture group
    )
  })

  mono_names <- purrr::map_chr(parsed_matches, "name")
  mono_names <- dplyr::recode_values(mono_names,
    "H" ~ "Hex",
    "N" ~ "HexNAc",
    "F" ~ "dHex",
    "S" ~ "NeuAc",
    "A" ~ "NeuAc",
    "G" ~ "NeuGc"
  )
  mono_counts <- purrr::map_int(parsed_matches, "count")

  # Create named vector
  comp <- mono_counts
  names(comp) <- mono_names
  comp
}

#' Reorder composition components
#'
#' Composition components include monosaccharides and substituents.
#' This function reorders the components based on the names.
#' It assumes that all components are valid.
#'
#' The order is determined by the order of the names in [available_monosaccharides()]
#' and [available_substituents()].
#' Monosaccharides are placed before substituents.
#'
#' @param components A named integer vector of composition components.
#' @returns A named integer vector of reordered composition components.
#' @noRd
.reorder_composition_components <- function(components) {
  mono_orders <- available_monosaccharides()
  sub_orders <- available_substituents()
  orders <- c(mono_orders, sub_orders)
  components[order(match(names(components), orders))]
}

# Helper function to format a subset of compositions
format_glycan_composition_subset <- function(x, indices, colored = TRUE) {
  if (!colored) {
    return(format(x)[indices])
  }

  # Format with colors if concrete type
  format_one_colored <- function(comp) {
    if (.is_na_composition_elem(comp)) {
      "<NA>"
    } else {
      mono_names <- names(comp)
      # Add colors to monosaccharides (generic monos automatically get black color)
      mono_names <- add_colors(mono_names, colored = TRUE)
      paste0(mono_names, "(", comp, ")", collapse = "")
    }
  }

  data <- vctrs::vec_data(x)
  comp_data <- vctrs::field(data, "data")[indices]

  purrr::map_chr(comp_data, format_one_colored)
}

#' Check if a composition element is NA
#'
#' NA elements are stored as NULL in the list_of data field.
#' @param x A single composition element (named integer vector or NULL)
#' @returns TRUE if NA (NULL), FALSE otherwise
#' @noRd
.is_na_composition_elem <- function(x) {
  is.null(x)
}

#' Validate a single composition element
#'
#' @param x A single composition element
#' @returns The validated and converted integer vector
#' @noRd
.valid_composition_element <- function(x) {
  if (!checkmate::test_integerish(x)) {
    cli::cli_abort("Each composition must be a named integer vector.")
  }
  if (length(x) == 0) {
    cli::cli_abort("Each composition must have at least one residue.")
  }
  if (!checkmate::test_named(x)) {
    cli::cli_abort("Each composition must be a named integer vector.")
  }
  if (!all(is_known_composition_component(names(x)))) {
    cli::cli_abort(c(
      "Must have only known monosaccharides",
      "i" = "Call {.fun available_monosaccharides} to see all known monosaccharides."
    ))
  }
  if (any(is.na(x)) || !all(x > 0)) {
    cli::cli_abort("Must have only positive numbers of residues.")
  }
  result <- as.integer(x)
  names(result) <- names(x)
  result
}
