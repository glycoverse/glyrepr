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
#' - Monosaccharides: either generic (e.g., "Hex", "HexNAc") or concrete 
#'   (e.g., "Glc", "Gal"). All monosaccharides in a composition must be 
#'   of the same type.
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
  if (!purrr::every(args, checkmate::test_integerish)) {
    cli::cli_abort(c(
      "Must be one or more named integer vectors.",
      "i" = "You might want to use {.fn as_glycan_composition} for more flexible input."
    ))
  }
  x <- purrr::map(args, ~ {
    result <- as.integer(.x)
    names(result) <- names(.x)
    # Sort by monosaccharides tibble order (top to bottom), with substituents at the end
    comp_order <- get_composition_component_order(names(result))
    result <- result[order(comp_order)]
    result
  })
  x <- new_glycan_composition(x)
  x <- valid_glycan_composition(x)
  x
}

new_glycan_composition <- function(x) {
  # Use vctrs::list_of instead of new_list_of
  x <- vctrs::list_of(!!!x, .ptype = integer())
  first_monos <- purrr::map_chr(x, function(.x) {
    if (length(.x) == 0) {
      return(NA_character_)
    } else {
      return(names(.x)[1])
    }
  })
  # Handle unknown monosaccharides gracefully by assigning "unknown" type
  mono_types <- purrr::map_chr(first_monos, function(mono) {
    if (is.na(mono)) {
      return("unknown")  # For empty compositions
    } else if (is_known_monosaccharide(mono)) {
      get_mono_type(mono)
    } else {
      "unknown"
    }
  })
  vctrs::new_rcrd(list(data = x, mono_type = mono_types), class = "glyrepr_composition")
}

# Helper function to check if a name is a known composition component (monosaccharide or substituent)
is_known_composition_component <- function(names) {
  is_known_monosaccharide(names) | names %in% available_substituents()
}

valid_glycan_composition <- function(x) {
  valid_one <- function(x) {
    # Reject empty compositions
    if (length(x) == 0) {
      cli::cli_abort("Each composition in {.arg ...} must have at least one residue.")
    }
    
    # Check if the composition is named
    if (is.null(names(x))) {
      cli::cli_abort("{.arg ...} must be named.")
    }
    # Check if the composition has only known monosaccharides and substituents
    if (!all(is_known_composition_component(names(x)))) {
      cli::cli_abort(c(
        "{.arg ...} must have only known monosaccharides.",
        "i" = "Call {.fun available_monosaccharides} to see all known monosaccharides."
      ))
    }
    # Check if all monosaccharides have the same type (generic or concrete)
    # Substituents are allowed to mix with either type
    local({
      component_names <- names(x)
      mono_names <- component_names[is_known_monosaccharide(component_names)]
      
      if (length(mono_names) > 0) {
        mono_types <- get_mono_type(mono_names)
        if (length(unique(mono_types)) > 1) {
          cli::cli_abort(c(
            "{.arg ...} must have only one type of monosaccharide.",
            "i" = "Call {.fun get_mono_type} to see the type of each monosaccharide."
          ))
        }
      }
    })
    # Check if all numbers of residues are positive
    if (!all(x > 0)) {
      cli::cli_abort("{.arg ...} must have only positive numbers of residues.")
    }
  }

  data <- vctrs::vec_data(x)
  # The data is stored in data$data, which is a list_of
  purrr::walk(vctrs::field(data, "data"), valid_one)
  x
}

#' Convert to Glycan Composition
#'
#' Convert an object to a glycan composition. The resulting composition can 
#' contain both monosaccharides and substituents.
#'
#' @param x An object to convert to a glycan composition.
#'   Can be a named integer vector, a list of named integer vectors,
#'   a glycan structure vector,
#'   or an existing glyrepr_composition object.
#'
#' @returns A glyrepr_composition object.
#'
#' @details
#' When converting from glycan structures, both monosaccharides and substituents 
#' are counted. Substituents are extracted from the `sub` attribute of each 
#' vertex in the structure. For example, a vertex with `sub = "3Me"` 
#' contributes one "Me" substituent to the composition.
#'
#' @examples
#' # Convert a named vector
#' as_glycan_composition(c(Hex = 5, HexNAc = 2))
#' 
#' # Convert a named vector with substituents
#' as_glycan_composition(c(Glc = 2, Gal = 1, Me = 1, S = 1))
#' 
#' # Convert a list of named vectors
#' as_glycan_composition(list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1)))
#' 
#' # Convert an existing composition (returns as-is)
#' comp <- glycan_composition(c(Hex = 5, HexNAc = 2))
#' as_glycan_composition(comp)
#' 
#' # Convert a glycan structure vector
#' strucs <- c(n_glycan_core(), o_glycan_core_1())
#' as_glycan_composition(strucs)
#' 
#' # Convert structures with substituents
#' # (This will count both monosaccharides and any substituents present)
#'
#' @export
as_glycan_composition <- function(x) {
  UseMethod("as_glycan_composition")
}

#' @rdname as_glycan_composition
#' @export
as_glycan_composition.glyrepr_composition <- function(x) {
  x
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

#' @rdname as_glycan_composition
#' @export
as_glycan_composition.glyrepr_structure <- function(x) {
  # Get mono_types directly from the structure data
  data <- vctrs::vec_data(x)
  structure_mono_types <- vctrs::field(data, "mono_type")
  
  # Use smap2 to convert each structure to composition with known mono_type
  compositions <- smap2(x, structure_mono_types, function(graph, mono_type) {
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
    comp_order <- get_composition_component_order(names(result))
    result <- result[order(comp_order)]
    result
  })
  
  # Create composition vector
  do.call(glycan_composition, compositions)
}

# Helper function to parse a single composition string
parse_single_composition <- function(char) {
  # Reject empty strings
  if (char == "") {
    return(list(composition = NULL, valid = FALSE))
  }
  
  # Try to parse the character string
  tryCatch({
    # Use regex to find all patterns like "MonoName(number)"
    pattern <- "([A-Za-z0-9]+)\\((\\d+)\\)"
    matches <- stringr::str_extract_all(char, pattern, simplify = FALSE)[[1]]

    if (length(matches) == 0) {
      return(list(composition = NULL, valid = FALSE))
    }

    # Check if the entire string was matched (no remaining characters)
    total_matched_length <- sum(stringr::str_length(matches))
    if (total_matched_length != stringr::str_length(char)) {
      return(list(composition = NULL, valid = FALSE))
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

    return(list(composition = comp, valid = TRUE))

  }, error = function(e) {
    return(list(composition = NULL, valid = FALSE))
  })
}

#' @rdname as_glycan_composition
#' @export
as_glycan_composition.character <- function(x) {
  # Handle empty character vector
  if (length(x) == 0) {
    return(glycan_composition())
  }

  # Handling NA
  if (any(is.na(x))) {
    cli::cli_abort("Cannot parse NA as glycan composition.")
  }
  
  # Parse each character string using the helper function
  parse_result <- purrr::map(x, parse_single_composition)
  
  # Extract validity and compositions
  valid_flags <- purrr::map_lgl(parse_result, "valid")
  compositions <- purrr::map(parse_result, "composition")
  
  # Find invalid indices
  invalid_indices <- which(!valid_flags)
  
  # Check for invalid characters
  if (length(invalid_indices) > 0) {
    cli::cli_abort(c(
      "Characters cannot be parsed as glycan compositions at index {invalid_indices}",
      "i" = "Expected format: 'Hex(5)HexNAc(2)' with monosaccharide names followed by counts in parentheses."
    ))
  }
  
  # Filter out NULL entries (from empty strings)
  compositions <- compositions[!purrr::map_lgl(compositions, is.null)]
  
  # Handle case where all strings were empty
  if (length(compositions) == 0) {
    return(glycan_composition())
  }
  
  # Create composition vector
  do.call(glycan_composition, compositions)
}

#' @rdname as_glycan_composition
#' @export
as_glycan_composition.default <- function(x) {
  if (is.null(names(x)) && is.list(x)) {
    # Handle list of named vectors - validate that all elements are named
    if (!all(purrr::map_lgl(x, ~ !is.null(names(.x))))) {
      cli::cli_abort(c(
        "All elements in the list must be named vectors.",
        "i" = "Each vector in the list should have names indicating monosaccharides."
      ))
    }
    do.call(glycan_composition, x)
  } else if (!is.null(names(x)) && is.numeric(x)) {
    # Handle single named vector
    glycan_composition(x)
  } else {
    cli::cli_abort(c(
      "Cannot convert object of class {.cls {class(x)}} to glyrepr_composition.",
      "i" = "Supported types: named integer vector, list of named integer vectors, or existing glyrepr_composition."
    ))
  }
}

#' @export
as.character.glyrepr_composition <- function(x, ...) {
  format(x)
}

# Helper function to get the order of composition components (monosaccharides + substituents)
get_composition_component_order <- function(component_names) {
  if (length(component_names) == 0) {
    return(integer(0))
  }
  
  # Check if all components are known, if not, return original order
  if (!all(is_known_composition_component(component_names))) {
    # Return sequential order for unknown components
    return(seq_along(component_names))
  }
  
  # Separate monosaccharides and substituents
  mono_mask <- is_known_monosaccharide(component_names)
  mono_names <- component_names[mono_mask]
  sub_names <- component_names[!mono_mask]
  
  # Get order for monosaccharides
  if (length(mono_names) > 0) {
    mono_type <- get_mono_type(mono_names[1])
    mono_order <- get_monosaccharide_order_with_type(mono_names, mono_type)
  } else {
    mono_order <- integer(0)
  }
  
  # Get order for substituents
  if (length(sub_names) > 0) {
    available_subs <- available_substituents()
    sub_order <- purrr::map_int(sub_names, function(sub) {
      idx <- which(available_subs == sub)
      if (length(idx) > 0) {
        return(1000 + idx)  # High numbers to put after monosaccharides
      } else {
        return(2000)  # Even higher for unknown substituents
      }
    })
  } else {
    sub_order <- integer(0)
  }
  
  # Combine orders
  result_order <- integer(length(component_names))
  result_order[mono_mask] <- mono_order
  result_order[!mono_mask] <- sub_order
  
  result_order
}

# Helper function to get the order of monosaccharides based on their position in the tibble
get_monosaccharide_order <- function(mono_names) {
  if (length(mono_names) == 0) {
    return(integer(0))
  }
  
  # Check if all monosaccharides are known, if not, return original order
  if (!all(is_known_monosaccharide(mono_names))) {
    # Return sequential order for unknown monosaccharides
    return(seq_along(mono_names))
  }
  
  # Determine mono_type from the first monosaccharide
  mono_type <- get_mono_type(mono_names[1])
  get_monosaccharide_order_with_type(mono_names, mono_type)
}

# Helper function to get the order of monosaccharides when mono_type is already known
get_monosaccharide_order_with_type <- function(mono_names, mono_type) {
  if (length(mono_names) == 0) {
    return(integer(0))
  }
  
  # Check if all monosaccharides are known, if not, return original order
  if (!all(is_known_monosaccharide(mono_names))) {
    # Return sequential order for unknown monosaccharides
    return(seq_along(mono_names))
  }
  
  # Get the corresponding column from monosaccharides tibble
  mono_column <- monosaccharides[[mono_type]]
  # Get the order for each monosaccharide based on its position in the specific column
  mono_order <- purrr::map_int(mono_names, function(mono) {
    idx <- which(mono_column == mono & !is.na(mono_column))
    if (length(idx) > 0) {
      return(idx[1])
    }
    # If not found, return a negative number
    return(-1)
  })
  mono_order
}

#' @export
vec_ptype_full.glyrepr_composition <- function(x, ...) "glycan_composition"

#' @export
vec_ptype_abbr.glyrepr_composition <- function(x, ...) "comp"

#' @export
format.glyrepr_composition <- function(x, ...) {
  format_one <- function(comp, mono_type) {
    paste0(names(comp), "(", comp, ")", collapse = "")
  }
  data <- vctrs::vec_data(x)
  purrr::map2_chr(vctrs::field(data, "data"), vctrs::field(data, "mono_type"), format_one)
}

#' @export
#' @rdname glycan_composition
is_glycan_composition <- function(x) {
  inherits(x, "glyrepr_composition")
}

#' @export
vec_ptype2.glyrepr_composition.glyrepr_composition <- function(x, y, ...) {
  new_glycan_composition(list())
}

#' @export
vec_cast.glyrepr_composition.glyrepr_composition <- function(x, to, ...) {
  x
}

# Helper function to format a subset of compositions
format_glycan_composition_subset <- function(x, indices, colored = TRUE) {
  if (!colored) {
    return(format(x)[indices])
  }
  
  # Format with colors if concrete type
  format_one_colored <- function(comp, mono_type) {
    mono_names <- names(comp)
    # Add colors to monosaccharides (generic monos automatically get black color)
    if (colored) {
      mono_names <- add_colors(mono_names, colored = TRUE)
    }
    paste0(mono_names, "(", comp, ")", collapse = "")
  }
  
  data <- vctrs::vec_data(x)
  comp_data <- vctrs::field(data, "data")[indices]
  mono_types <- vctrs::field(data, "mono_type")[indices]
  
  purrr::map2_chr(comp_data, mono_types, format_one_colored)
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
