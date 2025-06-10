#' Create a Glycan Composition
#'
#' Create a glycan composition from a list of named integer vectors.
#'
#' @param ... Named integer vectors.
#'   Names are monosaccharides, values are numbers of residues.
#' @param x A list of named integer vectors.
#'
#' @return A glyrepr_composition object.
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
#'
#' @export
glycan_composition <- function(...) {
  args <- list(...)
  x <- purrr::map(args, ~ {
    result <- as.integer(.x)
    names(result) <- names(.x)
    # Sort by monosaccharides tibble order (top to bottom)
    mono_order <- get_monosaccharide_order(names(result))
    result <- result[order(mono_order)]
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

valid_glycan_composition <- function(x) {
  valid_one <- function(x) {
    # Allow empty compositions (for conversion results)
    if (length(x) == 0) {
      return(TRUE)
    }
    
    # Check if the composition is named
    if (is.null(names(x))) {
      cli::cli_abort("{.arg ...} must be named.")
    }
    # Check if the composition has only known monosaccharides
    if (!all(is_known_monosaccharide(names(x)))) {
      cli::cli_abort(c(
        "{.arg ...} must have only known monosaccharides.",
        "i" = "Call {.fun available_monosaccharides} to see all known monosaccharides."
      ))
    }
    # Check if all residues have the same type (generic or concrete)
    local({
      mono_types <- get_mono_type(names(x))
      if (length(unique(mono_types)) > 1) {
        cli::cli_abort(c(
          "{.arg ...} must have only one type of monosaccharide.",
          "i" = "Call {.fun get_mono_type} to see the type of each monosaccharide."
        ))
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
#' Convert an object to a glycan composition.
#'
#' @param x An object to convert to a glycan composition.
#'   Can be a named integer vector, a list of named integer vectors,
#'   or an existing glyrepr_composition object.
#'
#' @return A glyrepr_composition object.
#'
#' @examples
#' # Convert a named vector
#' as_glycan_composition(c(Hex = 5, HexNAc = 2))
#' # Convert a list of named vectors
#' as_glycan_composition(list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1)))
#' # Convert an existing composition (returns as-is)
#' comp <- glycan_composition(c(Hex = 5, HexNAc = 2))
#' as_glycan_composition(comp)
#'
#' @export
as_glycan_composition <- function(x) {
  UseMethod("as_glycan_composition")
}

#' @export
as_glycan_composition.glyrepr_composition <- function(x) {
  x
}

#' @export
as_glycan_composition.glyrepr_structure <- function(x) {
  # Use structure_map to convert each structure to composition
  compositions <- structure_map(x, function(graph) {
    monos <- igraph::V(graph)$mono
    result_tb <- table(monos)
    result <- as.integer(result_tb)
    names(result) <- names(result_tb)
    # Sort by monosaccharides tibble order (top to bottom)
    mono_order <- get_monosaccharide_order(names(result))
    result <- result[order(mono_order)]
    result
  })
  
  # Create composition vector
  do.call(glycan_composition, compositions)
}

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

#' @export
obj_print_data.glyrepr_composition <- function(x, ..., max_n = 10, colored = TRUE) {
  if (length(x) == 0) {
    return()
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
  formatted <- purrr::map2_chr(vctrs::field(data, "data"), vctrs::field(data, "mono_type"), format_one_colored)
  
  n <- length(formatted)
  n_show <- min(n, max_n)
  
  # Print each composition on its own line with indexing, up to max_n
  for (i in seq_len(n_show)) {
    cat("[", i, "] ", formatted[i], "\n", sep = "")
  }
  if (n > max_n) {
    cat("... (", n - max_n, " more not shown)\n", sep = "")
  }
}
