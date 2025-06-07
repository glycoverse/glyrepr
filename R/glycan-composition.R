#' Get Composition of a Glycan
#'
#' Get the composition of a glycan graph.
#' The composition is a named integer vector.
#'
#' @param glycan A glycan graph.
#'
#' @return A named integer vector.
#'
#' @examples
#' glycan <- n_glycan_core()
#' get_composition(glycan)
#' 
#' @export
get_composition <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_graph")
  monos <- igraph::V(glycan)$mono
  result_tb <- table(monos)
  result <- as.integer(result_tb)
  names(result) <- names(result_tb)
  result
}

new_composition <- function(x) {
  # Use vctrs::list_of instead of new_list_of
  x <- vctrs::list_of(!!!x, .ptype = integer())
  first_monos <- purrr::map_chr(x, ~names(.x)[1])
  # Handle unknown monosaccharides gracefully by assigning "unknown" type
  mono_types <- purrr::map_chr(first_monos, function(mono) {
    if (is_known_monosaccharide(mono)) {
      decide_mono_type(mono)
    } else {
      "unknown"
    }
  })
  vctrs::new_rcrd(list(data = x, mono_type = mono_types), class = "glyrepr_composition")
}

valid_composition <- function(x) {
  valid_one <- function(x) {
    # Check if the composition is named
    if (is.null(names(x))) {
      cli::cli_abort("{.arg ...} must be named.")
    }
    # Check if the composition has only known monosaccharides
    if (!all(is_known_monosaccharide(names(x)))) {
      cli::cli_abort(c(
        "{.arg ...} must have only known monosaccharides.",
        "i" = "Call {.fun known_monosaccharides} to see all known monosaccharides."
      ))
    }
    # Check if all residues have the same type (simple, generic, or concrete)
    local({
      mono_types <- decide_mono_type(names(x))
      if (length(unique(mono_types)) > 1) {
        cli::cli_abort(c(
          "{.arg ...} must have only one type of monosaccharide.",
          "i" = "Call {.fun decide_mono_type} to see the type of each monosaccharide."
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

#' @export
composition <- function(...) {
  args <- list(...)
  x <- purrr::map(args, ~ {
    result <- as.integer(.x)
    names(result) <- names(.x)
    # Sort by monosaccharides tibble order (top to bottom)
    mono_order <- get_monosaccharide_order(names(result))
    result <- result[order(mono_order)]
    result
  })
  x <- new_composition(x)
  x <- valid_composition(x)
  x
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
  mono_type <- decide_mono_type(mono_names[1])
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
    if (mono_type == "simple") {
      paste0(names(comp), comp, collapse = "")
    } else {
      paste0(names(comp), "(", comp, ")", collapse = "")
    }
  }
  data <- vctrs::vec_data(x)
  purrr::map2_chr(vctrs::field(data, "data"), vctrs::field(data, "mono_type"), format_one)
}
