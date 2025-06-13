#' Convert Monosaccharide Types
#'
#' This function converts monosaccharide types of monosaccharide characters,
#' glycan compositions, or glycan structures.
#' Supported types: "concrete" and "generic" (see details below).
#' The conversion can only be done from "concrete" to "generic".
#'
#' @details
#'
#' # Two types of monosaccharides
#'
#' There are two types of monosaccharides:
#' - concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.
#' - generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.
#'
#' For the full list of monosaccharides, use [available_monosaccharides()].
#' 
#' @param x Either of these objects:
#'   - A character of monosaccharide;
#'   - A glycan composition vector ("glyrepr_composition" object);
#'   - A glycan structure vector ("glyrepr_structure" object).
#' @param to A character string specifying the target monosaccharide type.
#'
#' @returns A new object of the same class as `x` 
#' with monosaccharides converted to the target type.
#'
#' @examples
#' # Allowed
#' convert_mono_type(c("Gal", "GlcNAc"), to = "generic")  # concrete -> generic
#' 
#' # Not allowed
#' \dontrun{
#' convert_mono_type(c("Hex", "HexNAc"), to = "concrete")  # generic -> concrete
#' }
#'
#' # Convert glycan compositions
#' comps <- glycan_composition(
#'   c(Hex = 5, HexNAc = 2),
#'   c(Hex = 5, HexNAc = 4, dHex = 1)
#' )
#' convert_mono_type(comps, to = "generic")
#' 
#' # Convert glycan structures
#' strucs <- glycan_structure(
#'   n_glycan_core(),
#'   o_glycan_core_1()
#' )
#' convert_mono_type(strucs, to = "generic")
#'
#' @export
convert_mono_type <- function(x, to) {
  UseMethod("convert_mono_type")
}

#' @export
#' @rdname convert_mono_type
convert_mono_type.character <- function(x, to) {
  checkmate::assert_character(x)
  checkmate::assert_choice(to, c("concrete", "generic"))

  from <- get_mono_type(x)
  
  # Check if conversion is valid (no backward conversion)
  tryCatch(
    valid_from_to(from, to, strict = FALSE),
    error_backward_convert = function(e) {
      cli::cli_abort(c(
        "These monosaccharides cannot be converted to {.val {to}}: {.val {x[e$bad_index]}}",
        "i" = "Conversion could only be done in this direction: concrete -> generic"
      ), call = rlang::expr(convert_mono_type()))
    }
  )
  
  # If already the target type, return as-is
  if (all(from == to)) {
    return(x)
  }
  
  res <- convert_mono_type_(x, from, to)

  bad_monos <- x[is.na(res)]
  if (length(bad_monos) > 0) {
    cli::cli_warn(
      "Some monosaccharides cannot be converted to {.val {to}}: {.val {bad_monos}}.",
      call = rlang::expr(convert_mono_type())
    )
  }

  res
}

#' @export
#' @rdname convert_mono_type
convert_mono_type.glyrepr_structure <- function(x, to) {
  if (!is_glycan_structure(x)) {
    rlang::abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }
  
  checkmate::assert_choice(to, c("concrete", "generic"))
  
  # Performance optimization: Batch type checking for the entire vector
  types <- get_mono_type(x)
  
  # Early exit if all structures are already the target type
  if (all(types == to)) {
    return(x)
  }
  
  # Validate conversion direction only once for unique types
  unique_types <- unique(types)
  for (from_type in unique_types) {
    if (from_type != to) {
      tryCatch(
        valid_from_to(from_type, to, strict = FALSE),
        error_backward_convert = function(e) {
          cli::cli_abort(c(
            "Cannot convert from {.val {from_type}} to {.val {to}}.",
            "i" = "Can only convert in this order: concrete -> generic."
          ),
          call = rlang::expr(convert_mono_type()),
          class = "error_backward_convert")
        }
      )
    }
  }
  
  # Use smap_structure with optimized lambda function
  # Skip redundant type checking and validation since we already did it
  smap_structure(x, function(graph) {
    from <- get_graph_mono_type(graph)
    if (from == to) {
      return(graph)
    }
    convert_glycan_mono_type_impl(graph, from, to)
  })
}

#' @export
#' @rdname convert_mono_type
convert_mono_type.glyrepr_composition <- function(x, to) {
  if (!is_glycan_composition(x)) {
    rlang::abort(c(
      "Input must be a glyrepr_composition vector.",
      "i" = "Use `glycan_composition()` to create a glyrepr_composition from named vectors."
    ))
  }
  
  checkmate::assert_choice(to, c("concrete", "generic"))
  
  # Get current mono types
  data <- vctrs::vec_data(x)
  current_types <- vctrs::field(data, "mono_type")
  
  # Check if all are already the target type
  if (all(current_types == to)) {
    return(x)
  }
  
  # Convert each composition
  compositions <- vctrs::field(data, "data")
  new_compositions <- purrr::map2(compositions, current_types, function(comp, from_type) {
    # Check if conversion is valid
    tryCatch(
      valid_from_to(from_type, to, strict = FALSE),
      error_backward_convert = function(e) {
        cli::cli_abort(c(
          "Cannot convert composition from {.val {from_type}} to {.val {to}}.",
          "i" = "Conversion could only be done in this direction: concrete -> generic"
        ), call = rlang::expr(convert_mono_type()))
      }
    )
    
    # If already target type, return as-is
    if (from_type == to) {
      return(comp)
    }
    
    # Convert monosaccharide names
    old_names <- names(comp)
    new_names <- convert_mono_type_(old_names, from_type, to)
    
    # Check for unconvertible monosaccharides
    bad_names <- old_names[is.na(new_names)]
    if (length(bad_names) > 0) {
      cli::cli_abort(
        "Some monosaccharides in composition cannot be converted to {.val {to}}: {.val {bad_names}}.",
        call = rlang::expr(convert_mono_type())
      )
    }
    
    # Create new composition with converted names, aggregating counts for duplicates
    if (length(comp) == 0) {
      # Handle empty composition
      result <- integer(0)
      names(result) <- character(0)
    } else {
      result <- tapply(comp, new_names, sum, na.rm = TRUE)
      result_names <- names(result)  # Save names before conversion
      result <- as.integer(result)
      names(result) <- result_names  # Restore names
      
      # Sort by monosaccharide order
      mono_order <- get_monosaccharide_order(names(result))
      result <- result[order(mono_order)]
    }
    
    result
  })
  
  # Create new composition vector
  do.call(glycan_composition, new_compositions)
}

#' Get Monosaccharide Types
#'
#' This function determines the type of monosaccharides in character vectors,
#' glycan compositions, or glycan structures.
#' Supported types: "concrete" and "generic" (see details below).
#'
#' @details
#'
#' # Two types of monosaccharides
#'
#' There are two types of monosaccharides:
#' - concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.
#' - generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.
#'
#' For the full list of monosaccharides, use [available_monosaccharides()].
#'
#' @param x Either of these objects:
#'   - A character vector of monosaccharide names;
#'   - A glycan composition vector ("glyrepr_composition" object);
#'   - A glycan structure vector ("glyrepr_structure" object).
#'
#' @return A character vector specifying the monosaccharide type(s).
#'   For structures and compositions, returns the type for each element.
#'
#' @examples
#' # Character vector
#' get_mono_type(c("Gal", "Hex"))
#' 
#' # Glycan structures
#' get_mono_type(n_glycan_core(mono_type = "concrete"))
#' get_mono_type(n_glycan_core(mono_type = "generic"))
#' 
#' # Glycan compositions
#' comp <- glycan_composition(c(Glc = 2, GalNAc = 1))
#' get_mono_type(comp)
#'
#' @seealso [convert_mono_type()]
#'
#' @export
get_mono_type <- function(x) {
  UseMethod("get_mono_type")
}

#' @export
#' @rdname get_mono_type
get_mono_type.character <- function(x) {
  checkmate::assert_character(x)
  result <- vector("character", length = length(x))
  result[x %in% monosaccharides$concrete] <- "concrete"
  result[x %in% monosaccharides$generic] <- "generic"
  unknown <- x[result == ""]
  if (length(unknown) > 0) {
    cli::cli_abort("Unknown monosaccharide: {.val {unknown}}.")
  }
  result
}

#' @export
#' @rdname get_mono_type
get_mono_type.glyrepr_structure <- function(x) {
  if (!is_glycan_structure(x)) {
    rlang::abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }
  
  data <- vctrs::vec_data(x)
  vctrs::field(data, "mono_type")
}

#' @export
#' @rdname get_mono_type
get_mono_type.glyrepr_composition <- function(x) {
  if (!is_glycan_composition(x)) {
    rlang::abort(c(
      "Input must be a glyrepr_composition vector.",
      "i" = "Use `glycan_composition()` to create a glyrepr_composition from named vectors."
    ))
  }
  
  data <- vctrs::vec_data(x)
  vctrs::field(data, "mono_type")
}

valid_from_to <- function(from, to, strict) {
  factorize_types <- function(types) {
    factor(types, levels = c("generic", "concrete"), ordered = TRUE)
  }
  from <- factorize_types(from)
  to <- factorize_types(to)

  bad_order <- from < to
  if (any(bad_order)) {
    rlang::abort(
      class = "error_backward_convert",
      bad_index = which(bad_order)
    )
  }
  if (strict && any(from == to)) {
    rlang::abort(
      class = "error_convert_self",
      bad_index = which(from == to)
    )
  }
}

get_graph_mono_type <- function(graph) {
  # This function is the implementation of get_mono_type for a single glycan graph.
  first_mono <- igraph::vertex_attr(graph, "mono")[[1]]
  get_mono_type(first_mono)
}

convert_mono_type_ <- function(mono, from, to) {
  purrr::map2_chr(mono, from, convert_one_mono_type, to = to)
}

convert_one_mono_type <- function(mono, from, to) {
  from_ <- monosaccharides[[from]]
  to_ <- monosaccharides[[to]]
  to_[match(mono, from_)]  # it might be NA
}

convert_glycan_mono_type_impl <- function(glycan, from, to) {
  old_names <- igraph::V(glycan)$mono
  new_names <- convert_mono_type_(old_names, from, to)
  raise_error_for_na(old_names, new_names, to)
  igraph::set_vertex_attr(glycan, "mono", value = new_names)
}

raise_error_for_na <- function(old_names, new_names, to) {
  bad_names <- old_names[is.na(new_names)]
  if (length(bad_names) > 0) {
    cli::cli_abort(
      "Some monosaccharides cannot be converted to {.val {to}}: {.val {bad_names}}.",
      call = rlang::expr(convert_mono_type())
    )
  }
}
