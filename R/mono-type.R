#' Convert Monosaccharides to Generic Type
#'
#' This function converts monosaccharide types of monosaccharide characters,
#' glycan compositions, or glycan structures from concrete to generic type.
#' This is a simplified version that only supports conversion from "concrete"
#' to "generic" monosaccharides.
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
#'
#' @returns A new object of the same class as `x`
#' with monosaccharides converted to generic type.
#'
#' @examples
#' # Convert character vectors
#' convert_to_generic(c("Gal", "GlcNAc"))
#'
#' # Convert glycan compositions
#' comps <- glycan_composition(
#'   c(Gal = 5, GlcNAc = 2),
#'   c(Glc = 5, GalNAc = 4, Fuc = 1)
#' )
#' convert_to_generic(comps)
#'
#' # Convert glycan structures
#' strucs <- glycan_structure(
#'   n_glycan_core(),
#'   o_glycan_core_1()
#' )
#' convert_to_generic(strucs)
#'
#' @export
convert_to_generic <- function(x) {
  UseMethod("convert_to_generic")
}

#' @export
#' @rdname convert_to_generic
convert_to_generic.character <- function(x) {
  checkmate::assert_character(x)

  from <- get_mono_type(x)

  # Check if conversion is valid (no backward conversion)
  if (any(from == "generic")) {
    # Already generic, return as-is
    return(x)
  }

  res <- convert_mono_type_(x, from, "generic")

  bad_monos <- x[is.na(res)]
  if (length(bad_monos) > 0) {
    cli::cli_warn(
      "Some monosaccharides cannot be converted to generic: {.val {bad_monos}}.",
      call = rlang::expr(convert_to_generic())
    )
  }

  res
}

#' @export
#' @rdname convert_to_generic
convert_to_generic.glyrepr_structure <- function(x) {
  if (!is_glycan_structure(x)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }

  # Get current mono types
  from <- get_mono_type(x)

  # Use spmap_structure with optimized implementation
  spmap_structure(list(x, from), function(graph, from) {
    if (from == "generic") {
      return(graph)
    }
    convert_glycan_mono_type_impl(graph, from, "generic")
  })
}

#' @export
#' @rdname convert_to_generic
convert_to_generic.glyrepr_composition <- function(x) {
  if (!is_glycan_composition(x)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_composition vector.",
      "i" = "Use `glycan_composition()` to create a glyrepr_composition from named vectors."
    ))
  }

  # Get current mono types
  data <- vctrs::vec_data(x)
  current_types <- vctrs::field(data, "mono_type")

  # Check if all are already generic
  if (all(current_types == "generic")) {
    return(x)
  }

  # Convert each composition
  compositions <- vctrs::field(data, "data")
  new_compositions <- purrr::map2(compositions, current_types, function(comp, from_type) {
    # If already generic, return as-is
    if (from_type == "generic") {
      return(comp)
    }

    # Convert monosaccharide names
    old_names <- names(comp)
    new_names <- convert_mono_type_(old_names, from_type, "generic")

    # Check for unconvertible monosaccharides
    bad_names <- old_names[is.na(new_names)]
    if (length(bad_names) > 0) {
      cli::cli_abort(
        "Some monosaccharides in composition cannot be converted to generic: {.val {bad_names}}.",
        call = rlang::expr(convert_to_generic())
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
#' @seealso [convert_to_generic()]
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
    cli::cli_abort(c(
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
    cli::cli_abort(c(
      "Input must be a glyrepr_composition vector.",
      "i" = "Use `glycan_composition()` to create a glyrepr_composition from named vectors."
    ))
  }
  
  data <- vctrs::vec_data(x)
  vctrs::field(data, "mono_type")
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
      call = rlang::expr(convert_to_generic())
    )
  }
}
