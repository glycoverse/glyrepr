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
#' strucs <- c(n_glycan_core(), o_glycan_core_1())
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

  convert_mono_type_impl(x)
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

  if (length(x) == 0) {
    return(x)
  }

  # Get current mono types
  current_type <- get_mono_type.glyrepr_composition(x)
  if (current_type == "generic") {
    return(x)
  }

  # Convert each composition
  compositions <- vctrs::field(x, "data")
  new_compositions <- purrr::map(compositions, function(comp) {
    # Convert monosaccharide names
    old_names <- names(comp)
    new_names <- convert_mono_type_impl(old_names)

    # Create new composition with converted names, aggregating counts for duplicates
    result <- tapply(comp, new_names, sum, na.rm = TRUE)
    res_names <- names(result)
    result <- as.integer(result)
    names(result) <- res_names

    .reorder_composition_components(result)
  })

  # Create new composition vector
  new_glycan_composition(new_compositions)
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
#' # Special monosaccharides
#'
#' Some monosaccharides are special in that they have no generic names in database or literature.
#' For example, "Mur" is a rare monosaccharide that has no popular generic name.
#' In `glyrepr`, we assign a "g" prefix to these monosaccharides as their generic names.
#' This includes "gNeu", "gKdn", "gPse", "gLeg", "gAci", "g4eLeg", "gBac", "gKdo", "gMur".
#' These names might only be meaningful inside `glycoverse`.
#' Take care when you export results from `glycoverse` functions to other analysis tools.
#'
#' @param x Either of these objects:
#'   - A character vector of monosaccharide names;
#'   - A glycan composition vector ("glyrepr_composition" object);
#'   - A glycan structure vector ("glyrepr_structure" object).
#'
#' @returns
#'   - For character input, returns a character vector of the same length as `x`.
#'   - For `glyrepr_structure` and `glyrepr_composition` input, returns a character scalar.
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
  is_concrete <- x %in% monosaccharides$concrete
  is_generic <- x %in% monosaccharides$generic
  is_unknown <- !(is_concrete | is_generic)
  if (any(is_unknown)) {
    cli::cli_abort("Unknown monosaccharide: {.val {x[is_unknown]}}.")
  }
  result[is_concrete] <- "concrete"
  result[is_generic] <- "generic"
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

  if (length(x) == 0) {
    return(character())
  }

  graphs1 <- attr(x, "graphs")[[1]]
  monos1 <- igraph::V(graphs1)$mono
  get_mono_type_impl(monos1)
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

  if (length(x) == 0) {
    return(character())
  }

  # All compositions in a composition vector must be of the same type.
  # Therefore, we just need to check the first composition.
  # Filter out substituents before determining monosaccharide type.
  comp_names <- names(vctrs::field(x, "data")[[1]])
  mono_names <- comp_names[!comp_names %in% available_substituents()]
  get_mono_type_impl(mono_names)
}

#' Decide mono type from a vector of monosaccharide names
#'
#' This function handles the special cases where some monosaccharides
#' have the same name for both generic and concrete types.
#' For example, "Mur" is both a generic and concrete monosaccharide.
#'
#' This function is used internally when creating [glycan_composition()] and [glycan_structure()].
#'
#' @param x A character vector of monosaccharide names.
#' @returns A character scalar of monosaccharide type. Can be "concrete", "generic", "mixed", or "unknown".
#' @noRd
get_mono_type_impl <- function(x) {
  types <- tryCatch(
    get_mono_type.character(x),
    error = function(e) "unknown"
  )

  unique_types <- unique(types)
  if (length(unique_types) > 1) {
    "mixed"
  } else {
    unique_types[[1]]
  }
}

get_graph_mono_type <- function(graph) {
  # This function is the implementation of get_mono_type for a single glycan graph.
  monos <- igraph::vertex_attr(graph, "mono")
  get_mono_type_impl(monos)
}

#' Convert monosaccharide names to generic type
#'
#' This function converts concrete monosaccharide names to generic type.
#' Generic monosaccharides are returned as is.
#' It assumes that all monosaccharides are of the same type.
#' Substituents are not converted.
#'
#' @param monos A character vector of monosaccharide names.
#' @returns A character vector of monosaccharide names in generic type.
#' @noRd
convert_mono_type_impl <- function(monos) {
  mono_type <- get_mono_type_impl(monos)
  if (mono_type == "generic") {
    return(monos)
  }
  if (mono_type == "mixed") {
    cli::cli_abort("Mixed monosaccharide types are not supported.")
  }
  from_ <- monosaccharides[["concrete"]]
  to_ <- monosaccharides[["generic"]]
  dplyr::if_else(
    monos %in% available_substituents(),
    monos,
    to_[match(monos, from_)]
  )
}

convert_glycan_mono_type_impl <- function(glycan, from, to) {
  old_names <- igraph::V(glycan)$mono
  new_names <- convert_mono_type_impl(old_names)
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
