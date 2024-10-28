#' Convert the Type of Monosaacharides in a Glycan Graph
#'
#' @description
#' This function converts the all monosaccharides in a glycan graph
#' to a different type.
#' The types are: concrete, generic, and simple.
#' The conversion can only be done from "concrete" to "generic" or "simple",
#' and from "generic" to "simple".
#' Conversion in other orders is not allowed.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param glycan A glycan graph.
#' @param to A character string specifying the target monosaccharide type.
#' It can be "concrete", "generic", or "simple".
#'
#' @return A glycan graph with monosaccharides converted to the target type.
#'
#' @examples
#' concrete_glycan <- n_glycan_core(mono_type = "concrete")
#' convert_glycan_mono_type(concrete_glycan, to = "generic")
#' convert_glycan_mono_type(concrete_glycan, to = "simple")
#' generic_glycan <- n_glycan_core(mono_type = "generic")
#' convert_glycan_mono_type(generic_glycan, to = "simple")
#'
#' @seealso [convert_mono_type()], [decide_glycan_mono_type()], [decide_mono_type()],
#' [ensure_glycan_mono_type()]
#'
#' @export
convert_glycan_mono_type <- function(glycan, to) {
  check_to_arg(to)
  from <- decide_glycan_mono_type(glycan)
  valid_from_to_for_convert_glycan_mono_type(from, to)
  if (inherits(glycan, "ne_glycan_graph")) {
    convert_glycan_mono_type_ne(glycan, from, to)
  } else {  # "dn_glycan_graph"
    convert_glycan_mono_type_dn(glycan, from, to)
  }
}


valid_from_to_for_convert_glycan_mono_type <- function(from, to) {
  tryCatch(
    valid_from_to(from, to),
    error_backward_convert = function(e) {
      cli::cli_abort(c(
        "Cannot convert from {.val {from}} to {.val {to}}.",
        "i" = "Can only convert in this order: concrete -> generic -> simple."
      ),
      call = rlang::expr(convert_glycan_mono_type()),
      class = "error_backward_convert")
    },
    error_convert_self = function(e) {
      cli::cli_abort(
        "It is already {.val {to}}.",
        call = rlang::expr(convert_glycan_mono_type()),
        class = "error_convert_self"
      )
    }
  )
}


#' Ensure the Monosaacharides in a Glycan Graph are of a Specific Type
#'
#' This function ensures the glycan has the given type of monosaccharides.
#' It differs from `convert_glycan_mono_type()` in that it does not raise an error
#' if you try to "ensure the same mono type", e.g. from "generic" to "generic".
#' This is more handy in some cases.
#'
#' @param glycan A glycan graph.
#' @param to A character string specifying the target monosaccharide type.
#'  It can be "concrete", "generic", or "simple".
#'
#' @return A glycan graph with monosaccharides converted to the target type.
#'
#' @examples
#' concrete_glycan <- n_glycan_core(mono_type = "concrete")
#' ensure_glycan_mono_type(concrete_glycan, to = "generic")
#' ensure_glycan_mono_type(concrete_glycan, to = "concrete")
#'
#' @seealso [convert_glycan_mono_type()], [decide_glycan_mono_type()], [decide_mono_type()],
#' [convert_mono_type()]
#'
#' @export
ensure_glycan_mono_type <- function(glycan, to) {
  tryCatch(
    convert_glycan_mono_type(glycan, to),
    error_convert_self = function(e) glycan
  )
}


#' Convert a Monosaacharide to a Different Type
#'
#' @description
#' This function converts a monosaccharide to a different type.
#' The types are: concrete, generic, and simple.
#' The conversion can only be done from "concrete" to "generic" or "simple",
#' and from "generic" to "simple".
#' Conversion in other orders is not allowed.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param mono A character string specifying monosaccharide names.
#' @param to A character string specifying the target monosaccharide type.
#'  It can be "concrete", "generic", or "simple".
#'
#' @return A character string specifying the monosaccharide name in the target type.
#'
#' @examples
#' convert_mono_type(c("Gal", "Hex", "GlcNAc"), to = "simple")
#' convert_mono_type(c("Gal", "Man", "GlcNAc"), to = "generic")
#'
#' @seealso [convert_glycan_mono_type()], [decide_glycan_mono_type()], [decide_mono_type()],
#' [ensure_glycan_mono_type()]
#'
#' @export
convert_mono_type <- function(mono, to) {
  check_to_arg(to)
  from <- decide_mono_type(mono)
  valid_from_to_for_convert_mono_type(mono, from, to)
  convert_mono_type_(mono, from, to)
}


valid_from_to_for_convert_mono_type <- function(mono, from, to) {
  tryCatch(
    valid_from_to(from, to),
    error_backward_convert = function(e) {
      cli::cli_abort(c(
        "These monosaccharides cannot be converted to {.val {to}}: {.val {mono[e$bad_index]}}",
        "i" = "Conversion could only be done in this direction: concrete -> generic -> simple"
      ), call = rlang::expr(convert_mono_type()))
    },
    error_convert_self = function(e) {
      cli::cli_abort(
        "These monosaccharides are already {.val {to}}: {.val {mono[e$bad_index]}}",
        call = rlang::expr(convert_mono_type())
      )
    }
  )
}


#' Decide the Type of Monosaacharides in a Glycan Graph
#'
#' This function trys to decide the type of monosaccharides in a glycan graph.
#'
#' @details
#' By saying "trys" in the description, it means that the function only
#' checks the first monosaccharide in the graph to decide the type.
#' This is reasonable because the monosaccharides in a glycan graph are always
#' of the same type if it was created with the functions in this package.
#' The validation functions will ensure this.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param glycan A glycan graph.
#'
#' @return A character string specifying the monosaccharide type.
#'
#' @examples
#' decide_glycan_mono_type(n_glycan_core(mono_type = "concrete"))
#' decide_glycan_mono_type(n_glycan_core(mono_type = "generic"))
#' decide_glycan_mono_type(n_glycan_core(mono_type = "simple"))
#'
#' @seealso [convert_glycan_mono_type()], [convert_mono_type()], [decide_mono_type()],
#' [ensure_glycan_mono_type()]
#'
#' @export
decide_glycan_mono_type <- function(glycan) {
  stopifnot(is_glycan(glycan))
  if (inherits(glycan, "ne_glycan_graph")) {
    decide_glycan_mono_type_ne(glycan)
  } else {  # "dn_glycan_graph"
    decide_glycan_mono_type_dn(glycan)
  }
}


#' Decide the Type of a Monosaacharide
#'
#' This function decides the type of a monosaccharide name.
#'
#' @details
#' # Three types of monosaccharides
#' There are three types of monosaccharides:
#' - concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.
#' - generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.
#' - simple: e.g. "H", "N", "S", "F".
#'
#' For the full list of monosaccharides, see `glyrepr::monosaccharides`.
#'
#' @param monos A character string specifying monosaccharide names.
#'
#' @return A character string specifying the monosaccharide types.
#'
#' @examples
#' decide_mono_type(c("Gal", "Hex", "H"))
#'
#' @seealso [convert_glycan_mono_type()], [convert_mono_type()], [decide_glycan_mono_type()],
#' [ensure_glycan_mono_type()]
#'
#' @export
decide_mono_type <- function(monos) {
  if (!is.character(monos)) {
    rlang::abort("Mono must be a character string.")
  }
  result <- vector("character", length = length(monos))
  result[monos %in% monosaccharides$concrete] <- "concrete"
  result[monos %in% monosaccharides$generic] <- "generic"
  result[monos %in% monosaccharides$simple] <- "simple"
  unknown <- monos[result == ""]
  if (length(unknown) > 0) {
    cli::cli_abort("Unknown monosaccharide: {.val {unknown}}.")
  }
  result
}


check_to_arg <- function(to) {
  if (!length(to) == 1) {
    rlang::abort("Only one `to` mono type can be specified.")
  }
  if (!(to %in% c("concrete", "generic", "simple"))) {
    rlang::abort("Must be one of: concrete, generic, simple.")
  }
}


valid_from_to <- function(from, to) {
  factorize_types <- function(types) {
    factor(types, levels = c("simple", "generic", "concrete"), ordered = TRUE)
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
  if (any(from == to)) {
    rlang::abort(
      class = "error_convert_self",
      bad_index = which(from == to)
    )
  }
}


decide_glycan_mono_type_ne <- function(glycan) {
  first_mono <- igraph::vertex_attr(glycan, "mono")[[1]]
  decide_mono_type(first_mono)
}


decide_glycan_mono_type_dn <- function(glycan) {
  mono_nodes <- which(igraph::V(glycan)$type == "mono")
  mono_names <- igraph::vertex_attr(glycan, "mono")[mono_nodes]
  first_mono <- mono_names[[1]]
  decide_mono_type(first_mono)
}


convert_mono_type_ <- function(mono, from, to) {
  convert_one_mono_type <- function(mono, from, to) {
    from_ <- monosaccharides[[from]]
    to_ <- monosaccharides[[to]]
    to_[match(mono, from_)]
  }
  purrr::map2_chr(mono, from, convert_one_mono_type, to = to)
}


convert_glycan_mono_type_ne <- function(glycan, from, to) {
  new_names <- convert_mono_type_(igraph::V(glycan)$mono, from, to)
  igraph::set_vertex_attr(glycan, "mono", value = new_names)
}


convert_glycan_mono_type_dn <- function(glycan, from, to) {
  mono_nodes <- which(igraph::V(glycan)$type == "mono")
  old_names <- igraph::vertex_attr(glycan, "mono")[mono_nodes]
  new_names <- convert_mono_type_(old_names, from, to)
  igraph::set_vertex_attr(glycan, "mono", value = new_names, index = mono_nodes)
}
