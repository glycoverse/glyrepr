#' Convert the Type of Monosaacharides in a Glycan Structure
#'
#' @description
#' This function converts all monosaccharides in a glycan structure
#' to a different type.
#' The types are: "concrete", "generic", and "simple" (see details below).
#' The conversion can only be done from "concrete" to "generic" or "simple",
#' and from "generic" to "simple".
#' Conversion in other orders is not allowed.
#'
#' Converting monosaccharides types from higher to lower levels can lead to
#' information loss. In real cases, you should avoid doing this directly and
#' just use as high level as possible.
#' This function is for other packages in `glycoverse` that know what they are
#' doing for specific reasons.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#' @inheritSection convert_mono_type Conversion to "simple" type
#'
#' @param glycan A glycan structure.
#' @param to A character string specifying the target monosaccharide type.
#' It can be "concrete", "generic", or "simple".
#' @param strict If `TRUE`, the function will raise an error if the monosaccharides
#' are already in the target type. Default is `TRUE`.
#'
#' @return A glycan structure with monosaccharides converted to the target type.
#'
#' @examples
#' concrete_glycan <- n_glycan_core(mono_type = "concrete")
#' convert_glycan_mono_type(concrete_glycan, to = "generic")
#' convert_glycan_mono_type(concrete_glycan, to = "simple")
#' generic_glycan <- n_glycan_core(mono_type = "generic")
#' convert_glycan_mono_type(generic_glycan, to = "simple")
#'
#' @seealso [convert_mono_type()], [decide_glycan_mono_type()], [decide_mono_type()]
#'
#' @export
convert_glycan_mono_type <- function(glycan, to, strict = TRUE) {
  checkmate::assert_class(glycan, "glycan_structure")
  checkmate::assert_choice(to, c("concrete", "generic", "simple"))
  checkmate::assert_flag(strict)

  from <- decide_glycan_mono_type(glycan)
  valid_from_to_for_convert_glycan_mono_type(from, to, strict)
  convert_glycan_mono_type_impl(glycan, from, to)
}


valid_from_to_for_convert_glycan_mono_type <- function(from, to, strict) {
  tryCatch(
    valid_from_to(from, to, strict),
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
#' @section Conversion to "simple" type:
#'
#' When converting "concrete" to "simple", the following mapping is used:
#' - Glc, Man, Gal -> H
#' - GlcNAc, GalNAc -> N
#' - Fuc -> F
#' - Neu5Ac -> A
#' - Neu5Gc -> G
#' - Sia -> S
#'
#' When converting "generic" to "simple", the following mapping is used:
#' - Hex -> H
#' - HexNAc -> N
#' - dHex -> F
#' - NeuAc -> A
#' - NeuGc -> G
#' - Sia -> S
#'
#' Note that conversion from "generic" to "simple" uses an over-simplified
#' rule, but this mapping is sufficient for most cases.
#'
#' @param mono A character string specifying monosaccharide names.
#' @param to A character scalar specifying the target monosaccharide type.
#'  It can be "concrete", "generic", or "simple".
#' @param strict If `TRUE`, the function will raise an error if the monosaccharides
#' are already in the target type. Default is `TRUE`.
#'
#' @return A character string specifying the monosaccharide name in the target type.
#' If a monosaccharide cannot be converted, it will be `NA`.
#' If any NA values are present in the input, a warning will be thrown.
#'
#' @examples
#' convert_mono_type(c("Gal", "Hex", "GlcNAc"), to = "simple")
#' convert_mono_type(c("Gal", "Man", "GlcNAc"), to = "generic")
#'
#' @seealso [convert_glycan_mono_type()], [decide_glycan_mono_type()], [decide_mono_type()]
#'
#' @export
convert_mono_type <- function(mono, to, strict = TRUE) {
  checkmate::assert_character(mono)
  checkmate::assert_choice(to, c("concrete", "generic", "simple"))
  checkmate::assert_flag(strict)

  from <- decide_mono_type(mono)
  valid_from_to_for_convert_mono_type(mono, from, to, strict)
  res <- convert_mono_type_(mono, from, to)

  bad_monos <- mono[is.na(res)]
  if (length(bad_monos) > 0) {
    cli::cli_warn(
      "Some monosaccharides cannot be converted to {.val {to}}: {.val {bad_monos}}.",
      call = rlang::expr(convert_mono_type())
    )
  }

  res
}


valid_from_to_for_convert_mono_type <- function(mono, from, to, strict) {
  tryCatch(
    valid_from_to(from, to, strict),
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


#' Decide the Type of Monosaacharides in a Glycan Structure
#'
#' This function trys to decide the type of monosaccharides in a glycan structure.
#'
#' @details
#' By saying "trys" in the description, it means that the function only
#' checks the first monosaccharide in the graph to decide the type.
#' This is reasonable because the monosaccharides in a glycan structure are always
#' of the same type if it was created with the functions in this package.
#' The validation functions will ensure this.
#'
#' @inheritSection decide_mono_type Three types of monosaccharides
#'
#' @param glycan A glycan structure.
#'
#' @return A character string specifying the monosaccharide type.
#'
#' @examples
#' decide_glycan_mono_type(n_glycan_core(mono_type = "concrete"))
#' decide_glycan_mono_type(n_glycan_core(mono_type = "generic"))
#' decide_glycan_mono_type(n_glycan_core(mono_type = "simple"))
#'
#' @seealso [convert_glycan_mono_type()], [convert_mono_type()], [decide_mono_type()]
#'
#' @export
decide_glycan_mono_type <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_structure")
  decide_glycan_mono_type_impl(glycan)
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
#' For the full list of monosaccharides, use [available_monosaccharides()].
#'
#' @param monos A character string specifying monosaccharide names.
#'
#' @return A character string specifying the monosaccharide types.
#'
#' @examples
#' decide_mono_type(c("Gal", "Hex", "H"))
#'
#' @seealso [convert_glycan_mono_type()], [convert_mono_type()], [decide_glycan_mono_type()]
#'
#' @export
decide_mono_type <- function(monos) {
  checkmate::assert_character(monos)
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


valid_from_to <- function(from, to, strict) {
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
  if (strict && any(from == to)) {
    rlang::abort(
      class = "error_convert_self",
      bad_index = which(from == to)
    )
  }
}


decide_glycan_mono_type_impl <- function(glycan) {
  first_mono <- igraph::vertex_attr(glycan, "mono")[[1]]
  decide_mono_type(first_mono)
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
      call = rlang::expr(convert_glycan_mono_type())
    )
  }
}
