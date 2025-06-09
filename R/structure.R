new_glycan_structure <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  structure(graph, class = c("glycan_structure", "igraph"))
}

validate_glycan_structure <- function(glycan) {
  checkmate::assert_class(glycan, "glycan_structure")
  
  # Check if it is a directed graph
  if (!is_directed_graph(glycan)) {
    rlang::abort("Glycan structure must be directed.")
  }
  
  # Check if it is an out tree
  if (!is_out_tree(glycan)) {
    rlang::abort("Glycan structure must be an out tree.")
  }
  
  # Check if graph has a vertex attribute "mono"
  # This is the monosaccharide name, e.g. "GlcNAc", "Man", etc.
  if (!has_vertex_attrs(glycan, "mono")) {
    rlang::abort("Glycan structure must have a vertex attribute 'mono'.")
  }
  
  # Check if no NA in vertex attribute "mono"
  mono_names <- igraph::vertex_attr(glycan, "mono")
  if (any(is.na(mono_names))) {
    rlang::abort("Glycan structure must have no NA in vertex attribute 'mono'.")
  }
  
  # Check if all monosaccharides are known
  if (!all(is_known_mono(mono_names))) {
    unknown_monos <- unique(igraph::V(glycan)$mono[!is_known_mono(igraph::V(glycan)$mono)])
    msg <- glue::glue("Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}")
    rlang::abort(msg, monos = unknown_monos)
  }
  
  # Check if mixed use of generic and concrete monosaccharides
  if (mix_generic_concrete(mono_names)) {
    rlang::abort("Monosaccharides must be either all generic or all concrete.")
  }
  
  # Check if graph has a vertex attribute "sub"
  # This is the substituent name, e.g. "Ac", "S", "P", or "" (no).
  if (!has_vertex_attrs(glycan, "sub")) {
    rlang::abort("Glycan structure must have a vertex attribute 'sub'.")
  }
  
  # Check if no NA in vertex attribute "sub"
  subs <- igraph::vertex_attr(glycan, "sub")
  if (any(is.na(subs))) {
    rlang::abort("Glycan structure must have no NA in vertex attribute 'sub'.")
  }
  
  # Check if all substituents are valid
  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue("Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}")
    rlang::abort(msg, subs = invalid_subs)
  }
  
  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attrs(glycan, "linkage")) {
    rlang::abort("Glycan structure must have an edge attribute 'linkage'.")
  }
  
  # Check if no NA in edge attribute "linkage"
  linkages <- igraph::edge_attr(glycan, "linkage")
  if (any(is.na(linkages))) {
    rlang::abort("Glycan structure must have no NA in edge attribute 'linkage'.")
  }
  
  # Check if all linkages are valid
  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
  
  # Check if "anomer" attribute exists
  if (is.null(glycan$anomer)) {
    rlang::abort("Glycan structure must have a graph attribute 'anomer'.")
  }
  
  # Check if "anomer" attribute is valid
  if (!valid_anomer(glycan$anomer)) {
    rlang::abort(glue::glue("Invalid anomer: {glycan$anomer}"))
  }
  
  # Check if "alditol" attribute exists
  if (is.null(glycan$alditol)) {
    rlang::abort("Glycan structure must have a graph attribute 'alditol'.")
  }
  
  # Check if "alditol" attribute is valid
  if (!is.logical(glycan$alditol)) {
    rlang::abort("Glycan structure attribute 'alditol' must be logical.")
  }

  glycan
}

ensure_name_vertex_attr <- function(glycan) {
  if (!("name" %in% igraph::vertex_attr_names(glycan))) {
    names <- as.character(seq_len(igraph::vcount(glycan)))
    glycan <- igraph::set_vertex_attr(glycan, "name", value = names)
  }
  glycan
}

#' Convert igraph Graph to Glycan Structure
#'
#' @description
#' A glycan structure is a subclass of igraph graph with additional constraints.
#' This function checks these constraints and append the S3 class `glycan_structure`
#' to the graph object.
#'
#' Generally, **it is not recommended** to create a glycan structure using this
#' function manually.
#' Use the [glyparse](https://github.com/glycoverse/glyparse) package
#' to generate a glycan structure from a structure string.
#'
#' A glycan structure is a directed modeling of a glycan structure,
#' where nodes are monosaccharides and edges are linkages.
#'
#' @details
#' # S3 Classes
#'
#' A glycan structure is merely an igraph graph with an additional S3 class `glycan_structure`.
#' Therefore, `sloop::s3_class()` of a glycan structure is 
#' `c("glycan_structure", "igraph")`.
#'
#' Constraints:
#' - The graph must be directed and an outward tree (reducing end as root).
#' - The graph must have a graph attribute `anomer`, in the form of "a1".
#' - The graph must have a graph attribute `alditol`, a logical value.
#' - The graph must have a vertex attribute `mono` for monosaccharide names.
#' - The graph must have a vertex attribute `sub` for substituents.
#' - The graph must have an edge attribute `linkage` for linkages.
#' - Monosaccharide name must be known, either generic (Hex, HexNAc, etc.)
#'   or concrete (Glc, Gal, etc.), but not a mixture of both.
#'   NA is not allowed.
#' - Substituent must be valid, in the form of "xY", where x is position
#'   and Y is substituent name, e.g. "2Ac", "3S", etc.
#' - Linkages must be valid, in the form of "a/bX-Y", where X and Y are integers,
#'   e.g. "b1-4", "a2-3", etc.
#'   NA is not allowed.
#'
#' @param graph An igraph graph object.
#' @param x An object to check.
#'
#' @return A glycan structure object.
#'
#' @examples
#' library(igraph)
#'
#' # A simple glycan structure: GlcNAc(b1-4)GlcNAc
#' graph <- make_graph(~ 1-+2)  # 1 and 2 are monosaccharides
#' V(graph)$mono <- c("GlcNAc", "GlcNAc")
#' V(graph)$sub <- ""
#' E(graph)$linkage <- "b1-4"
#' graph$anomer <- "a1"
#' graph$alditol <- FALSE
#' glycan_structure(graph)
#'
#' @importFrom magrittr %>%
#' @export
glycan_structure <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  
  graph %>%
    new_glycan_structure() %>%
    validate_glycan_structure() %>%
    ensure_name_vertex_attr()
}


#' @export 
#' @rdname glycan_structure
is_glycan_structure <- function(x) {
  inherits(x, "glycan_structure")
}
