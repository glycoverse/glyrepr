# New Node-Edge Glycan Graph
new_ne_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  structure(graph, class = c("ne_glycan_graph", "glycan_graph", "igraph"))
}


validate_ne_glycan_graph <- function(glycan) {
  checkmate::assert_class(glycan, "ne_glycan_graph")
  # Check if it is a directed graph
  if (!is_directed_graph(glycan)) {
    rlang::abort("Glycan graph must be directed.")
  }
  # Check if it is an out tree
  if (!is_out_tree(glycan)) {
    rlang::abort("Glycan graph must be an out tree.")
  }
  # Check if graph has a vertex attribute "mono"
  # This is the monosaccharide name, e.g. "GlcNAc", "Man", etc.
  if (!has_vertex_attrs(glycan, "mono")) {
    rlang::abort("Glycan graph must have a vertex attribute 'mono'.")
  }
  # Check if no NA in vertex attribute "mono"
  mono_names <- igraph::vertex_attr(glycan, "mono")
  if (any(is.na(mono_names))) {
    rlang::abort("Glycan graph must have no NA in vertex attribute 'mono'.")
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
    rlang::abort("Glycan graph must have a vertex attribute 'sub'.")
  }
  # Check if no NA in vertex attribute "sub"
  subs <- igraph::vertex_attr(glycan, "sub")
  if (any(is.na(subs))) {
    rlang::abort("Glycan graph must have no NA in vertex attribute 'sub'.")
  }
  # Check if all substituents are valid
  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue("Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}")
    rlang::abort(msg, subs = invalid_subs)
  }
  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attrs(glycan, "linkage")) {
    rlang::abort("Glycan graph must have an edge attribute 'linkage'.")
  }
  # Check if no NA in edge attribute "linkage"
  linkages <- igraph::edge_attr(glycan, "linkage")
  if (any(is.na(linkages))) {
    rlang::abort("Glycan graph must have no NA in edge attribute 'linkage'.")
  }
  # Check if all linkages are valid
  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
  # Check if "anomer" attribute exists
  if (is.null(glycan$anomer)) {
    rlang::abort("Glycan graph must have a graph attribute 'anomer'.")
  }
  # Check if "anomer" attribute is valid
  if (!valid_anomer(glycan$anomer)) {
    rlang::abort(glue::glue("Invalid anomer: {glycan$anomer}"))
  }
  # Check if "alditol" attribute exists
  if (is.null(glycan$alditol)) {
    rlang::abort("Glycan graph must have a graph attribute 'alditol'.")
  }
  # Check if "alditol" attribute is valid
  if (!is.logical(glycan$alditol)) {
    rlang::abort("Glycan graph attribute 'alditol' must be logical.")
  }

  glycan
}
