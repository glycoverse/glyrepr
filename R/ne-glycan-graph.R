# New Node-Edge Glycan Graph
new_ne_glycan_graph <- function(graph) {
  stopifnot(igraph::is_igraph(graph))
  structure(graph, class = c("ne_glycan_graph", "glycan_graph", "igraph"))
}


validate_ne_glycan_graph <- function(glycan) {
  stopifnot(inherits(glycan, "ne_glycan_graph"))
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
  if (na_in_vertex_attribute(glycan, "mono")) {
    rlang::abort("Glycan graph must have no NA in vertex attribute 'mono'.")
  }
  # Check if all monosaccharides are known
  if (!all(is_known_mono(igraph::vertex_attr(glycan, "mono")))) {
    unknown_monos <- unique(igraph::V(glycan)$mono[!is_known_mono(igraph::V(glycan)$mono)])
    msg <- glue::glue("Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}")
    rlang::abort(msg, monos = unknown_monos)
  }
  # Check if mixed use of generic and concrete monosaccharides
  if (mix_generic_concrete(igraph::V(glycan)$mono)) {
    rlang::abort("Monosaccharides must be either all generic or all concrete.")
  }
  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attrs(glycan, "linkage")) {
    rlang::abort("Glycan graph must have an edge attribute 'linkage'.")
  }
  # Check if all linkages are valid
  if (!all(valid_linkages(igraph::edge_attr(glycan, "linkage")))) {
    invalid_linkages <- unique(igraph::E(glycan)$linkage[!valid_linkages(igraph::E(glycan)$linkage)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
}
