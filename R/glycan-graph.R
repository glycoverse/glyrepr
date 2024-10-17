new_glycan_graph <- function(graph) {
  stopifnot(igraph::is_igraph(graph))
  structure(graph, class = c("glycan_graph", "igraph"))
}


validate_glycan_graph <- function(glycan) {
  stopifnot(inherits(glycan, "glycan_graph"))
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
  if (!has_vertex_attribute(glycan, "mono")) {
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
  # Check if graph has an edge attribute "linkage"
  if (!has_edge_attribute(glycan, "linkage")) {
    rlang::abort("Glycan graph must have an edge attribute 'linkage'.")
  }
  # Check if all linkages are valid
  if (!all(is_valid_linkage(igraph::edge_attr(glycan, "linkage")))) {
    invalid_linkages <- unique(igraph::E(glycan)$linkage[!is_valid_linkage(igraph::E(glycan)$linkage)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
}


is_directed_graph <- function(graph) {
  igraph::is_directed(graph)
}


is_out_tree <- function(graph) {
  igraph::is_tree(graph, mode = "out")
}


has_vertex_attribute <- function(graph, attr) {
  attr %in% igraph::vertex_attr_names(graph)
}


na_in_vertex_attribute <- function(graph, attr) {
  any(is.na(igraph::vertex_attr(graph, attr)))
}


has_edge_attribute <- function(graph, attr) {
  attr %in% igraph::edge_attr_names(graph)
}


is_known_mono <- function(monos) {
  known_monos <- c(unique(monosaccharides$generic), monosaccharides$concrete)
  monos %in% known_monos
}


is_valid_linkage <- function(linkage) {
  stringr::str_detect(linkage, "^[ab]\\d+-\\d+$")
}
