new_glycan_graph <- function(graph) {
  stopifnot(igraph::is_igraph(graph))
  structure(graph, class = c("glycan_graph", "igraph"))
}


validate_glycan_graph <- function(glycan) {
  stopifnot(inherits(glycan, "glycan_graph"))
  # Check if it is a directed graph
  if (!igraph::is_directed(glycan)) {
    rlang::abort("Glycan graph must be directed.")
  }
  # Check if it is an out tree
  if (!igraph::is_tree(glycan, mode = "out")) {
    rlang::abort("Glycan graph must be an out tree.")
  }
  # Check if graph has a vertex attribute "mono"
  # This is the monosaccharide name, e.g. "GlcNAc", "Man", etc.
  if (is.null(igraph::V(glycan)$mono)) {
    rlang::abort("Glycan graph must have a vertex attribute 'mono'.")
  }
  # Check if no NA in vertex attribute "mono"
  if (any(is.na(igraph::V(glycan)$mono))) {
    rlang::abort("Glycan graph must have no NA in vertex attribute 'mono'.")
  }
  # Check if all monosaccharides are known
  if (!all(is_known_mono(igraph::V(glycan)$mono))) {
    unknown_monos <- unique(igraph::V(glycan)$mono[!is_known_mono(igraph::V(glycan)$mono)])
    abort_unknown_mono(unknown_monos)
  }
  # Check if graph has an edge attribute "linkage"
  if (is.null(igraph::E(glycan)$linkage)) {
    rlang::abort("Glycan graph must have an edge attribute 'linkage'.")
  }
  # Check if all linkages are valid
  if (!all(is_valid_linkage(igraph::E(glycan)$linkage))) {
    invalid_linkages <- unique(igraph::E(glycan)$linkage[!is_valid_linkage(igraph::E(glycan)$linkage)])
    abort_invalid_linkage(invalid_linkages)
  }
}


abort_unknown_mono <- function(monos) {
  monos_str <- paste0(monos, collapse = ", ")
  msg <- glue::glue("Unknown monosaccharide: {monos_str}")
  rlang::abort("error_bad_mono", message = msg, monos = monos)
}


is_known_mono <- function(monos) {
  known_monos <- c(unique(monosaccharides$generic), monosaccharides$concrete)
  monos %in% known_monos
}


abort_invalid_linkage <- function(linkages) {
  linkages_str <- paste0(linkages, collapse = ", ")
  msg <- glue::glue("Invalid linkage: {linkages_str}")
  rlang::abort("error_bad_linkage", message = msg, linkages = linkages)
}


is_valid_linkage <- function(linkage) {
  stringr::str_detect(linkage, "^[ab]\\d+-\\d+$")
}
