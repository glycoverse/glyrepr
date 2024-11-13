# New Dual-Node Glycan Graph
new_dn_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  structure(graph, class = c("dn_glycan_graph", "glycan_graph", "igraph"))
}


validate_dn_glycan_graph <- function(glycan) {
  checkmate::assert_class(glycan, "dn_glycan_graph")
  # Check if it is a directed graph
  if (!is_directed_graph(glycan)) {
    rlang::abort("Glycan graph must be directed.")
  }
  # Check if it is an out tree
  if (!is_out_tree(glycan)) {
    rlang::abort("Glycan graph must be an out tree.")
  }
  # Check if graph has the right vertex attributes
  if (!has_vertex_attrs(glycan, c("type", "mono", "sub", "linkage"))) {
    rlang::abort("Glycan graph must have vertex attributes 'type', 'mono', 'sub', and 'linkage'.")
  }
  # Check if no NA in "type" attribute
  if (any(is.na(igraph::vertex_attr(glycan, "type")))) {
    rlang::abort("Glycan graph must have no NA in 'type' attribute.")
  }
  # Check if only "mono" and "linkage" are in "type" attribute
  if (!all(igraph::V(glycan)$type %in% c("mono", "linkage"))) {
    rlang::abort("The 'type' of a node could either be 'mono' or 'linkage'.")
  }
  # Check if "mono" and "linkage" nodes are alternating
  if (!check_alternating(glycan)) {
    rlang::abort("The 'mono' and 'linkage' nodes must be alternating.")
  }
  # Check if no NA in "mono" attribute for "mono" nodes
  mono_names <- get_vertex_attr(glycan, "mono", "mono")
  if (any(is.na(mono_names))) {
    rlang::abort("Mono nodes must have no NA in 'mono' attribute.")
  }
  # Check if all monosaccharides are known
  if (!all(is_known_mono(mono_names))) {
    unknown_monos <- mono_names[!is_known_mono(mono_names)]
    msg <- glue::glue("Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}")
    rlang::abort(msg, monos = unknown_monos)
  }
  # Check if mixed use of generic and concrete monosaccharides
  if (mix_generic_concrete(mono_names)) {
    rlang::abort("Mono nodes must not mix generic and concrete monosaccharides.")
  }
  # Chekc if no NA in "sub" attribute
  subs <- get_vertex_attr(glycan, "mono", "sub")
  if (any(is.na(subs))) {
    rlang::abort("Mono nodes must have no NA in 'sub' attribute.")
  }
  # Check if all substituents are valid
  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue("Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}")
    rlang::abort(msg, subs = invalid_subs)
  }
  # Check if no NA in "linkage" attribute
  linkages <- get_vertex_attr(glycan, "linkage", "linkage")
  if (any(is.na(linkages))) {
    rlang::abort("Linkage nodes must have no NA in 'linkage' attribute.")
  }
  # Check if all linkages are valid
  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue("Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}")
    rlang::abort(msg, linkages = invalid_linkages)
  }
  # Check if no exposed "linkage" nodes
  linkage_nodes <- which(igraph::V(glycan)$type == "linkage")
  out_degrees <- igraph::degree(glycan, v = linkage_nodes, mode = "out")
  if (any(out_degrees == 0)) {
    rlang::abort("Linkage nodes must not have an out-degree of zero (no exposed linkages)")
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


check_alternating <- function(graph) {
  edge_list <- igraph::as_edgelist(graph, names = FALSE)
  types <- igraph::V(graph)$type
  u <- edge_list[,1]
  v <- edge_list[,2]
  all(types[u] != types[v])
}


get_vertex_attr <- function(graph, type, attr) {
  igraph::vertex_attr(graph, attr)[igraph::vertex_attr(graph, "type") == type]
}
