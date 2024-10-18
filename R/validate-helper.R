# Is the graph directed?
is_directed_graph <- function(graph) {
  igraph::is_directed(graph)
}


# Is the graph an outward tree?
is_out_tree <- function(graph) {
  igraph::is_tree(graph, mode = "out")
}


# Does the graph have these vertex attributes?
has_vertex_attrs <- function(graph, attrs) {
  all(attrs %in% igraph::vertex_attr_names(graph))
}


# Does the graph have these edge attributes?
has_edge_attrs <- function(graph, attrs) {
  all(attrs %in% igraph::edge_attr_names(graph))
}


# Are all monosaaccharides known?
is_known_mono <- function(monos) {
  known_monos <- c(unique(monosaccharides$generic), monosaccharides$concrete)
  monos %in% known_monos
}


# Are generic and concrete monosaccharides mixed?
mix_generic_concrete <- function(monos) {
  generic <- unique(monosaccharides$generic)
  concrete <- monosaccharides$concrete
  any(monos %in% generic) && any(monos %in% concrete)
}


# Returns a logical vector indicating whether the linkages are valid
valid_linkages <- function(linkages) {
  stringr::str_detect(linkages, "^[ab]\\d+-\\d+$")
}
