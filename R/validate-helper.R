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


.unique_no_na <- function(x) unique(x[!is.na(x)])


# Are all monosaaccharides known?
is_known_mono <- function(monos) {
  known_monos <- c(
    .unique_no_na(monosaccharides$simple),
    .unique_no_na(monosaccharides$generic),
    monosaccharides$concrete
  )
  monos %in% known_monos
}


# Are generic and concrete monosaccharides mixed?
mix_generic_concrete <- function(monos) {
  has_simple <- any(monos %in% .unique_no_na(monosaccharides$simple))
  has_generic <- any(monos %in% .unique_no_na(monosaccharides$generic))
  has_concrete <- any(monos %in% monosaccharides$concrete)
  sum(as.integer(c(has_simple, has_generic, has_concrete))) > 1
}


# Returns a logical vector indicating whether the linkages are valid
valid_linkages <- function(linkages) {
  stringr::str_detect(linkages, "^[ab]\\d+-\\d+$")
}
