args <- commandArgs(TRUE)
.libPaths(c(args[[1]], .libPaths()))
library(glyrepr)
options(cli.num_colors = 1L, cli.width = 80L)
output <- args[[2]]
graph_signature <- function(g) {
  list(
    directed = igraph::is_directed(g),
    edges = unname(igraph::as_edgelist(g, names = FALSE)),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
normalize <- function(x) {
  if (inherits(x, "igraph")) {
    return(graph_signature(x))
  }
  if (inherits(x, "glyrepr_structure")) {
    return(list(
      class = class(x),
      values = as.character(x),
      names = names(x),
      graphs = lapply(attr(x, "graphs"), graph_signature)
    ))
  }
  if (is.list(x)) {
    out <- lapply(x, normalize)
    attributes(out) <- attributes(x)
    return(out)
  }
  x
}
condition_signature <- function(cnd) {
  list(
    class = class(cnd),
    message = conditionMessage(cnd),
    call = deparse(conditionCall(cnd)),
    parent = if (is.null(cnd$parent)) NULL else condition_signature(cnd$parent)
  )
}
capture <- function(expr, env) {
  warnings <- list()
  result <- tryCatch(
    withCallingHandlers(eval(expr, env), warning = function(cnd) {
      warnings[[length(warnings) + 1L]] <<- condition_signature(cnd)
      invokeRestart("muffleWarning")
    }),
    error = function(cnd) cnd
  )
  list(
    result = if (inherits(result, "error")) {
      condition_signature(result)
    } else {
      normalize(result)
    },
    warnings = warnings
  )
}
inputs <- list(
  named = c(
    a = "Gal6S",
    b = "{6S|1}Gal",
    missing = NA_character_,
    floating = "{Neu5Ac(a2-3)}Gal(b1-4)Glc"
  ),
  empty = character(),
  named_empty = setNames(character(), character()),
  missing = c(a = NA_character_, b = NA_character_),
  unusual_names = setNames(c("Glc", NA, "Glc"), c("", NA, "duplicate")),
  invalid = c(a = "Glc", bad = "bad", again = "bad", missing = NA_character_),
  invalid_only = c(bad = "bad", empty = ""),
  parse_priority = c("Gal(b1-4)[Man(a1-4)]Glc", "bad"),
  graph_validation = c("Gal(b1-4)[Man(a1-4)]Glc", "Glc"),
  floating_invalid = c("{Gal(a1-3)|1}Glc", "{6S|1}Glc", "{6S|1}Glc6S"),
  list_strings = list("Glc"),
  logical = NA,
  null = NULL,
  numeric = 1,
  factor = factor("Glc"),
  matrix = matrix(c("Gal", "Glc"), nrow = 1)
)
g <- get_structure_graphs(as_glycan_structure("Gal(b1-4)Glc"))
g <- igraph::set_graph_attr(g, "source", list(id = "retained"))
g <- igraph::set_vertex_attr(g, "label", value = c("leaf", "root"))
g <- igraph::set_edge_attr(g, "weight", value = 2.5)
bad <- igraph::set_vertex_attr(g, "mono", value = c("Invalid", "Glc"))
inputs$graph <- g
inputs$list_graphs <- list(a = g, b = g)
inputs$list_mixed <- list(a = g, missing = NA, invalid = bad, missing2 = NULL)
inputs$bad_graph <- bad
inputs$vector <- as_glycan_structure(inputs$named)
cases <- list()
for (name in names(inputs)) {
  for (policy in c("error", "na")) {
    env <- list2env(
      list(x = inputs[[name]], policy = policy),
      parent = globalenv()
    )
    cases[[paste(name, policy)]] <- capture(
      quote(as_glycan_structure(x, on_failure = policy)),
      env
    )
  }
}
env <- new.env(parent = globalenv())
env$x <- as_glycan_structure(inputs$named)
env$g <- g
env$bad <- bad
expressions <- list(
  default = quote(as_glycan_structure(c("bad", "Glc"))),
  invalid_policy = quote(as_glycan_structure("Glc", on_failure = "skip")),
  partial_policy = quote(as_glycan_structure("Glc", on_failure = "n")),
  cast_character = quote(vctrs::vec_cast(
    c(one = "Glc", two = NA_character_),
    glycan_structure()
  )),
  cast_graph = quote(vctrs::vec_cast(g, glycan_structure())),
  cast_list = quote(vctrs::vec_cast(list(g, g), glycan_structure())),
  cast_back = quote(vctrs::vec_cast(x, character())),
  constructor = quote(glycan_structure(g, NA, g)),
  constructor_error = quote(glycan_structure("Glc")),
  concatenate = quote(c(x, x[1], NA)),
  slice = quote(x[c(4, 1, NA)]),
  element = quote(x[[1]]),
  unique = quote(unique(x)),
  list = quote(as.list(x)),
  format = quote(format(x)),
  table_roundtrip = quote(structure_from_tibbles(
    structure_nodes(x),
    structure_edges(x),
    get_anomer(x),
    structure_floating_parts(x),
    structure_floating_substituents(x),
    get_alditol(x)
  )),
  trusted_constructor = quote(new_glycan_structure(
    as.character(x),
    attr(x, "graphs")
  )),
  validated_graph = quote(validate_glycan_graph(g)),
  canonical_graph = quote(canonicalize_glycan_graph(g)),
  graph_iupac = quote(graph_to_iupac(g)),
  invalid_graph = quote(validate_glycan_graph(bad)),
  graph_list_validation = quote(validate_glycan_graph_vector(list(g))),
  mapped = quote(smap_structure(x, identity)),
  mapped_chr = quote(smap_chr(x, graph_to_iupac)),
  enumerate = quote(enumerate_floating_localizations(x[4], max_variants = 256)),
  enumerate_graph = quote(enumerate_floating_graph_localizations(get_structure_graphs(x[
    4
  ]))),
  localize = quote(localize_floating_parts(
    x[4],
    enumerate_floating_localizations(x[4])$assignments[[1]]
  ))
)
for (name in names(expressions)) {
  cases[[name]] <- capture(expressions[[name]], env)
}
accessors <- c(
  "get_structure_graphs",
  "structure_to_iupac",
  "structure_nodes",
  "structure_edges",
  "structure_floating_parts",
  "structure_floating_substituents",
  "structure_candidate_edges",
  "structure_component_membership",
  "structure_floating_candidates",
  "get_anomer",
  "get_alditol",
  "get_mono_type",
  "get_structure_level",
  "has_linkages",
  "has_floating_parts",
  "has_floating_substituents",
  "count_mono",
  "as_glycan_composition",
  "convert_to_generic",
  "remove_linkages",
  "remove_substituents",
  "fill_anomer_pos"
)
for (fn in accessors) {
  for (kind in c("x", "g")) {
    cases[[paste(fn, kind)]] <- capture(
      as.call(list(as.name(fn), as.name(kind))),
      env
    )
  }
}
api <- lapply(sort(getNamespaceExports("glyrepr")), function(n) {
  value <- getExportedValue("glyrepr", n)
  if (is.function(value)) formals(value) else class(value)
})
names(api) <- sort(getNamespaceExports("glyrepr"))
saveRDS(
  list(api = api, cases = cases, session = capture.output(sessionInfo())),
  output
)
cat(
  "Captured",
  length(api),
  "exported signatures and",
  length(cases),
  "API cases\n"
)
