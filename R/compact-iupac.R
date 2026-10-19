# Private character-input backend. Public graph APIs retain their validation,
# arbitrary attributes, and trusted-input contracts.
.compact_bliss_labels <- function(edges, colors) {
  graph <- .compact_blank_graph_native(
    edges,
    length(colors),
    FALSE,
    .compact_graph_template()
  )
  if (is.null(graph)) {
    graph <- igraph::make_graph(edges, n = length(colors), directed = FALSE)
  }
  as.integer(
    igraph::canonical_permutation(graph, colors = colors, sh = "fm")$labeling
  )
}

.compact_iupac_arrays <- function(x, progress = NULL) {
  residues <- c(monosaccharides$concrete, monosaccharides$generic)
  residues <- unique(residues[!is.na(residues)])
  .compact_parse_native(
    unname(x),
    residues,
    infer_anomer_pos(residues),
    names(unusual_configuration_monosaccharides),
    unname(unusual_configuration_monosaccharides),
    base::order,
    .compact_bliss_labels,
    byte_order = identical(Sys.getlocale("LC_COLLATE"), "C"),
    progress = .structure_progress_stage(
      progress,
      "Parsing unique structures",
      length(x)
    )
  )
}

.compact_build_graph <- function(a) {
  graph <- igraph::make_graph(a$edges, n = length(a$mono), directed = TRUE)
  igraph::vertex_attr(graph) <- list(
    name = as.character(seq_along(a$mono)),
    mono = a$mono,
    sub = a$sub
  )
  igraph::edge_attr(graph) <- list(linkage = a$linkage)
  igraph::graph_attr(graph) <- list(anomer = a$anomer, alditol = a$alditol)
  if (!is.null(a$floating_parts)) {
    graph <- igraph::set_graph_attr(graph, "floating_parts", a$floating_parts)
  }
  if (!is.null(a$floating_substituents)) {
    graph <- igraph::set_graph_attr(
      graph,
      "floating_substituents",
      a$floating_substituents
    )
  }
  graph
}

.compact_structure_from_arrays <- function(
  x,
  unique_x,
  arrays,
  progress = NULL
) {
  keys <- vapply(arrays, `[[`, character(1), "iupac")
  unique_keys <- !duplicated(keys)
  graphs <- .compact_build_graphs(arrays[unique_keys], progress)
  names(graphs) <- keys[unique_keys]
  # This private constructor has historically returned unnamed values; the
  # existing vctrs method restores names, including its empty-input semantics.
  new_glycan_structure(unname(keys[match(x, unique_x)]), graphs)
}

.compact_iupac_outcomes <- function(x, progress = NULL) {
  arrays <- .compact_iupac_arrays(x, progress)
  graph_cache <- new.env(hash = TRUE, parent = emptyenv())
  .structure_progress_map(
    seq_along(arrays),
    function(i) {
      tryCatch(
        {
          a <- arrays[[i]]
          if (!identical(a$status, "ok")) {
            # Only the failing element takes the reference path under recovery.
            # This preserves its condition class, message, and normalization rules.
            return(process_glycan_structure_element(.parse_iupac_condensed_single(x[[
              i
            ]])))
          }
          key <- a$iupac
          if (!exists(key, graph_cache, inherits = FALSE)) {
            graph_cache[[key]] <- .compact_build_graph(a)
          }
          list(graph = graph_cache[[key]], iupac = key)
        },
        error = \(cnd) cnd
      )
    },
    progress,
    "Building structures"
  )
}

# Cache only the empty format template, never parsed structures or user inputs.
.compact_graph_template <- local({
  model <- NULL
  pid <- NULL
  function() {
    current_pid <- Sys.getpid()
    if (is.null(model) || !identical(pid, current_pid)) {
      pid <<- current_pid
      candidate <- igraph::make_empty_graph(n = 0, directed = TRUE)
      model <<- if (is.null(.compact_build_graphs_native(list(), candidate))) {
        FALSE
      } else {
        candidate
      }
    }
    if (identical(model, FALSE)) NULL else model
  }
})

.compact_build_graphs <- function(arrays, progress = NULL) {
  model <- .compact_graph_template()
  if (is.null(model)) {
    return(.structure_progress_map(
      arrays,
      .compact_build_graph,
      progress,
      "Building structures"
    ))
  }
  graphs <- .compact_build_graphs_native(
    arrays,
    model,
    .structure_progress_stage(progress, "Building structures", length(arrays))
  )
  if (!is.null(graphs)) {
    return(graphs)
  }
  .structure_progress_map(
    arrays,
    .compact_build_graph,
    progress,
    "Building structures"
  )
}
