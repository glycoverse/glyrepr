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

.compact_iupac_outcomes <- function(x, progress = NULL, arrays = NULL) {
  if (is.null(arrays)) {
    arrays <- .compact_iupac_arrays(x, progress)
  }
  ok <- vapply(arrays, function(a) identical(a$status, "ok"), logical(1))
  keys <- vapply(arrays[ok], `[[`, character(1), "iupac")
  unique_keys <- !duplicated(keys)
  graphs <- .compact_build_graphs(arrays[ok][unique_keys], progress)
  outcomes <- vector("list", length(arrays))
  outcomes[ok] <- Map(
    function(graph, iupac) list(graph = graph, iupac = iupac),
    graphs[match(keys, keys[unique_keys])],
    keys
  )
  outcomes[!ok] <- .structure_progress_map(
    which(!ok),
    function(i) {
      tryCatch(
        {
          # Only unsupported or failing elements need reference diagnostics.
          process_glycan_structure_element(.parse_iupac_condensed_single(x[[
            i
          ]]))
        },
        error = \(cnd) cnd
      )
    },
    progress,
    "Building structures"
  )
  outcomes
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
        # Keep ordinary vectors and no native pointer in the cached closure.
        # covr may serialize this environment during R shutdown, after igraph's
        # finalizer has freed a live graph's pointer-backed endpoint vectors.
        .compact_blank_graph_native(integer(), 0L, TRUE, candidate)
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
