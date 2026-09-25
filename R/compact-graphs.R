# Graph bridge uses only public igraph accessors. A NULL outcome requests the
# existing R path, including its condition classes, messages and failure order.
.compact_graph_record <- function(graph, validate = FALSE) {
  if (!inherits(graph, "igraph") || !igraph::is_directed(graph)) {
    return(NULL)
  }
  attrs <- igraph::graph_attr(graph)
  vertices <- igraph::vertex_attr(graph)
  edges <- igraph::edge_attr(graph)
  if (
    is.null(vertices$mono) ||
      is.null(vertices$sub) ||
      is.null(edges$linkage) ||
      is.null(attrs$anomer)
  ) {
    return(NULL)
  }
  # Attribute-bearing standard vectors use the reference path; the native
  # representation intentionally stores plain scalar values only.
  standard <- list(
    vertices$mono,
    vertices$sub,
    edges$linkage,
    attrs$anomer,
    attrs$alditol
  )
  if (any(vapply(standard, function(x) !is.null(attributes(x)), logical(1)))) {
    return(NULL)
  }
  a <- list(
    mono = vertices$mono,
    sub = vertices$sub,
    edges = as.integer(t(igraph::as_edgelist(graph, names = FALSE))),
    linkage = edges$linkage,
    anomer = attrs$anomer,
    alditol = if (is.null(attrs$alditol)) FALSE else attrs$alditol
  )
  if (!is.null(attrs$floating_parts)) {
    a$floating_parts <- normalize_floating_parts(graph)
  }
  if (!is.null(attrs$floating_substituents)) {
    a$floating_substituents <- normalize_floating_substituents(graph)
  }
  # The compact linkage representation has a single donor-position character.
  # Unusual valid graph linkages keep the established R implementation.
  links <- c(a$linkage, vapply(a$floating_parts, `[[`, character(1), "linkage"))
  if (anyNA(links) || any(!grepl("^[ab?][12?]-([1-9](/[1-9])*|[?])$", links))) {
    return(NULL)
  }
  if (validate) {
    return(.validate_structure_arrays(a))
  }
  # Trusted public modes skip chemical validation. These representation checks
  # still protect the native importer from malformed low-level input.
  if (
    !is.character(a$mono) ||
      !length(a$mono) ||
      anyNA(a$mono) ||
      !is.character(a$sub) ||
      length(a$sub) != length(a$mono) ||
      anyNA(a$sub) ||
      !is.character(a$linkage) ||
      length(a$linkage) * 2L != length(a$edges) ||
      !is.character(a$anomer) ||
      length(a$anomer) != 1L ||
      anyNA(a$anomer) ||
      !is.logical(a$alditol) ||
      length(a$alditol) != 1L ||
      anyNA(a$alditol)
  ) {
    return(NULL)
  }
  a
}

.compact_graph_results <- function(
  graphs,
  mode = "canonicalize",
  validate = FALSE,
  operation = "",
  progress = NULL
) {
  records <- .structure_progress_map(
    graphs,
    function(graph) {
      tryCatch(.compact_graph_record(graph, validate), error = function(e) NULL)
    },
    progress,
    "Validating graphs"
  )
  from <- to <- character()
  if (operation == "generic") {
    from <- unique(monosaccharides$concrete[!is.na(monosaccharides$concrete)])
    to <- convert_mono_type_impl(from)
  } else if (operation == "fill_anomer_pos") {
    from <- unique(unlist(lapply(records, `[[`, "mono"), use.names = FALSE))
    to <- as.character(infer_anomer_pos(from))
  }
  .compact_graphs_native(
    records,
    mode,
    validate,
    operation,
    from,
    to,
    base::order,
    .compact_bliss_labels,
    identical(Sys.getlocale("LC_COLLATE"), "C"),
    progress = .structure_progress_stage(
      progress,
      "Canonicalizing graphs",
      length(graphs)
    )
  )
}

.compact_graph_ok <- function(result) {
  !is.null(result) && identical(result$status, "ok")
}

.compact_restore_graph <- function(graph, a) {
  vertex_attrs <- igraph::vertex_attr(graph)
  edge_attrs <- igraph::edge_attr(graph)
  graph_attrs <- igraph::graph_attr(graph)
  vertex_order <- a$vertex_order
  edge_order <- a$edge_order
  if (
    identical(vertex_order, seq_len(igraph::vcount(graph))) &&
      identical(edge_order, seq_len(igraph::ecount(graph))) &&
      !any(c("floating_parts", "floating_substituents") %in% names(graph_attrs))
  ) {
    if (!identical(vertex_attrs$name, as.character(vertex_order))) {
      graph <- igraph::set_vertex_attr(
        graph,
        "name",
        value = as.character(vertex_order)
      )
    }
    if (!identical(vertex_attrs$mono, a$mono)) {
      graph <- igraph::set_vertex_attr(graph, "mono", value = a$mono)
    }
    if (!identical(vertex_attrs$sub, a$sub)) {
      graph <- igraph::set_vertex_attr(graph, "sub", value = a$sub)
    }
    if (!identical(edge_attrs$linkage, a$linkage)) {
      graph <- igraph::set_edge_attr(graph, "linkage", value = a$linkage)
    }
    if (!identical(graph_attrs$anomer, a$anomer)) {
      graph <- igraph::set_graph_attr(graph, "anomer", value = a$anomer)
    }
    return(normalize_alditol_attr(graph))
  }
  vertex_attrs <- lapply(vertex_attrs, `[`, vertex_order)
  vertex_attrs$name <- as.character(seq_along(a$mono))
  vertex_attrs$mono <- a$mono
  vertex_attrs$sub <- a$sub
  vertex_attrs <- vertex_attrs[c("name", setdiff(names(vertex_attrs), "name"))]
  edge_attrs <- lapply(edge_attrs, `[`, edge_order)
  edge_attrs$linkage <- a$linkage
  graph_attrs$anomer <- a$anomer
  graph_attrs$alditol <- a$alditol
  graph_attrs$floating_parts <- NULL
  graph_attrs$floating_substituents <- NULL
  if (!is.null(a$floating_parts)) {
    graph_attrs$floating_parts <- a$floating_parts
  }
  if (!is.null(a$floating_substituents)) {
    graph_attrs$floating_substituents <- a$floating_substituents
  }
  result <- igraph::make_graph(a$edges, n = length(a$mono), directed = TRUE)
  igraph::vertex_attr(result) <- vertex_attrs
  igraph::edge_attr(result) <- edge_attrs
  igraph::graph_attr(result) <- graph_attrs
  result
}

.compact_transform_structure <- function(x, operation) {
  iupacs <- glycan_structure_iupac_data(x)
  keys <- unique(iupacs[!is.na(iupacs)])
  if (!length(keys)) {
    return(x)
  }
  graphs <- attr(x, "graphs")[keys]
  convertible <- monosaccharides$concrete[
    !is.na(monosaccharides$concrete) &
      !is.na(monosaccharides$generic) &
      monosaccharides$concrete != monosaccharides$generic
  ]
  changed <- vapply(
    graphs,
    function(g) {
      if (operation == "generic") {
        return(any(igraph::vertex_attr(g, "mono") %in% convertible))
      }
      if (operation == "remove_substituents") {
        return(
          any(nzchar(igraph::vertex_attr(g, "sub"))) ||
            length(igraph::graph_attr(g, "floating_substituents")) > 0L
        )
      }
      links <- c(
        igraph::graph_attr(g, "anomer"),
        igraph::edge_attr(g, "linkage"),
        vapply(
          igraph::graph_attr(g, "floating_parts"),
          `[[`,
          character(1),
          "linkage"
        )
      )
      if (operation == "fill_anomer_pos") {
        return(any(substr(links, 2L, 2L) == "?"))
      }
      !identical(links[[1]], "??") || any(links[-1] != "??-?")
    },
    logical(1)
  )
  if (!any(changed)) {
    return(x)
  }
  native <- .compact_graph_results(graphs[changed], operation = operation)
  if (!all(vapply(native, .compact_graph_ok, logical(1)))) {
    return(NULL)
  }
  new_keys <- unname(keys)
  new_keys[changed] <- vapply(native, `[[`, character(1), "iupac")
  if (identical(unname(keys), new_keys)) {
    return(x)
  }
  keep <- !duplicated(new_keys)
  restore <- which(changed & keep)
  for (i in restore) {
    graphs[[i]] <- .compact_restore_graph(
      graphs[[i]],
      native[[match(i, which(changed))]]
    )
  }
  graphs <- graphs[keep]
  names(graphs) <- new_keys[keep]
  values <- new_keys[match(iupacs, keys)]
  names(values) <- names(x)
  new_glycan_structure(values, graphs)
}

.compact_table_record <- function(nodes, edges, parts, subs, anomer, alditol) {
  n <- nrow(nodes)
  parts <- lapply(seq_len(nrow(parts)), function(i) {
    p <- list(
      root = parts$root_node[[i]],
      linkage = parts$linkage[[i]],
      parents = parts$parents[[i]]
    )
    if ("nodes" %in% names(parts)) {
      p$nodes <- parts$nodes[[i]]
    }
    p <- normalize_floating_part(p, n)
    if (is.null(p$nodes)) {
      seen <- pending <- p$root
      while (length(pending)) {
        next_nodes <- edges$to_node[edges$from_node %in% pending]
        pending <- setdiff(next_nodes, seen)
        seen <- c(seen, pending)
      }
      p$nodes <- sort(as.integer(seen))
    }
    p
  })
  subs <- lapply(seq_len(nrow(subs)), function(i) {
    normalize_floating_substituent(
      list(substituent = subs$substituent[[i]], parents = subs$parents[[i]]),
      n
    )
  })
  links <- c(edges$linkage, vapply(parts, `[[`, character(1), "linkage"))
  if (any(!grepl("^[ab?][12?]-([1-9](/[1-9])*|[?])$", links))) {
    stop("unsupported native linkage")
  }
  .validate_structure_arrays(list(
    mono = nodes$mono,
    sub = nodes$sub,
    edges = as.integer(rbind(edges$from_node, edges$to_node)),
    linkage = edges$linkage,
    anomer = anomer,
    alditol = alditol,
    floating_parts = parts,
    floating_substituents = subs
  ))
}

.compact_table_structure <- function(records, input_names) {
  native <- .compact_arrays_native(
    records,
    base::order,
    .compact_bliss_labels,
    identical(Sys.getlocale("LC_COLLATE"), "C")
  )
  ok <- vapply(native, function(a) a$status %in% c("ok", "missing"), logical(1))
  if (!all(ok)) {
    return(NULL)
  }
  keys <- vapply(
    native,
    function(a) if (a$status == "missing") NA_character_ else a$iupac,
    character(1)
  )
  keep <- !is.na(keys) & !duplicated(keys)
  graphs <- lapply(native[keep], .compact_build_graph)
  names(graphs) <- keys[keep]
  names(keys) <- input_names
  new_glycan_structure(keys, graphs)
}

.compact_enumerate_localizations <- function(
  graph,
  domains,
  combinations,
  np,
  ns,
  input_id
) {
  if (combinations > .Machine$integer.max) {
    return(NULL)
  }
  record <- tryCatch(.compact_graph_record(graph), error = function(e) NULL)
  if (is.null(record)) {
    return(NULL)
  }
  result <- .compact_localizations_native(
    record,
    domains,
    as.integer(combinations)
  )
  if (!.compact_graph_ok(result)) {
    return(NULL)
  }
  graphs <- lapply(result$records, function(a) {
    out <- graph
    ne <- igraph::ecount(graph)
    if (length(a$linkage) > ne) {
      added <- seq.int(ne + 1L, length(a$linkage))
      out <- igraph::add_edges(
        out,
        a$edges[as.vector(rbind(2L * added - 1L, 2L * added))],
        linkage = a$linkage[added]
      )
    }
    out <- igraph::set_vertex_attr(out, "sub", value = a$sub)
    out <- delete_floating_parts_attr(out)
    delete_floating_substituents_attr(out)
  })
  assignments <- lapply(result$assignments, function(parents) {
    tibble::tibble(
      glycan_id = rep(as.integer(input_id), np + ns),
      part_id = c(seq_len(np), rep(NA_integer_, ns)),
      parent_node = parents,
      substituent_id = c(rep(NA_integer_, np), seq_len(ns))
    )
  })
  list(graphs = graphs, assignments = assignments)
}

.compact_localize_parts <- function(graph, assignments) {
  record <- tryCatch(
    .compact_graph_record(graph, validate = TRUE),
    error = function(e) NULL
  )
  if (is.null(record)) {
    return(NULL)
  }
  a <- .compact_localize_parts_native(
    record,
    assignments$part_id,
    assignments$parent_node
  )
  if (!.compact_graph_ok(a)) {
    return(NULL)
  }
  ne <- igraph::ecount(graph)
  added <- seq.int(ne + 1L, length(a$linkage))
  graph <- igraph::add_edges(
    graph,
    a$edges[as.vector(rbind(2L * added - 1L, 2L * added))],
    linkage = a$linkage[added]
  )
  set_floating_parts_attr(graph, a$floating_parts)
}
