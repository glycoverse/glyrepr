#' Convert Glycan Structures to Graph Tables
#'
#' @description
#' `structure_nodes()`, `structure_edges()`, `structure_floating_parts()`, and
#' `structure_floating_substituents()` convert
#' a glycan structure vector or one glycan `igraph` to normalized graph tables.
#' `structure_from_tibbles()` rebuilds a `glyrepr_structure` vector from those
#' tables, a vector of reducing-end anomers, and optional alditol status.
#'
#' The `glycan_id` column is the integer position of each glycan in the input
#' vector. Duplicate structures are expanded to their original vector positions.
#' Missing structures have no node or edge rows and are reconstructed from
#' missing values in `anomers`.
#' If `x` is named, all four tibbles also contain a `glycan_name` column.
#' `structure_from_tibbles()` uses `glycan_name` as output names when that
#' column is present.
#'
#' For structure-vector input, `structure_nodes()$node_id` follows residue order
#' in the complete canonical IUPAC-condensed string. For graph input, it follows
#' the graph's current numeric vertex positions without canonicalizing or
#' renumbering them. A graph is represented with `glycan_id = 1L` and no
#' `glycan_name` column.
#'
#' In `structure_floating_parts()`, `root_node` and every integer in the `nodes`
#' and `parents` list-columns refer to `structure_nodes()$node_id` for the same
#' glycan. `nodes` contains every node in the floating component. An empty
#' `parents` vector means all feasible nodes outside that component are
#' candidates. The `linkage` column describes the virtual attachment to an
#' unresolved parent; this attachment is intentionally absent from
#' `structure_edges()`. During
#' reconstruction, a row with exactly one effective candidate parent is
#' normalized to an ordinary edge and is therefore absent from the resulting
#' `structure_floating_parts()` table.
#'
#' Parent indices written after `|` in an IUPAC-condensed floating part are
#' complete-sequence node IDs, identical to `structure_nodes()$node_id` for a
#' canonical structure. Residues in floating blocks precede the main tree, and
#' substituent blocks contribute no nodes. `structure_from_tibbles()` expects
#' these same global node IDs and preserves cross-component domains.
#'
#' In `structure_floating_substituents()`, each row describes one unresolved
#' substituent. `substituent` is its canonical position-and-name token, and the
#' `parents` list-column contains candidate global node IDs. An empty vector
#' means all feasible residue nodes are candidates. A singleton candidate is
#' normalized into `structure_nodes()$sub`, so it does not remain in the
#' floating-substituent table.
#'
#' @param x A glycan structure vector or one glycan `igraph`.
#' @param nodes A data frame with columns `glycan_id`, `node_id`, `mono`, and
#'   `sub`, and optionally `glycan_name`.
#' @param edges A data frame with columns `glycan_id`, `edge_id`, `from_node`,
#'   `to_node`, and `linkage`, and optionally `glycan_name`.
#' @param anomers A character vector of reducing-end anomers, one per glycan.
#' @param floating_parts A data frame returned by
#'   `structure_floating_parts()`, or `NULL` when no floating parts are present.
#' @param floating_substituents A data frame returned by
#'   `structure_floating_substituents()`, or `NULL` when no floating
#'   substituents are present.
#' @param alditols A logical vector indicating alditol status, either one value
#'   or one per glycan. Missing values are allowed only for missing glycans.
#'   Defaults to `FALSE` for backward compatibility.
#'
#' @returns
#' - `structure_nodes()` returns a tibble with columns `glycan_id`, `node_id`,
#'   `mono`, and `sub`.
#' - `structure_edges()` returns a tibble with columns `glycan_id`, `edge_id`,
#'   `from_node`, `to_node`, and `linkage`.
#' - `structure_floating_parts()` returns a tibble with columns `glycan_id`,
#'   `part_id`, `root_node`, the list-column `nodes`, `linkage`, and the
#'   list-column `parents`.
#' - `structure_floating_substituents()` returns a tibble with columns
#'   `glycan_id`, `substituent_id`, `substituent`, and the list-column
#'   `parents`.
#' - `structure_from_tibbles()` returns a `glyrepr_structure` vector.
#'
#' @examples
#' glycans <- c(o_glycan_core_1(), o_glycan_core_1())
#' nodes <- structure_nodes(glycans)
#' edges <- structure_edges(glycans)
#' structure_from_tibbles(nodes, edges, get_anomer(glycans))
#'
#' floating <- as_glycan_structure(
#'   "{6S|1,2}{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' floating_parts <- structure_floating_parts(floating)
#' floating_substituents <- structure_floating_substituents(floating)
#' structure_from_tibbles(
#'   structure_nodes(floating),
#'   structure_edges(floating),
#'   get_anomer(floating),
#'   floating_parts,
#'   floating_substituents
#' )
#'
#' @name structure_tables
NULL

structure_table_input <- function(x) {
  if (inherits(x, "igraph")) {
    return(list(graphs = list(x), glycan_names = NULL))
  }

  checkmate::assert_class(x, "glyrepr_structure")
  list(graphs = as.list(x), glycan_names = names(x))
}


#' @rdname structure_tables
#' @export
structure_nodes <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_nodes(has_glycan_names))
  }

  node_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_nodes_one(glycan_id, graph, glycan_name)
    }
  )

  dplyr::bind_rows(node_tables)
}


#' @rdname structure_tables
#' @export
structure_edges <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_edges(has_glycan_names))
  }

  edge_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_edges_one(glycan_id, graph, glycan_name)
    }
  )

  dplyr::bind_rows(edge_tables)
}


#' @rdname structure_tables
#' @export
structure_floating_parts <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_floating_parts(has_glycan_names))
  }

  part_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_floating_parts_one(glycan_id, graph, glycan_name)
    }
  )

  dplyr::bind_rows(part_tables)
}


#' @rdname structure_tables
#' @export
structure_floating_substituents <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_floating_substituents(has_glycan_names))
  }

  substituent_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_floating_substituents_one(
        glycan_id,
        graph,
        glycan_name
      )
    }
  )

  dplyr::bind_rows(substituent_tables)
}


#' List Candidate Parents for Floating Metadata
#'
#' @description
#' `structure_floating_candidates()` expands floating-part and
#' floating-substituent metadata to one row per candidate parent. This provides
#' a uniform representation for explicitly restricted and unrestricted
#' floating metadata.
#'
#' For an unrestricted floating part, every feasible node outside its own
#' component is returned; for an unrestricted floating substituent, every
#' feasible residue node in the complete structure is returned. These rows use
#' `scope = "all"`. For metadata written with an explicit `|<parents>` suffix,
#' only the declared parent nodes are returned and `scope` is `"explicit"`.
#'
#' Floating-part rows have a non-missing `part_id`, `root_node`, and `linkage`.
#' Floating-substituent rows instead have a non-missing `substituent_id` and
#' `substituent`. Exactly one of `part_id` and `substituent_id` is non-missing
#' in each row.
#'
#' Node indices refer to `structure_nodes()$node_id` for the same glycan. Missing
#' structures and structures without floating metadata contribute no rows.
#' Duplicate structures are expanded to their original vector positions. For
#' graph input, node indices are current numeric vertex positions and
#' `glycan_id` is `1L`. If vector input is named, the result also contains a
#' `glycan_name` column.
#'
#' @param x A glycan structure vector or one glycan `igraph`.
#'
#' @returns A tibble with columns `glycan_id`, `part_id`, `root_node`,
#'   `parent_node`, `linkage`, `scope`, `substituent_id`, and `substituent`,
#'   plus `glycan_name` when `x` is named.
#'
#' @examples
#' glycans <- as_glycan_structure(c(
#'   unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
#'   restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
#'   substituent = "{6S}Gal(a1-3)Gal(a1-"
#' ))
#' structure_floating_candidates(glycans)
#'
#' @export
structure_floating_candidates <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_floating_candidates(has_glycan_names))
  }

  candidate_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_floating_candidates_one(glycan_id, graph, glycan_name)
    }
  )

  dplyr::bind_rows(candidate_tables)
}


#' Identify Main and Floating Structure Components
#'
#' @description
#' `structure_component_membership()` identifies whether each glycan node
#' belongs to the main tree or to a floating part. It provides a public,
#' normalized alternative to reading the private floating-part graph
#' attribute.
#'
#' Main-tree nodes have `component_type = "main"` and a missing `part_id`.
#' Nodes in a floating subtree have `component_type = "floating"` and the
#' corresponding `part_id` from [structure_floating_parts()].
#' Floating substituents do not introduce vertices, so they do not create
#' additional component-membership rows.
#'
#' Node indices refer to `structure_nodes()$node_id` for the same glycan.
#' Missing structures contribute no rows. Duplicate structures are expanded to
#' their original vector positions. For graph input, node indices are current
#' numeric vertex positions and `glycan_id` is `1L`. If vector input is named,
#' the result also contains a `glycan_name` column.
#'
#' @param x A glycan structure vector or one glycan `igraph`.
#'
#' @returns A tibble with columns `glycan_id`, `node_id`, `component_type`,
#'   and `part_id`, plus `glycan_name` when `x` is named.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' structure_component_membership(glycan)
#'
#' @export
structure_component_membership <- function(x) {
  input <- structure_table_input(x)
  graphs <- input$graphs
  glycan_names <- input$glycan_names
  has_glycan_names <- !is.null(glycan_names)
  if (length(graphs) == 0) {
    return(empty_structure_component_membership(has_glycan_names))
  }

  membership_tables <- purrr::map2(
    seq_along(graphs),
    graphs,
    function(glycan_id, graph) {
      glycan_name <- if (has_glycan_names) {
        glycan_names[[glycan_id]]
      } else {
        NULL
      }
      structure_component_membership_one(glycan_id, graph, glycan_name)
    }
  )

  dplyr::bind_rows(membership_tables)
}


#' List Potential Virtual Edges for Floating Parts
#'
#' @description
#' `structure_candidate_edges()` represents every potential floating-part
#' attachment as an explicit virtual edge. `from_node` is a candidate parent in
#' another floating component or the main tree, and `to_node` is the root of
#' the floating part.
#'
#' The rows correspond one-to-one with the floating-part rows from
#' [structure_floating_candidates()]. Floating substituents do not create
#' virtual graph edges. For unrestricted `{<floating>}` parts, every feasible
#' node outside the part's own component is returned and `scope` is `"all"`.
#' For explicitly restricted parts, only the declared parent nodes are returned
#' and `scope` is `"explicit"`.
#'
#' Node indices refer to `structure_nodes()$node_id` for the same glycan.
#' Missing structures and structures without floating parts contribute no
#' rows. Duplicate structures are expanded to their original vector positions.
#' For graph input, node indices are current numeric vertex positions and
#' `glycan_id` is `1L`. If vector input is named, the result also contains a
#' `glycan_name` column.
#'
#' @param x A glycan structure vector or one glycan `igraph`.
#'
#' @returns A tibble with columns `glycan_id`, `part_id`, `from_node`,
#'   `to_node`, `linkage`, and `scope`, plus `glycan_name` when `x` is named.
#'
#' @examples
#' glycan <- as_glycan_structure(
#'   "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
#' )
#' structure_candidate_edges(glycan)
#'
#' @export
structure_candidate_edges <- function(x) {
  candidates <- structure_floating_candidates(x)
  candidates <- candidates[!is.na(candidates$part_id), , drop = FALSE]

  candidates <- dplyr::rename(
    candidates,
    from_node = "parent_node",
    to_node = "root_node"
  )
  columns <- c(
    "glycan_id",
    if ("glycan_name" %in% names(candidates)) "glycan_name",
    "part_id",
    "from_node",
    "to_node",
    "linkage",
    "scope"
  )

  candidates[, columns]
}


#' @rdname structure_tables
#' @export
structure_from_tibbles <- function(
  nodes,
  edges,
  anomers,
  floating_parts = NULL,
  floating_substituents = NULL,
  alditols = FALSE
) {
  nodes <- validate_structure_nodes_table(nodes)
  edges <- validate_structure_edges_table(edges)
  anomers <- validate_structure_anomers(anomers)
  floating_parts <- validate_structure_floating_parts_table(floating_parts)
  floating_substituents <- validate_structure_floating_substituents_table(
    floating_substituents
  )
  alditols <- validate_structure_alditols(alditols, anomers)

  validate_structure_table_glycan_ids(
    nodes$glycan_id,
    length(anomers),
    "nodes"
  )
  validate_structure_table_glycan_ids(
    edges$glycan_id,
    length(anomers),
    "edges"
  )
  validate_structure_table_glycan_ids(
    floating_parts$glycan_id,
    length(anomers),
    "floating_parts"
  )
  validate_structure_table_glycan_ids(
    floating_substituents$glycan_id,
    length(anomers),
    "floating_substituents"
  )
  glycan_names <- structure_table_glycan_names(
    nodes,
    edges,
    floating_parts,
    floating_substituents,
    anomers
  )

  if (length(anomers) == 0) {
    return(glycan_structure())
  }

  graphs <- purrr::map(seq_along(anomers), function(glycan_id) {
    build_structure_graph_from_table_rows(
      nodes[nodes$glycan_id == glycan_id, , drop = FALSE],
      edges[edges$glycan_id == glycan_id, , drop = FALSE],
      floating_parts[
        floating_parts$glycan_id == glycan_id,
        ,
        drop = FALSE
      ],
      floating_substituents[
        floating_substituents$glycan_id == glycan_id,
        ,
        drop = FALSE
      ],
      anomers[[glycan_id]],
      alditols[[glycan_id]],
      glycan_id
    )
  })

  out <- do.call(glycan_structure, graphs)
  names(out) <- glycan_names
  out
}


#' Create an empty structure node tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the `structure_nodes()` columns.
#' @noRd
empty_structure_nodes <- function(has_glycan_name = FALSE) {
  out <- tibble::tibble(
    glycan_id = integer(),
    node_id = integer(),
    mono = character(),
    sub = character()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Create an empty structure edge tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the `structure_edges()` columns.
#' @noRd
empty_structure_edges <- function(has_glycan_name = FALSE) {
  out <- tibble::tibble(
    glycan_id = integer(),
    edge_id = integer(),
    from_node = integer(),
    to_node = integer(),
    linkage = character()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Create an empty floating-part tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the `structure_floating_parts()` columns.
#' @noRd
empty_structure_floating_parts <- function(has_glycan_name = FALSE) {
  out <- tibble::tibble(
    glycan_id = integer(),
    part_id = integer(),
    root_node = integer(),
    nodes = list(),
    linkage = character(),
    parents = list()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Create an empty floating-substituent tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the
#'   `structure_floating_substituents()` columns.
#' @noRd
empty_structure_floating_substituents <- function(
  has_glycan_name = FALSE
) {
  out <- tibble::tibble(
    glycan_id = integer(),
    substituent_id = integer(),
    substituent = character(),
    parents = list()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Create an empty floating-candidate tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the
#'   `structure_floating_candidates()` columns.
#' @noRd
empty_structure_floating_candidates <- function(
  has_glycan_name = FALSE
) {
  out <- tibble::tibble(
    glycan_id = integer(),
    part_id = integer(),
    root_node = integer(),
    parent_node = integer(),
    linkage = character(),
    scope = character(),
    substituent_id = integer(),
    substituent = character()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Create an empty component-membership tibble
#'
#' @param has_glycan_name Whether to include a `glycan_name` column.
#' @returns A zero-row tibble with the
#'   `structure_component_membership()` columns.
#' @noRd
empty_structure_component_membership <- function(
  has_glycan_name = FALSE
) {
  out <- tibble::tibble(
    glycan_id = integer(),
    node_id = integer(),
    component_type = character(),
    part_id = integer()
  )

  if (has_glycan_name) {
    out <- tibble::add_column(
      out,
      glycan_name = character(),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a structure node tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with node rows for one glycan.
#' @noRd
structure_nodes_one <- function(glycan_id, graph, glycan_name = NULL) {
  if (is.null(graph)) {
    return(empty_structure_nodes(!is.null(glycan_name)))
  }

  node_count <- igraph::vcount(graph)
  out <- tibble::tibble(
    glycan_id = rep(as.integer(glycan_id), node_count),
    node_id = seq_len(node_count),
    mono = igraph::vertex_attr(graph, "mono"),
    sub = igraph::vertex_attr(graph, "sub")
  )

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, node_count),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a structure edge tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with edge rows for one glycan.
#' @noRd
structure_edges_one <- function(glycan_id, graph, glycan_name = NULL) {
  if (is.null(graph) || igraph::ecount(graph) == 0) {
    return(empty_structure_edges(!is.null(glycan_name)))
  }

  edge_count <- igraph::ecount(graph)
  edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  out <- tibble::tibble(
    glycan_id = rep(as.integer(glycan_id), edge_count),
    edge_id = seq_len(edge_count),
    from_node = as.integer(edge_ends[, 1]),
    to_node = as.integer(edge_ends[, 2]),
    linkage = igraph::edge_attr(graph, "linkage")
  )

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, edge_count),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a floating-part tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with floating-part rows for one glycan.
#' @noRd
structure_floating_parts_one <- function(
  glycan_id,
  graph,
  glycan_name = NULL
) {
  if (is.null(graph)) {
    return(empty_structure_floating_parts(!is.null(glycan_name)))
  }

  parts <- normalize_floating_parts(graph)
  if (length(parts) == 0) {
    return(empty_structure_floating_parts(!is.null(glycan_name)))
  }

  out <- tibble::tibble(
    glycan_id = rep(as.integer(glycan_id), length(parts)),
    part_id = seq_along(parts),
    root_node = purrr::map_int(parts, "root"),
    nodes = purrr::map(parts, "nodes"),
    linkage = purrr::map_chr(parts, "linkage"),
    parents = purrr::map(parts, "parents")
  )

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, length(parts)),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a floating-substituent tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with floating-substituent rows for one glycan.
#' @noRd
structure_floating_substituents_one <- function(
  glycan_id,
  graph,
  glycan_name = NULL
) {
  if (is.null(graph)) {
    return(empty_structure_floating_substituents(!is.null(glycan_name)))
  }

  substituents <- normalize_floating_substituents(graph)
  if (length(substituents) == 0) {
    return(empty_structure_floating_substituents(!is.null(glycan_name)))
  }

  out <- tibble::tibble(
    glycan_id = rep(as.integer(glycan_id), length(substituents)),
    substituent_id = seq_along(substituents),
    substituent = purrr::map_chr(substituents, "substituent"),
    parents = purrr::map(substituents, "parents")
  )

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, length(substituents)),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a floating-candidate tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with candidate attachment rows for one glycan.
#' @noRd
structure_floating_candidates_one <- function(
  glycan_id,
  graph,
  glycan_name = NULL
) {
  if (is.null(graph)) {
    return(empty_structure_floating_candidates(!is.null(glycan_name)))
  }

  parts <- normalize_floating_parts(graph)
  substituents <- normalize_floating_substituents(graph)
  if (length(parts) == 0 && length(substituents) == 0) {
    return(empty_structure_floating_candidates(!is.null(glycan_name)))
  }

  part_tables <- purrr::map2(
    parts,
    seq_along(parts),
    function(part, part_id) {
      unrestricted <- length(part$parents) == 0
      parent_nodes <- if (unrestricted) {
        floating_part_candidate_parents(graph, part)
      } else {
        part$parents
      }

      tibble::tibble(
        glycan_id = rep(as.integer(glycan_id), length(parent_nodes)),
        part_id = rep(as.integer(part_id), length(parent_nodes)),
        root_node = rep(as.integer(part$root), length(parent_nodes)),
        parent_node = as.integer(parent_nodes),
        linkage = rep(part$linkage, length(parent_nodes)),
        scope = rep(
          if (unrestricted) "all" else "explicit",
          length(parent_nodes)
        ),
        substituent_id = rep(NA_integer_, length(parent_nodes)),
        substituent = rep(NA_character_, length(parent_nodes))
      )
    }
  )
  substituent_tables <- purrr::map2(
    substituents,
    seq_along(substituents),
    function(substituent, substituent_id) {
      unrestricted <- length(substituent$parents) == 0
      parent_nodes <- if (unrestricted) {
        floating_substituent_candidate_parents(graph, substituent)
      } else {
        substituent$parents
      }

      tibble::tibble(
        glycan_id = rep(as.integer(glycan_id), length(parent_nodes)),
        part_id = rep(NA_integer_, length(parent_nodes)),
        root_node = rep(NA_integer_, length(parent_nodes)),
        parent_node = as.integer(parent_nodes),
        linkage = rep(NA_character_, length(parent_nodes)),
        scope = rep(
          if (unrestricted) "all" else "explicit",
          length(parent_nodes)
        ),
        substituent_id = rep(
          as.integer(substituent_id),
          length(parent_nodes)
        ),
        substituent = rep(
          substituent$substituent,
          length(parent_nodes)
        )
      )
    }
  )
  out <- dplyr::bind_rows(c(part_tables, substituent_tables))

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, nrow(out)),
      .after = "glycan_id"
    )
  }

  out
}


#' Convert one graph to a component-membership tibble
#'
#' @param glycan_id Integer position of the glycan.
#' @param graph An igraph object or `NULL` for a missing structure.
#' @param glycan_name Optional glycan name.
#' @returns A tibble with component-membership rows for one glycan.
#' @noRd
structure_component_membership_one <- function(
  glycan_id,
  graph,
  glycan_name = NULL
) {
  if (is.null(graph)) {
    return(empty_structure_component_membership(!is.null(glycan_name)))
  }

  node_count <- igraph::vcount(graph)
  parts <- normalize_floating_parts(graph)
  if (length(parts) == 0) {
    component_type <- rep("main", node_count)
    part_id <- rep(NA_integer_, node_count)
  } else {
    part_id <- rep(NA_integer_, node_count)
    for (floating_id in seq_along(parts)) {
      part_id[parts[[floating_id]]$nodes] <- floating_id
    }
    component_type <- ifelse(is.na(part_id), "main", "floating")
  }

  out <- tibble::tibble(
    glycan_id = rep(as.integer(glycan_id), node_count),
    node_id = seq_len(node_count),
    component_type = component_type,
    part_id = as.integer(part_id)
  )

  if (!is.null(glycan_name)) {
    out <- tibble::add_column(
      out,
      glycan_name = rep(glycan_name, node_count),
      .after = "glycan_id"
    )
  }

  out
}


#' Validate a structure node table
#'
#' @param nodes A candidate node table.
#' @returns A tibble with the required node columns.
#' @noRd
validate_structure_nodes_table <- function(nodes) {
  required_cols <- c("glycan_id", "node_id", "mono", "sub")
  nodes <- validate_structure_table(
    nodes,
    required_cols,
    "nodes",
    optional_cols = "glycan_name"
  )

  validate_integerish_structure_column(nodes, "glycan_id", "nodes")
  validate_integerish_structure_column(nodes, "node_id", "nodes")
  validate_character_structure_column(nodes, "mono", "nodes")
  validate_character_structure_column(nodes, "sub", "nodes")
  if ("glycan_name" %in% names(nodes)) {
    validate_character_structure_column(
      nodes,
      "glycan_name",
      "nodes",
      any.missing = TRUE
    )
  }

  nodes$glycan_id <- as.integer(nodes$glycan_id)
  nodes$node_id <- as.integer(nodes$node_id)
  nodes
}


#' Validate a structure edge table
#'
#' @param edges A candidate edge table.
#' @returns A tibble with the required edge columns.
#' @noRd
validate_structure_edges_table <- function(edges) {
  required_cols <- c("glycan_id", "edge_id", "from_node", "to_node", "linkage")
  edges <- validate_structure_table(
    edges,
    required_cols,
    "edges",
    optional_cols = "glycan_name"
  )

  validate_integerish_structure_column(edges, "glycan_id", "edges")
  validate_integerish_structure_column(edges, "edge_id", "edges")
  validate_integerish_structure_column(edges, "from_node", "edges")
  validate_integerish_structure_column(edges, "to_node", "edges")
  validate_character_structure_column(edges, "linkage", "edges")
  if ("glycan_name" %in% names(edges)) {
    validate_character_structure_column(
      edges,
      "glycan_name",
      "edges",
      any.missing = TRUE
    )
  }

  edges$glycan_id <- as.integer(edges$glycan_id)
  edges$edge_id <- as.integer(edges$edge_id)
  edges$from_node <- as.integer(edges$from_node)
  edges$to_node <- as.integer(edges$to_node)
  edges
}


#' Validate a floating-part table
#'
#' @param floating_parts A candidate floating-part table or `NULL`.
#' @returns A tibble with the required floating-part columns.
#' @noRd
validate_structure_floating_parts_table <- function(floating_parts) {
  if (is.null(floating_parts)) {
    return(empty_structure_floating_parts())
  }

  required_cols <- c(
    "glycan_id",
    "part_id",
    "root_node",
    "linkage",
    "parents"
  )
  floating_parts <- validate_structure_table(
    floating_parts,
    required_cols,
    "floating_parts",
    optional_cols = c("glycan_name", "nodes")
  )

  validate_integerish_structure_column(
    floating_parts,
    "glycan_id",
    "floating_parts"
  )
  validate_integerish_structure_column(
    floating_parts,
    "part_id",
    "floating_parts"
  )
  validate_integerish_structure_column(
    floating_parts,
    "root_node",
    "floating_parts"
  )
  validate_character_structure_column(
    floating_parts,
    "linkage",
    "floating_parts"
  )
  if ("nodes" %in% names(floating_parts)) {
    if (!is.list(floating_parts$nodes)) {
      cli::cli_abort(
        "{.arg floating_parts} column {.field nodes} must be a list-column."
      )
    }
    valid_nodes <- purrr::map_lgl(
      floating_parts$nodes,
      function(nodes) {
        length(nodes) > 0 &&
          checkmate::test_integerish(
            nodes,
            lower = 1,
            any.missing = FALSE
          )
      }
    )
    if (!all(valid_nodes)) {
      cli::cli_abort(
        "{.arg floating_parts} column {.field nodes} must contain non-empty positive integer vectors."
      )
    }
  }
  if (!is.list(floating_parts$parents)) {
    cli::cli_abort(
      "{.arg floating_parts} column {.field parents} must be a list-column."
    )
  }
  valid_parents <- purrr::map_lgl(
    floating_parts$parents,
    ~ checkmate::test_integerish(
      .x,
      lower = 1,
      any.missing = FALSE
    )
  )
  if (!all(valid_parents)) {
    cli::cli_abort(
      "{.arg floating_parts} column {.field parents} must contain positive integer vectors."
    )
  }
  if ("glycan_name" %in% names(floating_parts)) {
    validate_character_structure_column(
      floating_parts,
      "glycan_name",
      "floating_parts",
      any.missing = TRUE
    )
  }

  floating_parts$glycan_id <- as.integer(floating_parts$glycan_id)
  floating_parts$part_id <- as.integer(floating_parts$part_id)
  floating_parts$root_node <- as.integer(floating_parts$root_node)
  if ("nodes" %in% names(floating_parts)) {
    floating_parts$nodes <- purrr::map(
      floating_parts$nodes,
      as.integer
    )
  }
  floating_parts$parents <- purrr::map(
    floating_parts$parents,
    as.integer
  )
  floating_parts
}


#' Validate a floating-substituent table
#'
#' @param floating_substituents A candidate floating-substituent table or
#'   `NULL`.
#' @returns A tibble with the required floating-substituent columns.
#' @noRd
validate_structure_floating_substituents_table <- function(
  floating_substituents
) {
  if (is.null(floating_substituents)) {
    return(empty_structure_floating_substituents())
  }

  required_cols <- c(
    "glycan_id",
    "substituent_id",
    "substituent",
    "parents"
  )
  floating_substituents <- validate_structure_table(
    floating_substituents,
    required_cols,
    "floating_substituents",
    optional_cols = "glycan_name"
  )

  validate_integerish_structure_column(
    floating_substituents,
    "glycan_id",
    "floating_substituents"
  )
  validate_integerish_structure_column(
    floating_substituents,
    "substituent_id",
    "floating_substituents"
  )
  validate_character_structure_column(
    floating_substituents,
    "substituent",
    "floating_substituents"
  )
  if (!is.list(floating_substituents$parents)) {
    cli::cli_abort(
      "{.arg floating_substituents} column {.field parents} must be a list-column."
    )
  }
  valid_parents <- purrr::map_lgl(
    floating_substituents$parents,
    ~ checkmate::test_integerish(
      .x,
      lower = 1,
      any.missing = FALSE
    )
  )
  if (!all(valid_parents)) {
    cli::cli_abort(
      "{.arg floating_substituents} column {.field parents} must contain positive integer vectors."
    )
  }
  if ("glycan_name" %in% names(floating_substituents)) {
    validate_character_structure_column(
      floating_substituents,
      "glycan_name",
      "floating_substituents",
      any.missing = TRUE
    )
  }

  floating_substituents$glycan_id <- as.integer(
    floating_substituents$glycan_id
  )
  floating_substituents$substituent_id <- as.integer(
    floating_substituents$substituent_id
  )
  floating_substituents$parents <- purrr::map(
    floating_substituents$parents,
    as.integer
  )
  floating_substituents
}


#' Validate a graph-table input
#'
#' @param table A candidate graph table.
#' @param required_cols Required column names.
#' @param arg_name User-facing argument name.
#' @param optional_cols Optional column names to keep.
#' @returns A tibble with required columns.
#' @noRd
validate_structure_table <- function(
  table,
  required_cols,
  arg_name,
  optional_cols = character()
) {
  if (!is.data.frame(table)) {
    cli::cli_abort("{.arg {arg_name}} must be a data frame.")
  }

  missing_cols <- setdiff(required_cols, names(table))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg {arg_name}} must contain the required columns.",
      "x" = "Missing column{?s}: {.field {missing_cols}}."
    ))
  }

  selected_cols <- c(
    required_cols[[1]],
    intersect(optional_cols, names(table)),
    required_cols[-1]
  )

  tibble::as_tibble(table[selected_cols])
}


#' Validate an integer-like graph-table column
#'
#' @param table A graph table.
#' @param column Column name.
#' @param arg_name User-facing argument name.
#' @returns Invisible `NULL`.
#' @noRd
validate_integerish_structure_column <- function(table, column, arg_name) {
  if (
    !checkmate::test_integerish(
      table[[column]],
      lower = 1,
      any.missing = FALSE
    )
  ) {
    cli::cli_abort(
      "{.arg {arg_name}} column {.field {column}} must contain positive integer values."
    )
  }

  invisible(NULL)
}


#' Validate a character graph-table column
#'
#' @param table A graph table.
#' @param column Column name.
#' @param arg_name User-facing argument name.
#' @param any.missing Whether missing values are allowed.
#' @returns Invisible `NULL`.
#' @noRd
validate_character_structure_column <- function(
  table,
  column,
  arg_name,
  any.missing = FALSE
) {
  if (!checkmate::test_character(table[[column]], any.missing = any.missing)) {
    cli::cli_abort(
      "{.arg {arg_name}} column {.field {column}} must contain character values."
    )
  }

  invisible(NULL)
}


#' Derive output names from graph tables or anomers
#'
#' @param nodes A validated node table.
#' @param edges A validated edge table.
#' @param floating_parts A validated floating-part table.
#' @param floating_substituents A validated floating-substituent table.
#' @param anomers A validated anomer vector.
#' @returns A character vector of names, or `NULL`.
#' @noRd
structure_table_glycan_names <- function(
  nodes,
  edges,
  floating_parts,
  floating_substituents,
  anomers
) {
  table_names <- structure_table_glycan_names_from_rows(
    nodes,
    edges,
    floating_parts,
    floating_substituents,
    anomers
  )
  if (is.null(table_names)) {
    return(names(anomers))
  }

  out_names <- names(anomers)
  if (is.null(out_names)) {
    out_names <- rep("", length(anomers))
  }

  has_table_name <- !is.na(table_names)
  out_names[has_table_name] <- table_names[has_table_name]
  out_names
}


#' Derive glycan names from node and edge rows
#'
#' @param nodes A validated node table.
#' @param edges A validated edge table.
#' @param floating_parts A validated floating-part table.
#' @param floating_substituents A validated floating-substituent table.
#' @param anomers A validated anomer vector.
#' @returns A character vector with `NA` where no table name is available, or
#'   `NULL` when neither table contains `glycan_name`.
#' @noRd
structure_table_glycan_names_from_rows <- function(
  nodes,
  edges,
  floating_parts,
  floating_substituents,
  anomers
) {
  name_rows <- dplyr::bind_rows(
    structure_table_name_rows(nodes),
    structure_table_name_rows(edges),
    structure_table_name_rows(floating_parts),
    structure_table_name_rows(floating_substituents)
  )

  if (nrow(name_rows) == 0) {
    return(NULL)
  }

  table_names <- rep(NA_character_, length(anomers))
  for (glycan_id in unique(name_rows$glycan_id)) {
    glycan_names <- unique(name_rows$glycan_name[
      name_rows$glycan_id == glycan_id
    ])
    glycan_names <- glycan_names[!is.na(glycan_names)]

    if (length(glycan_names) > 1) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has multiple {.field glycan_name} values."
      )
    }

    if (length(glycan_names) == 1) {
      table_names[[glycan_id]] <- glycan_names[[1]]
    }
  }

  table_names
}


#' Extract glycan-name rows from a graph table
#'
#' @param table A validated graph table.
#' @returns A tibble with `glycan_id` and `glycan_name`.
#' @noRd
structure_table_name_rows <- function(table) {
  if (!"glycan_name" %in% names(table)) {
    return(tibble::tibble(
      glycan_id = integer(),
      glycan_name = character()
    ))
  }

  unique(table[c("glycan_id", "glycan_name")])
}


#' Validate anomers for graph-table reconstruction
#'
#' @param anomers Candidate reducing-end anomer values.
#' @returns The validated anomer vector.
#' @noRd
validate_structure_anomers <- function(anomers) {
  checkmate::assert_character(anomers, any.missing = TRUE)

  valid <- is.na(anomers) | valid_anomer(anomers)
  if (!all(valid)) {
    invalid_anomers <- unique(anomers[!valid])
    cli::cli_abort(
      "Invalid anomer{?s}: {.val {invalid_anomers}}."
    )
  }

  anomers
}


#' Validate alditol status for graph-table reconstruction
#'
#' @param alditols Candidate alditol values.
#' @param anomers Validated reducing-end anomers.
#' @returns A logical vector with one value per anomer.
#' @noRd
validate_structure_alditols <- function(alditols, anomers) {
  checkmate::assert_logical(alditols, any.missing = TRUE)

  n_glycans <- length(anomers)
  if (!length(alditols) %in% c(1L, n_glycans)) {
    cli::cli_abort(
      "{.arg alditols} must have length one or the same length as {.arg anomers}."
    )
  }

  alditols <- vctrs::vec_recycle(alditols, n_glycans)
  invalid_missing <- !is.na(anomers) & is.na(alditols)
  if (any(invalid_missing)) {
    cli::cli_abort(
      "{.arg alditols} cannot be missing for a non-missing glycan."
    )
  }

  alditols
}


#' Validate graph-table glycan IDs against anomer length
#'
#' @param glycan_id Integer glycan IDs from a graph table.
#' @param n_anomers Length of the anomer vector.
#' @param arg_name User-facing argument name.
#' @returns Invisible `NULL`.
#' @noRd
validate_structure_table_glycan_ids <- function(
  glycan_id,
  n_anomers,
  arg_name
) {
  if (length(glycan_id) == 0) {
    return(invisible(NULL))
  }

  if (any(glycan_id > n_anomers)) {
    cli::cli_abort(
      "{.arg {arg_name}} contains {.field glycan_id} values outside the range defined by {.arg anomers}."
    )
  }

  invisible(NULL)
}


#' Build a glycan graph from one glycan's graph-table rows
#'
#' @param node_rows Node rows for one glycan.
#' @param edge_rows Edge rows for one glycan.
#' @param floating_rows Floating-part rows for one glycan.
#' @param floating_substituent_rows Floating-substituent rows for one glycan.
#' @param anomer Reducing-end anomer for the glycan.
#' @param alditol Alditol status for the glycan.
#' @param glycan_id Integer glycan ID used in error messages.
#' @returns An igraph object, or `NA` for a missing glycan.
#' @noRd
build_structure_graph_from_table_rows <- function(
  node_rows,
  edge_rows,
  floating_rows,
  floating_substituent_rows,
  anomer,
  alditol,
  glycan_id
) {
  if (nrow(node_rows) == 0) {
    if (nrow(edge_rows) > 0) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has edges without nodes."
      )
    }
    if (nrow(floating_rows) > 0) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has floating parts without nodes."
      )
    }
    if (nrow(floating_substituent_rows) > 0) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has floating substituents without nodes."
      )
    }
    if (!is.na(anomer)) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has an anomer but no nodes."
      )
    }

    return(NA)
  }

  if (is.na(anomer)) {
    cli::cli_abort(
      "Glycan {.val {glycan_id}} has nodes but a missing anomer."
    )
  }

  node_rows <- node_rows[order(node_rows$node_id), , drop = FALSE]
  edge_rows <- edge_rows[order(edge_rows$edge_id), , drop = FALSE]
  floating_rows <- floating_rows[
    order(floating_rows$part_id),
    ,
    drop = FALSE
  ]
  floating_substituent_rows <- floating_substituent_rows[
    order(floating_substituent_rows$substituent_id),
    ,
    drop = FALSE
  ]

  validate_consecutive_structure_ids(
    node_rows$node_id,
    "node_id",
    glycan_id
  )
  validate_consecutive_structure_ids(
    edge_rows$edge_id,
    "edge_id",
    glycan_id
  )
  validate_consecutive_structure_ids(
    floating_rows$part_id,
    "part_id",
    glycan_id
  )
  validate_consecutive_structure_ids(
    floating_substituent_rows$substituent_id,
    "substituent_id",
    glycan_id
  )
  validate_structure_edge_nodes(edge_rows, nrow(node_rows), glycan_id)

  graph <- igraph::make_empty_graph(n = nrow(node_rows), directed = TRUE)
  if (nrow(edge_rows) > 0) {
    graph <- igraph::add_edges(
      graph,
      c(rbind(edge_rows$from_node, edge_rows$to_node))
    )
  }

  igraph::V(graph)$name <- as.character(seq_len(igraph::vcount(graph)))
  igraph::V(graph)$mono <- node_rows$mono
  igraph::V(graph)$sub <- node_rows$sub
  igraph::E(graph)$linkage <- edge_rows$linkage
  graph$anomer <- anomer
  graph$alditol <- alditol
  parts <- lapply(seq_len(nrow(floating_rows)), function(row_id) {
    part <- list(
      root = floating_rows$root_node[[row_id]],
      linkage = floating_rows$linkage[[row_id]],
      parents = floating_rows$parents[[row_id]]
    )
    if ("nodes" %in% names(floating_rows)) {
      part$nodes <- floating_rows$nodes[[row_id]]
    }
    part
  })
  graph <- set_floating_parts_attr(graph, parts)
  substituents <- lapply(
    seq_len(nrow(floating_substituent_rows)),
    function(row_id) {
      list(
        substituent = floating_substituent_rows$substituent[[row_id]],
        parents = floating_substituent_rows$parents[[row_id]]
      )
    }
  )
  graph <- set_floating_substituents_attr(graph, substituents)

  graph
}


#' Validate that row IDs are consecutive within one glycan
#'
#' @param ids Integer IDs from one glycan's node or edge table.
#' @param column User-facing column name.
#' @param glycan_id Integer glycan ID used in error messages.
#' @returns Invisible `NULL`.
#' @noRd
validate_consecutive_structure_ids <- function(ids, column, glycan_id) {
  if (!identical(ids, seq_len(length(ids)))) {
    cli::cli_abort(
      "Glycan {.val {glycan_id}} column {.field {column}} must contain consecutive IDs starting at 1."
    )
  }

  invisible(NULL)
}


#' Validate edge endpoint IDs for one glycan
#'
#' @param edge_rows Edge rows for one glycan.
#' @param node_count Number of nodes for the glycan.
#' @param glycan_id Integer glycan ID used in error messages.
#' @returns Invisible `NULL`.
#' @noRd
validate_structure_edge_nodes <- function(edge_rows, node_count, glycan_id) {
  if (nrow(edge_rows) == 0) {
    return(invisible(NULL))
  }

  invalid_refs <- edge_rows$from_node > node_count |
    edge_rows$to_node > node_count
  if (any(invalid_refs)) {
    cli::cli_abort(
      "Glycan {.val {glycan_id}} edges must reference existing node IDs."
    )
  }

  invisible(NULL)
}
