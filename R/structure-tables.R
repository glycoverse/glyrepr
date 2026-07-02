#' Convert Glycan Structures to Graph Tables
#'
#' @description
#' `structure_nodes()` and `structure_edges()` convert a glycan structure vector
#' to node and edge tibbles. `structure_from_tibbles()` rebuilds a
#' `glyrepr_structure` vector from those tibbles and a vector of reducing-end
#' anomers.
#'
#' The `glycan_id` column is the integer position of each glycan in the input
#' vector. Duplicate structures are expanded to their original vector positions.
#' Missing structures have no node or edge rows and are reconstructed from
#' missing values in `anomers`.
#' If `x` is named, the node and edge tibbles also contain a `glycan_name`
#' column. `structure_from_tibbles()` uses `glycan_name` as output names when
#' that column is present.
#'
#' @param x A glycan structure vector.
#' @param nodes A data frame with columns `glycan_id`, `node_id`, `mono`, and
#'   `sub`, and optionally `glycan_name`.
#' @param edges A data frame with columns `glycan_id`, `edge_id`, `from_node`,
#'   `to_node`, and `linkage`, and optionally `glycan_name`.
#' @param anomers A character vector of reducing-end anomers, one per glycan.
#'
#' @returns
#' - `structure_nodes()` returns a tibble with columns `glycan_id`, `node_id`,
#'   `mono`, and `sub`.
#' - `structure_edges()` returns a tibble with columns `glycan_id`, `edge_id`,
#'   `from_node`, `to_node`, and `linkage`.
#' - `structure_from_tibbles()` returns a `glyrepr_structure` vector.
#'
#' @examples
#' glycans <- c(o_glycan_core_1(), o_glycan_core_1())
#' nodes <- structure_nodes(glycans)
#' edges <- structure_edges(glycans)
#' structure_from_tibbles(nodes, edges, get_anomer(glycans))
#'
#' @name structure_tables
NULL


#' @rdname structure_tables
#' @export
structure_nodes <- function(x) {
  checkmate::assert_class(x, "glyrepr_structure")

  graphs <- as.list(x)
  glycan_names <- names(x)
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
  checkmate::assert_class(x, "glyrepr_structure")

  graphs <- as.list(x)
  glycan_names <- names(x)
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
structure_from_tibbles <- function(nodes, edges, anomers) {
  nodes <- validate_structure_nodes_table(nodes)
  edges <- validate_structure_edges_table(edges)
  anomers <- validate_structure_anomers(anomers)

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
  glycan_names <- structure_table_glycan_names(nodes, edges, anomers)

  if (length(anomers) == 0) {
    return(glycan_structure())
  }

  graphs <- purrr::map(seq_along(anomers), function(glycan_id) {
    build_structure_graph_from_table_rows(
      nodes[nodes$glycan_id == glycan_id, , drop = FALSE],
      edges[edges$glycan_id == glycan_id, , drop = FALSE],
      anomers[[glycan_id]],
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
#' @param anomers A validated anomer vector.
#' @returns A character vector of names, or `NULL`.
#' @noRd
structure_table_glycan_names <- function(nodes, edges, anomers) {
  table_names <- structure_table_glycan_names_from_rows(nodes, edges, anomers)
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
#' @param anomers A validated anomer vector.
#' @returns A character vector with `NA` where no table name is available, or
#'   `NULL` when neither table contains `glycan_name`.
#' @noRd
structure_table_glycan_names_from_rows <- function(nodes, edges, anomers) {
  name_rows <- dplyr::bind_rows(
    structure_table_name_rows(nodes),
    structure_table_name_rows(edges)
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
#' @param anomer Reducing-end anomer for the glycan.
#' @param glycan_id Integer glycan ID used in error messages.
#' @returns An igraph object, or `NA` for a missing glycan.
#' @noRd
build_structure_graph_from_table_rows <- function(
  node_rows,
  edge_rows,
  anomer,
  glycan_id
) {
  if (nrow(node_rows) == 0) {
    if (nrow(edge_rows) > 0) {
      cli::cli_abort(
        "Glycan {.val {glycan_id}} has edges without nodes."
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
