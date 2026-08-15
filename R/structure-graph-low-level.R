#' Validate a Glycan Graph
#'
#' Validate that a single `igraph` satisfies the structural and biochemical
#' requirements of a glycan graph. The graph is returned unchanged.
#'
#' This function does not canonicalize the graph or generate an
#' IUPAC-condensed string. Use [canonicalize_glycan_graph()] and
#' [graph_to_iupac()] for those operations.
#'
#' @param graph A single `igraph` glycan graph.
#'
#' @returns `graph`, unchanged. An error is thrown when `graph` is invalid.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
validate_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")

  if (!is_directed_graph(graph)) {
    cli::cli_abort("Glycan structure must be directed.")
  }

  validate_floating_graph_shape(graph)
  validate_floating_substituent_parents(graph)

  if (!has_vertex_attrs(graph, "mono")) {
    cli::cli_abort("Glycan structure must have a vertex attribute 'mono'.")
  }

  mono_names <- igraph::vertex_attr(graph, "mono")
  if (any(is.na(mono_names))) {
    cli::cli_abort(
      "Glycan structure must have no NA in vertex attribute 'mono'."
    )
  }

  if (!all(is_known_mono(mono_names))) {
    unknown_monos <- unique(igraph::V(graph)$mono[
      !is_known_mono(igraph::V(graph)$mono)
    ])
    msg <- glue::glue(
      "Unknown monosaccharide: {stringr::str_c(unknown_monos, collapse = ', ')}"
    )
    cli::cli_abort(msg, monos = unknown_monos)
  }

  if (!has_vertex_attrs(graph, "sub")) {
    cli::cli_abort("Glycan structure must have a vertex attribute 'sub'.")
  }

  subs <- igraph::vertex_attr(graph, "sub")
  if (any(is.na(subs))) {
    cli::cli_abort(
      "Glycan structure must have no NA in vertex attribute 'sub'."
    )
  }

  if (!all(valid_substituent(subs))) {
    invalid_subs <- unique(subs[!valid_substituent(subs)])
    msg <- glue::glue(
      "Unknown substituent: {stringr::str_c(invalid_subs, collapse = ', ')}"
    )
    cli::cli_abort(msg, subs = invalid_subs)
  }

  if (!has_edge_attrs(graph, "linkage")) {
    cli::cli_abort("Glycan structure must have an edge attribute 'linkage'.")
  }

  linkages <- igraph::edge_attr(graph, "linkage")
  if (any(is.na(linkages))) {
    cli::cli_abort(
      "Glycan structure must have no NA in edge attribute 'linkage'."
    )
  }

  if (!all(valid_linkages(linkages))) {
    invalid_linkages <- unique(linkages[!valid_linkages(linkages)])
    msg <- glue::glue(
      "Invalid linkage: {stringr::str_c(invalid_linkages, collapse = ', ')}"
    )
    cli::cli_abort(msg, linkages = invalid_linkages)
  }

  if (any_dup_linkage_pos(graph, linkages)) {
    cli::cli_abort("Duplicated linkage positions.")
  }

  validate_floating_metadata_assignments(graph)

  if (is.null(graph$anomer)) {
    cli::cli_abort("Glycan structure must have a graph attribute 'anomer'.")
  }

  if (!valid_anomer(graph$anomer)) {
    cli::cli_abort(glue::glue("Invalid anomer: {graph$anomer}"))
  }

  alditol <- igraph::graph_attr(graph, "alditol")
  if (
    !is.null(alditol) &&
      (!is.logical(alditol) || length(alditol) != 1 || is.na(alditol))
  ) {
    cli::cli_abort(
      "Glycan structure graph attribute {.field alditol} must be one non-missing logical value."
    )
  }

  graph
}


# Compatibility alias retained for downstream packages that used the former
# internal scalar validator.
validate_single_glycan_structure <- function(glycan) {
  validate_glycan_graph(glycan)
}


#' Canonicalize a Glycan Graph
#'
#' Add a vertex `name` attribute when needed and reorder the vertices and edges
#' of one glycan graph to match its IUPAC-condensed sequence.
#'
#' This function assumes that `graph` has already passed
#' [validate_glycan_graph()]. It performs no semantic validation.
#'
#' @inheritParams validate_glycan_graph
#'
#' @returns A canonicalized `igraph` glycan graph.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
canonicalize_glycan_graph <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  graph <- normalize_alditol_attr(graph)
  graph <- ensure_name_vertex_attr(graph)
  canonical_names <- as.character(seq_len(igraph::vcount(graph)))
  if (!identical(igraph::V(graph)$name, canonical_names)) {
    igraph::V(graph)$name <- canonical_names
  }
  .reorder_one_graph(graph)
}


#' Validate a List of Glycan Graphs
#'
#' Check the container used to store individually valid glycan graphs in one
#' glycan structure vector. Generic and concrete residues may coexist within a
#' graph and across graphs.
#'
#' This function assumes that every element has already passed
#' [validate_glycan_graph()]. It does not repeat scalar graph validation.
#'
#' @param graphs A list of individually valid `igraph` glycan graphs.
#' @param label An optional label retained for backward compatibility.
#'
#' @returns `NULL`, invisibly.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
validate_glycan_graph_vector <- function(graphs, label = NULL) {
  checkmate::assert_list(graphs, types = "igraph")
  checkmate::assert_string(label, null.ok = TRUE)
  invisible(NULL)
}


#' Generate IUPAC-Condensed from a Glycan Graph
#'
#' Generate one IUPAC-condensed string directly from one glycan graph.
#'
#' This low-level function assumes that `graph` is valid and canonical. It
#' performs no semantic validation or canonicalization.
#'
#' @inheritParams validate_glycan_graph
#'
#' @returns A single, unnamed IUPAC-condensed string.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
graph_to_iupac <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  raw_parts <- igraph::graph_attr(graph, "floating_parts")
  raw_substituents <- igraph::graph_attr(graph, "floating_substituents")
  if (is.null(raw_parts) && is.null(raw_substituents)) {
    seq_cache <- build_seq_cache(graph)
    root <- seq_cache$root
    return(format_reducing_end_iupac(
      seq_glycan_iupac(root, seq_cache),
      graph
    ))
  }

  parts <- normalize_floating_parts(graph)
  substituents <- normalize_floating_substituents(graph)

  if (length(parts) == 0 && length(substituents) == 0) {
    seq_cache <- build_seq_cache(graph)
    root <- seq_cache$root
    return(format_reducing_end_iupac(
      seq_glycan_iupac(root, seq_cache),
      graph
    ))
  }

  main_vertices <- floating_metadata_main_vertices(graph, parts)
  main <- igraph::induced_subgraph(graph, main_vertices)
  main <- delete_floating_parts_attr(main)
  main <- delete_floating_substituents_attr(main)
  seq_cache <- build_seq_cache(main)
  root <- seq_cache$root
  main_iupac <- format_reducing_end_iupac(
    seq_glycan_iupac(root, seq_cache),
    graph
  )
  floating_iupac <- purrr::map_chr(
    parts,
    floating_part_iupac,
    graph = graph
  )
  floating_sub_iupac <- purrr::map_chr(
    substituents,
    floating_substituent_iupac
  )
  floating_iupac <- c(floating_sub_iupac, floating_iupac)

  paste0(paste0(floating_iupac, collapse = ""), main_iupac)
}


#' Construct a Glycan Structure Vector from Trusted Data
#'
#' Assemble a glycan structure vector from IUPAC-condensed values and a graph
#' lookup table without graph validation, canonicalization, IUPAC generation,
#' vector-level compatibility checks, or graph deduplication.
#'
#' `graphs` must be a named list containing one graph for every distinct,
#' non-missing value in `iupac`. Additional named graphs are allowed so that
#' vctrs prototypes can retain their graph lookup tables. Graph names must be
#' unique and non-missing. This function checks these inexpensive representation
#' invariants but trusts that each graph matches its key.
#'
#' @param iupac A character vector of canonical IUPAC-condensed strings. Missing
#'   values are allowed, and names are preserved exactly.
#' @param graphs A named list of valid, canonical `igraph` glycan graphs keyed
#'   by IUPAC-condensed strings.
#'
#' @returns A `glyrepr_structure` vector.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
new_glycan_structure <- function(iupac = character(), graphs = list()) {
  checkmate::assert_character(iupac, any.missing = TRUE)
  checkmate::assert_list(graphs, types = "igraph")

  graph_names <- names(graphs)
  if (
    length(graphs) > 0 &&
      (is.null(graph_names) ||
        anyNA(graph_names) ||
        any(!nzchar(graph_names)) ||
        anyDuplicated(graph_names) > 0)
  ) {
    cli::cli_abort(
      "{.arg graphs} must have unique, non-missing IUPAC-condensed names."
    )
  }

  required_graphs <- unique(unname(iupac[!is.na(iupac)]))
  missing_graphs <- setdiff(required_graphs, graph_names)
  if (length(missing_graphs) > 0) {
    cli::cli_abort(c(
      "{.arg graphs} does not contain every structure in {.arg iupac}.",
      "x" = "Missing graph{?s}: {.val {missing_graphs}}."
    ))
  }

  nms <- names(iupac)
  out <- vctrs::new_vctr(
    as.list(unname(iupac)),
    graphs = graphs,
    class = "glyrepr_structure"
  )
  attr(out, "names") <- nms
  out
}
