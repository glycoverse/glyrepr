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

  if (!is_out_tree(graph)) {
    cli::cli_abort("Glycan structure must be an out tree.")
  }

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

  if (mix_generic_concrete(mono_names)) {
    cli::cli_abort(
      "Monosaccharides must be either all generic or all concrete."
    )
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

  if (any_dup_linkage_pos(graph)) {
    cli::cli_abort("Duplicated linkage positions.")
  }

  if (is.null(graph$anomer)) {
    cli::cli_abort("Glycan structure must have a graph attribute 'anomer'.")
  }

  if (!valid_anomer(graph$anomer)) {
    cli::cli_abort(glue::glue("Invalid anomer: {graph$anomer}"))
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
  graph <- ensure_name_vertex_attr(graph)
  .reorder_one_graph(graph)
}


#' Validate Compatibility Across Glycan Graphs
#'
#' Check that a list of individually valid glycan graphs can coexist in one
#' glycan structure vector. All graphs must use the same monosaccharide type:
#' either concrete or generic.
#'
#' This function assumes that every element has already passed
#' [validate_glycan_graph()]. It does not repeat scalar graph validation.
#'
#' @param graphs A list of individually valid `igraph` glycan graphs.
#' @param label An optional label used in error messages.
#'
#' @returns `NULL`, invisibly. An error is thrown when the graphs are
#'   incompatible.
#'
#' @template low-level-structure-pipeline
#'
#' @family low-level glycan structure functions
#' @export
validate_glycan_graph_vector <- function(graphs, label = NULL) {
  checkmate::assert_list(graphs, types = "igraph")
  checkmate::assert_string(label, null.ok = TRUE)

  if (length(graphs) <= 1) {
    return(invisible(NULL))
  }

  mono_types <- purrr::map_chr(graphs, get_graph_mono_type)

  if (any(mono_types == "mixed")) {
    cli::cli_abort(c(
      "All structures must have a single monosaccharide type.",
      "x" = "{.val {label}} contains structures with mixed generic and concrete monosaccharides."
    ))
  }

  unique_types <- unique(mono_types)
  if (length(unique_types) > 1) {
    concrete_count <- sum(mono_types == "concrete")
    generic_count <- sum(mono_types == "generic")

    if (is.null(label)) {
      cli::cli_abort(c(
        "All structures must have the same monosaccharide type.",
        "x" = "Found {.val {concrete_count}} concrete and {.val {generic_count}} generic structure(s) in the same vector.",
        "i" = "Use {.fn convert_to_generic} to convert concrete structures to generic type."
      ))
    } else {
      cli::cli_abort(c(
        "All structures must have the same monosaccharide type.",
        "x" = "{.val {label}} has mixed types: {.val {concrete_count}} concrete and {.val {generic_count}} generic structure(s)."
      ))
    }
  }

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
  root <- which(igraph::degree(graph, mode = "in") == 0)
  seq_cache <- build_seq_cache(graph, root)
  paste0(seq_glycan_iupac(root, seq_cache), "(", graph$anomer, "-")
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
#' @param graphs A named list of valid, canonical, mutually compatible `igraph`
#'   glycan graphs keyed by IUPAC-condensed strings.
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
