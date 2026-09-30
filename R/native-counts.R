.native_count_schema <- local({
  schema <- NULL
  function() {
    if (is.null(schema)) {
      known <- available_monosaccharides()
      schema <<- list(
        known = known,
        components = .composition_component_order(),
        generic = unique(monosaccharides$generic),
        converted = convert_mono_type_impl(known),
        concrete = setdiff(monosaccharides$concrete, monosaccharides$generic)
      )
    }
    schema
  }
})

.native_graph_counts <- function(
  graphs,
  mono = NULL,
  include_subs = FALSE,
  count = FALSE
) {
  schema <- .native_count_schema()
  records <- lapply(graphs, function(graph) {
    residues <- igraph::vertex_attr(graph, "mono")
    subs <- igraph::vertex_attr(graph, "sub")
    if (
      !is.character(residues) ||
        !length(residues) ||
        !is.null(attributes(residues)) ||
        anyNA(residues) ||
        !all(residues %in% schema$known) ||
        (!is.null(subs) && (!is.character(subs) || !is.null(attributes(subs))))
    ) {
      return(NULL)
    }
    floating <- normalize_floating_substituents(graph)
    subs <- c(subs, vapply(floating, `[[`, character(1), "substituent"))
    list(residues, subs)
  })
  if (any(vapply(records, is.null, logical(1)))) {
    return(NULL)
  }
  known <- schema$known
  generic <- schema$generic
  targets <- if (is.null(mono)) {
    if (include_subs) schema$components else known
  } else if (mono %in% generic) {
    known[schema$converted == mono]
  } else {
    mono
  }
  uncertain <- !is.null(mono) &&
    mono %in% schema$concrete
  .graph_counts_native(
    records,
    schema$components,
    targets,
    unique(generic),
    count,
    uncertain
  )
}

.native_structure_compositions <- function(
  x,
  mono = NULL,
  include_subs = FALSE,
  count = FALSE
) {
  codes <- glycan_structure_iupac_data(x)
  used <- unique(codes[!is.na(codes)])
  result <- .native_graph_counts(
    unname(attr(x, "graphs")[used]),
    mono,
    include_subs,
    count
  )
  if (is.null(result)) {
    return(NULL)
  }
  index <- match(codes, used)
  result <- result[index]
  names(result) <- names(x)
  result
}
