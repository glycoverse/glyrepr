pkgload::load_all(quiet = TRUE)
for (path in Sys.glob(
  "benchmarks/array-backend/known-reference/failure-*.rds"
)) {
  d <- readRDS(path)
  cat("INDEX", d$index, "\n")
  g <- .parse_iupac_condensed_single(d$input)
  cat(
    "raw valid:",
    !inherits(tryCatch(validate_glycan_graph(g), error = identity), "error"),
    "\n"
  )
  cat("canonical graph validation:\n")
  print(tryCatch(
    {
      validate_glycan_graph(d$reference$graph)
      "valid"
    },
    error = conditionMessage
  ))
  a <- list(
    mono = igraph::V(g)$mono,
    sub = igraph::V(g)$sub,
    edges = as.integer(t(igraph::as_edgelist(g, names = FALSE))),
    linkage = igraph::E(g)$linkage,
    anomer = g$anomer,
    alditol = isTRUE(g$alditol),
    floating_parts = g$floating_parts,
    floating_substituents = g$floating_substituents
  )
  result <- structure_from_arrays(list(a))
  cat(
    "array vs reference key:",
    identical(as.character(result)[[1]], d$reference$iupac),
    "\n"
  )
  cat(d$reference$iupac, "\n")
}
