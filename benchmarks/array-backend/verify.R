# Run from package root: Rscript benchmarks/array-backend/verify.R
# Uses the frozen corpus already tracked with the compact-IUPAC benchmark.
pkgload::load_all(quiet = TRUE)
outdir <- "benchmarks/array-backend/results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
inputs <- unique(c(
  readRDS("benchmarks/compact-cpp/results-v2/corpus.rds"),
  readRDS("benchmarks/compact-cpp/results-v2/extended-inputs.rds")
))
record <- function(g) {
  list(
    mono = igraph::V(g)$mono,
    sub = igraph::V(g)$sub,
    edges = as.integer(t(igraph::as_edgelist(g, names = FALSE))),
    linkage = igraph::E(g)$linkage,
    anomer = g$anomer,
    alditol = isTRUE(g$alditol),
    floating_parts = g$floating_parts,
    floating_substituents = g$floating_substituents
  )
}
signature <- function(g) {
  list(
    vertices = igraph::vertex_attr(g),
    edges = igraph::as_edgelist(g, names = FALSE),
    edge_attributes = igraph::edge_attr(g),
    attributes = igraph::graph_attr(g)
  )
}
rows <- vector("list", length(inputs))
valid_records <- list()
for (i in seq_along(inputs)) {
  raw <- tryCatch(.parse_iupac_condensed_single(inputs[[i]]), error = identity)
  reference <- if (inherits(raw, "error")) {
    raw
  } else {
    tryCatch(process_glycan_structure_element(raw), error = identity)
  }
  if (inherits(reference, "error")) {
    rows[[i]] <- data.frame(
      index = i,
      reference_valid = FALSE,
      native_status = NA_character_,
      parity = NA
    )
    next
  }
  # Arbitrary input ordering: reverse the parser graph, remapping metadata
  # explicitly because igraph::permute does not understand graph attributes.
  a <- record(raw)
  n <- length(a$mono)
  a$mono <- rev(a$mono)
  a$sub <- rev(a$sub)
  a$edges <- n + 1L - a$edges
  a$floating_parts <- lapply(a$floating_parts, function(p) {
    p$root <- n + 1L - p$root
    p$nodes <- n + 1L - p$nodes
    p$parents <- n + 1L - p$parents
    p
  })
  a$floating_substituents <- lapply(a$floating_substituents, function(p) {
    p$parents <- n + 1L - p$parents
    p
  })
  native <- .compact_arrays_native(
    list(.validate_structure_arrays(a)),
    base::order,
    .compact_bliss_labels,
    identical(Sys.getlocale("LC_COLLATE"), "C")
  )[[1]]
  value <- tryCatch(structure_from_arrays(list(a)), error = identity)
  if (inherits(value, "error")) {
    saveRDS(
      list(
        index = i,
        input = inputs[[i]],
        record = a,
        reference = reference,
        native = native
      ),
      file.path(outdir, paste0("failure-", i, ".rds"))
    )
    cat("FAILURE", i, inputs[[i]], "\n")
    rows[[i]] <- data.frame(
      index = i,
      reference_valid = TRUE,
      native_status = native$status,
      parity = FALSE
    )
    next
  }
  parity <- identical(as.character(value)[[1]], reference$iupac) &&
    isTRUE(all.equal(
      signature(get_structure_graphs(value)),
      signature(reference$graph)
    ))
  rows[[i]] <- data.frame(
    index = i,
    reference_valid = TRUE,
    native_status = native$status,
    parity = parity
  )
  if (length(valid_records) < 500L) {
    valid_records[[length(valid_records) + 1L]] <- a
  }
  if (i %% 500L == 0L) cat(i, "checked\n")
}
audit <- do.call(rbind, rows)
write.csv(audit, file.path(outdir, "parity.csv"), row.names = FALSE)
saveRDS(valid_records, file.path(outdir, "workload.rds"))
writeLines(
  c(
    capture.output(sessionInfo()),
    paste("inputs", length(inputs)),
    paste("valid", sum(audit$reference_valid)),
    paste("mismatches", sum(!audit$parity, na.rm = TRUE)),
    paste("git_head", system("git rev-parse HEAD", intern = TRUE))
  ),
  file.path(outdir, "session.txt")
)
stopifnot(all(audit$parity[audit$reference_valid]))
cat("PARITY PASSED", sum(audit$reference_valid), "valid records\n")
