pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
result_dir <- file.path(compact_dir, "results-v2")
dir.create(result_dir, showWarnings = FALSE)
x <- as.character(glydb::glydb_structures())
saveRDS(x, file.path(result_dir, "corpus.rds"))
a <- compact_arrays(x)
status <- vapply(a, `[[`, "", "status")
reason <- vapply(a, function(z) if (is.null(z$reason)) "" else z$reason, "")
write.csv(
  data.frame(index = seq_along(x), status, reason),
  file.path(result_dir, "coverage.csv"),
  row.names = FALSE
)
print(table(status, reason))
stopifnot(all(status == "ok"))
graph_signature <- function(g) {
  list(
    directed = igraph::is_directed(g),
    edges = unname(igraph::as_edgelist(g, names = FALSE)),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
check_pair <- function(x) {
  baseline <- as_glycan_structure(x)
  candidate <- compact_structure(x, fallback = FALSE)
  stopifnot(
    identical(as.character(baseline), as.character(candidate)),
    identical(names(baseline), names(candidate)),
    identical(names(attr(baseline, "graphs")), names(attr(candidate, "graphs")))
  )
  bg <- attr(baseline, "graphs")
  cg <- attr(candidate, "graphs")
  for (key in names(bg)) {
    if (!identical(graph_signature(bg[[key]]), graph_signature(cg[[key]]))) {
      stop("graph mismatch: ", key)
    }
  }
  invisible(TRUE)
}
for (start in seq(1L, length(x), by = 250L)) {
  ids <- start:min(length(x), start + 249L)
  check_pair(x[ids])
  cat("Corpus parity:", max(ids), "/", length(x), "\n")
}
# Singletons, generic/concrete names, alditols, missing donor, and vector semantics.
check_pair(compact_residues)
check_pair(paste0(compact_residues, "-ol"))
check_pair(setNames(
  c(
    x[c(4, 4, 7)],
    NA_character_,
    "Man(a1-6)[Man(a1-3)]Man",
    "Man(a1-3)[Man(a1-6)]Man"
  ),
  letters[1:6]
))
check_pair(character())
check_pair(c(a = NA_character_, b = NA_character_))
fixtures <- c(
  "Gal(?-?)Glc",
  "Gal(b1-3/?)Glc",
  "Gal(b1-6/3)Glc",
  "Gal(b1-3/6)Glc",
  "Gal(b1-?)[Man(a1-?)]Glc",
  "Gal3S(b1-4)Glc",
  "Neu5,9Ac2",
  "Gal(b1-4)GlcNAc-ol(?1-"
)
for (s in fixtures) {
  ref <- tryCatch(as_glycan_structure(s), error = identity)
  if (!inherits(ref, "error")) check_pair(s)
}
# Seeded noncanonical trees: repeated acceptor positions are intentionally avoided.
set.seed(20260921)
make_tree <- function(n) {
  if (n == 1L) {
    return(sample(
      c("Glc", "Gal", "Man", "Fuc", "Hex", "Neu5Ac", "Galf", "D-Fuc"),
      1L
    ))
  }
  nk <- sample(seq_len(min(4L, n - 1L)), 1L)
  sizes <- rep(1L, nk)
  if (n - 1L - nk > 0L) {
    for (j in seq_len(n - 1L - nk)) {
      k <- sample.int(nk, 1L)
      sizes[k] <- sizes[k] + 1L
    }
  }
  kids <- vapply(sizes, make_tree, "")
  kids <- paste0(
    kids,
    "(",
    sample(c("a1-", "b1-", "?1-"), nk, replace = TRUE),
    sample(2:9, nk),
    ")"
  )
  paste0(
    kids[1],
    paste0(if (nk > 1) paste0("[", kids[-1], "]") else "", collapse = ""),
    make_tree(1L)
  )
}
random <- vapply(sample(2:40, 500, replace = TRUE), make_tree, "")
check_pair(random)
saveRDS(random, file.path(result_dir, "random-trees.rds"))
invalid <- c(
  "",
  " ",
  "Glc ",
  " Glc",
  "Gal(b1-4)[Man(a1-4)]Glc",
  "Gal(c1-4)Glc",
  "Gal(a3-4)Glc",
  "Gal(a1-0)Glc",
  "Gal(a1-4)",
  "[Gal(a1-4)]",
  "Glc[]",
  "GalGlc",
  "Gal(a1-4)Glc-ol-ol",
  "Gal-ol(a1-4)Glc",
  "Gal(a1-4))Glc",
  "[Gal(a1-4)Glc",
  "Gal(a1-4)]Glc"
)
for (s in invalid) {
  ref_error <- inherits(
    tryCatch(as_glycan_structure(s), error = identity),
    "error"
  )
  cpp_error <- inherits(
    tryCatch(compact_structure(s), error = identity),
    "error"
  )
  stopifnot(identical(ref_error, cpp_error))
}
# Test both pure C ordering and the active locale bridge on the same corpus.
original_locale <- Sys.getlocale("LC_COLLATE")
for (locale in unique(c("C", original_locale, "en_US.UTF-8"))) {
  chosen <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale))
  if (nzchar(chosen)) {
    check_pair(random)
    cat("Locale parity:", chosen, "\n")
  }
}
Sys.setlocale("LC_COLLATE", original_locale)
writeLines(
  c(
    paste("corpus", length(x)),
    paste("fast", sum(status == "ok")),
    paste("fallback", sum(status == "fallback")),
    paste("random", length(random)),
    paste("invalid", length(invalid)),
    "PASS: canonical strings, ordered graph attributes/endpoints, names, NA, duplicates, locales"
  ),
  file.path(result_dir, "parity.txt")
)
cat("PARITY PASS\n")
