pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
result_dir <- file.path(compact_dir, "results-v2")
dir.create(result_dir, showWarnings = FALSE)
files <- list.files("tests/testthat", "^test-.*\\.R$", full.names = TRUE)
files <- files[grepl("floating|substituent|iupac|alditol", files)]
literals <- unique(unlist(lapply(files, function(path) {
  d <- getParseData(parse(path, keep.source = TRUE))
  text <- d$text[d$token == "STR_CONST"]
  vapply(text, function(s) eval(parse(text = s)), "")
})))
modified <- as.vector(outer(
  compact_residues,
  c("6S", "6/3S", "?Me", "3S6Ac", "3/6S3/6Ac", "3S3Ac"),
  paste0
))
set.seed(20260922)
floating <- replicate(400, {
  count <- sample(1:3, 1)
  # Include valid and impossible assignments, cycles, singleton resolution,
  # links between floating parts, and conflicts with modified main residues.
  blocks <- vapply(
    seq_len(count),
    function(i) {
      parents <- sample(setdiff(seq_len(count + 3L), i), sample(1:3, 1))
      suffix <- if (sample(c(TRUE, FALSE), 1)) {
        paste0("|", paste(parents, collapse = ","))
      } else {
        ""
      }
      paste0(
        "{",
        sample(c("Gal", "Fuc", "Neu5Ac", "Gal6S"), 1),
        "(",
        sample(c("a1-3", "a1-?", "b1-3/6"), 1),
        ")",
        suffix,
        "}"
      )
    },
    ""
  )
  subs <- if (sample(c(TRUE, FALSE), 1)) {
    paste0(
      "{",
      sample(c("6S", "?S", "3/6Ac"), 1),
      "|",
      paste(sample(seq_len(count + 3L), 2), collapse = ","),
      "}"
    )
  } else {
    ""
  }
  paste0(
    subs,
    paste(blocks, collapse = ""),
    sample(c("Gal(b1-4)[Gal(b1-?)]Glc", "Gal6S(b1-4)[Gal(b1-?)]Glc"), 1)
  )
})
special <- c(
  "{Gal(a1-3)|2}{Fuc(a1-6)|3}Glc", # chained singleton resolution
  "{Gal(a1-3)|2}{Fuc(a1-6)|1}Glc", # closed cycle
  "{Gal(a1-3)|2,3}{Fuc(a1-6)|1,3}Glc",
  "{Gal(a1-3)|2}{Fuc(a1-6)}Glc",
  "{6S|1}Glc",
  "{6S}Glc",
  "{?S}{?S}Glc",
  "{6S|1,2}{Gal(a1-3)|2}Glc",
  "{Gal(a1-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc",
  "{Gal(a1-3)|3,4}Gal(a1-?)[Gal(a1-?)]Glc",
  "{Fuc(a1-?)}{Fuc(a1-?)}Gal",
  "{Gal(a1-3)|99}Glc",
  "{Gal(a1-3)|0}Glc",
  "{Gal(a1-3)|1}Glc",
  "{Gal(a1-3)|2,2}Glc",
  "{Gal(a1-3)|2147483648}Glc",
  "{6S|1}Gal6S",
  "{Gal-ol(a1-3)}Glc",
  "{6S|1,2}Gal(b1-4)Glc4S",
  "{6S|1,2}Gal(b1-4)Glc6S",
  "{Gal(b1-3/6)|2}Gal3S6Ac",
  "{Gal(??-?)|2}Gal3S6Ac"
)
x <- unique(c(literals, modified, floating, special))
saveRDS(x, file.path(result_dir, "extended-inputs.rds"))
a <- compact_arrays(x)
signature <- function(g) {
  list(
    edges = unname(igraph::as_edgelist(g, names = FALSE)),
    vertex = igraph::vertex_attr(g),
    edge = igraph::edge_attr(g),
    graph = igraph::graph_attr(g)
  )
}
rows <- vector("list", length(x))
for (i in seq_along(x)) {
  ref <- tryCatch(as_glycan_structure(x[i]), error = identity)
  ref_ok <- !inherits(ref, "error")
  native_ok <- a[[i]]$status == "ok"
  match <- identical(ref_ok, native_ok)
  if (ref_ok && native_ok) {
    cg <- compact_graph(a[[i]])
    match <- identical(a[[i]]$iupac, unname(as.character(ref))) &&
      identical(signature(cg), signature(attr(ref, "graphs")[[1]]))
  }
  rows[[i]] <- data.frame(
    index = i,
    reference_ok = ref_ok,
    native_ok = native_ok,
    match = match,
    native_reason = if (native_ok) "" else a[[i]]$reason
  )
  if (!match) {
    cat("MISMATCH", i, x[i], "reference", ref_ok, "native", native_ok, "\n")
  }
  if (i %% 250 == 0) cat("Extended audit", i, "/", length(x), "\n")
}
rows <- do.call(rbind, rows)
write.csv(rows, file.path(result_dir, "extended-audit.csv"), row.names = FALSE)
stopifnot(all(rows$match))
# Floating symmetry must preserve graph metadata as well as printed strings.
for (locale in c("C", "C.UTF-8", "en_US.UTF-8")) {
  if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) {
    for (s in special[rows$reference_ok[match(special, x)]]) {
      ref <- as_glycan_structure(s)
      actual <- compact_structure(s, fallback = FALSE)
      stopifnot(
        identical(as.character(ref), as.character(actual)),
        identical(
          signature(attr(ref, "graphs")[[1]]),
          signature(attr(actual, "graphs")[[1]])
        )
      )
    }
  }
}
writeLines(
  c(
    paste("inputs", nrow(rows)),
    paste("accepted", sum(rows$reference_ok)),
    paste("rejected", sum(!rows$reference_ok)),
    "PASS: native acceptance/rejection and exact canonical graph parity; no reference fallback"
  ),
  file.path(result_dir, "extended-audit.txt")
)
