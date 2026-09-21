pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
# Mine literal strings from existing tests as an independent grammar/error corpus.
files <- c(
  "tests/testthat/test-iupac-to-structure.R",
  "tests/testthat/test-structure-to-iupac.R"
)
literals <- unique(unlist(lapply(files, function(path) {
  d <- getParseData(parse(path, keep.source = TRUE))
  text <- d$text[d$token == "STR_CONST"]
  vapply(text, function(s) eval(parse(text = s)), "")
})))
a <- compact_arrays(literals)
ok <- which(vapply(a, function(z) z$status == "ok", TRUE))
for (i in ok) {
  reference <- as_glycan_structure(literals[i])
  stopifnot(identical(unname(as.character(reference)), a[[i]]$iupac))
}
writeLines(
  c(
    paste("test_literals", length(literals)),
    paste("fast_accepted", length(ok)),
    "PASS: every fast-accepted literal accepted by reference with identical canonical string"
  ),
  file.path(compact_dir, "results", "syntax-audit.txt")
)
