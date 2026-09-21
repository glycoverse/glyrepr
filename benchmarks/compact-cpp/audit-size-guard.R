pkgload::load_all(quiet = TRUE)
source("benchmarks/compact-cpp/prototype.R")
x <- paste0(paste(rep("{Gal(a1-?)}", 2049), collapse = ""), "Glc")
a <- compact_arrays(x)[[1]]
stopifnot(a$status == "error", a$reason == "size guard")
writeLines(
  "PASS: complete floating forest exceeding 2048 residues hits native size guard",
  file.path(compact_dir, "results-v2", "size-guard.txt")
)
