dir <- file.path("benchmarks", "compact-cpp", "results-api")
files <- list.files(dir, "^timings-(R|compact)-[1-3]\\.csv$", full.names = TRUE)
x <- do.call(rbind, lapply(files, read.csv))
stopifnot(length(files) == 6L, nrow(x) == 36L)
write.csv(x, file.path(dir, "timings.csv"), row.names = FALSE)
s <- do.call(
  rbind,
  lapply(split(x, interaction(x$workload, x$method)), function(g) {
    data.frame(
      workload = g$workload[1],
      method = g$method[1],
      n = nrow(g),
      median = median(g$elapsed),
      min = min(g$elapsed),
      max = max(g$elapsed)
    )
  })
)
write.csv(s, file.path(dir, "summary.csv"), row.names = FALSE)
pairs <- merge(
  subset(x, method == "R"),
  subset(x, method == "compact"),
  by = c("worker", "workload", "round")
)
pairs$speedup <- pairs$elapsed.x / pairs$elapsed.y
write.csv(pairs, file.path(dir, "paired-speedups.csv"), row.names = FALSE)
a <- readRDS(file.path(dir, "public-baseline.rds"))
b <- readRDS(file.path(dir, "public-candidate.rds"))
stopifnot(identical(a$api, b$api), identical(a$cases, b$cases))
writeLines(
  sprintf(
    "PASS: %d exported signatures and %d behavior/condition cases identical",
    length(a$api),
    length(a$cases)
  ),
  file.path(dir, "public-parity.txt")
)
paths <- c(
  list.files("src", "\\.(cpp|win)$|Makevars$", full.names = TRUE),
  c(
    "DESCRIPTION",
    "NAMESPACE",
    "R/compact-iupac.R",
    "R/structure.R",
    "R/RcppExports.R"
  ),
  list.files(dir, "\\.rds$", full.names = TRUE)
)
write.csv(
  data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
  file.path(dir, "provenance-md5.csv"),
  row.names = FALSE
)
print(s, row.names = FALSE)
