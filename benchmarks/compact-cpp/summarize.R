result_dir <- file.path("benchmarks", "compact-cpp", "results")
files <- list.files(result_dir, "^timings-[0-9]+\\.csv$", full.names = TRUE)
x <- do.call(rbind, lapply(files, read.csv))
write.csv(x, file.path(result_dir, "timings.csv"), row.names = FALSE)
groups <- split(x, interaction(x$workload, x$stage, x$method, drop = TRUE))
s <- do.call(
  rbind,
  lapply(groups, function(g) {
    data.frame(
      workload = g$workload[1],
      stage = g$stage[1],
      method = g$method[1],
      n = nrow(g),
      median = median(g$elapsed),
      min = min(g$elapsed),
      max = max(g$elapsed)
    )
  })
)
write.csv(s, file.path(result_dir, "summary.csv"), row.names = FALSE)
pairs <- merge(
  subset(x, method == "R"),
  subset(x, method == "compact"),
  by = c("worker", "workload", "stage", "round")
)
pairs$speedup <- pairs$elapsed.x / pairs$elapsed.y
write.csv(
  pairs,
  file.path(result_dir, "paired-speedups.csv"),
  row.names = FALSE
)
print(s, row.names = FALSE)
print(aggregate(speedup ~ workload + stage, pairs, median))
workloads <- readRDS(file.path(result_dir, "workloads.rds"))
corpus <- readRDS(file.path(result_dir, "corpus.rds"))
coverage <- read.csv(file.path(result_dir, "coverage.csv"))
workload_coverage <- do.call(
  rbind,
  lapply(names(workloads), function(name) {
    ids <- match(workloads[[name]], corpus)
    data.frame(
      workload = name,
      rows = length(ids),
      distinct = length(unique(ids)),
      fast = sum(coverage$status[ids] == "ok"),
      floating = sum(coverage$reason[ids] == "floating metadata"),
      modified = sum(coverage$reason[ids] == "substituted or unknown residue")
    )
  })
)
write.csv(
  workload_coverage,
  file.path(result_dir, "workload-coverage.csv"),
  row.names = FALSE
)
paths <- c(
  list.files(
    file.path("benchmarks", "compact-cpp"),
    "\\.(R|cpp)$",
    full.names = TRUE
  ),
  file.path(result_dir, c("corpus.rds", "workloads.rds", "random-trees.rds"))
)
write.csv(
  data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
  file.path(result_dir, "provenance-md5.csv"),
  row.names = FALSE
)
