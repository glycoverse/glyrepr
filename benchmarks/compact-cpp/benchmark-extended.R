# Run one worker in a fresh R process: Rscript benchmarks/compact-cpp/benchmark.R 1
pkgload::load_all(quiet = TRUE)
compile_time <- system.time(source("benchmarks/compact-cpp/prototype.R"))[[
  "elapsed"
]]
worker <- as.integer(commandArgs(TRUE)[1])
result_dir <- file.path(compact_dir, "results-v2")
x <- readRDS(file.path(compact_dir, "results", "corpus.rds"))
coverage <- read.csv(file.path(compact_dir, "results", "coverage.csv"))
workloads <- readRDS(file.path(compact_dir, "results", "workloads.rds"))
set.seed(20260922)
workloads$modified_500 <- x[sample(
  which(coverage$reason == "substituted or unknown residue"),
  500
)]
workloads$floating_564 <- x[which(coverage$reason == "floating metadata")]
if (worker == 1L) {
  saveRDS(workloads, file.path(result_dir, "workloads.rds"))
}
rows <- list()
record <- function(workload, stage, method, round, f) {
  gc()
  elapsed <- system.time(value <- f())[["elapsed"]]
  rows[[length(rows) + 1L]] <<- data.frame(
    worker,
    workload,
    stage,
    method,
    round,
    elapsed
  )
  invisible(value)
}
for (name in names(workloads)) {
  input <- workloads[[name]]
  # Warm both implementations outside measured regions.
  invisible(as_glycan_structure(input[1:3]))
  invisible(compact_structure(input[1:3], fallback = FALSE))
  for (round in 1:2) {
    methods <- if ((worker + round) %% 2L) {
      c("R", "compact")
    } else {
      c("compact", "R")
    }
    for (method in methods) {
      record(
        name,
        "public_constructor",
        method,
        round,
        if (method == "R") {
          function() as_glycan_structure(input)
        } else {
          function() compact_structure(input, fallback = FALSE)
        }
      )
    }
  }
  cat("worker", worker, name, "complete\n")
}
input <- workloads$fast_500
arrays <- compact_arrays(input)
for (round in 1:2) {
  # Arrays exclude graph creation; graph_materialization excludes parsing. No speedup
  # ratio is asserted between these different-output stages.
  record("fast_500", "arrays_only", "compact", round, function() {
    compact_arrays(input)
  })
  record("fast_500", "graph_materialization", "compact", round, function() {
    lapply(arrays, compact_graph)
  })
  for (method in if ((worker + round) %% 2L) {
    c("R", "compact")
  } else {
    c("compact", "R")
  }) {
    record(
      "fast_500",
      "parse_validate_canonical_pairs",
      method,
      round,
      if (method == "R") {
        function() {
          lapply(input, function(s) {
            process_glycan_structure_element(.parse_iupac_condensed_single(s))
          })
        }
      } else {
        function() {
          lapply(compact_arrays(input), function(a) {
            list(graph = compact_graph(a), iupac = a$iupac)
          })
        }
      }
    )
  }
}
write.csv(
  do.call(rbind, rows),
  file.path(result_dir, paste0("timings-", worker, ".csv")),
  row.names = FALSE
)
writeLines(
  c(
    paste("compile_and_source_seconds", compile_time),
    capture.output(sessionInfo()),
    paste("git_head", system("git rev-parse HEAD", intern = TRUE)),
    paste("system", paste(Sys.info(), collapse = " "))
  ),
  file.path(result_dir, paste0("session-", worker, ".txt"))
)
if (worker == 1L) {
  for (method in c("R", "compact")) {
    path <- file.path(result_dir, paste0("profile-", method, ".out"))
    Rprof(path, interval = 0.001)
    invisible(
      if (method == "R") {
        as_glycan_structure(input)
      } else {
        compact_structure(input, fallback = FALSE)
      }
    )
    Rprof(NULL)
    write.csv(
      summaryRprof(path)$by.total,
      file.path(result_dir, paste0("profile-", method, ".csv"))
    )
  }
}
