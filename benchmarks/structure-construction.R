# Usage:
# Rscript benchmarks/structure-construction.R prepare SOURCE INPUT.rds
# Rscript benchmarks/structure-construction.R measure SOURCE INPUT.rds OUTPUT.rds
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L, args[1] %in% c("prepare", "measure"))
devtools::load_all(args[2], quiet = TRUE)

if (args[1] == "prepare") {
  strings <- unique(as.character(glydb::glydb_structures()))
  ordinary <- strings[!startsWith(strings, "{") & !is.na(strings)]
  set.seed(20260908)
  strings <- sample(ordinary, min(500L, length(ordinary)))
  parse_one <- getFromNamespace(".parse_iupac_condensed_single", "glyrepr")
  canonical <- as_glycan_structure(strings)
  input <- list(
    strings = strings,
    graphs = lapply(strings, parse_one),
    canonical_graphs = as.list(canonical),
    iupac = as.character(canonical),
    lookup = attr(canonical, "graphs"),
    corpus_size = length(ordinary),
    glydb_version = as.character(packageVersion("glydb"))
  )
  saveRDS(input, args[3])
} else {
  stopifnot(length(args) == 4L)
  input <- readRDS(args[3])
  cases <- list(
    character = function() as_glycan_structure(input$strings),
    raw_graphs = function() as_glycan_structure(input$graphs),
    canonical_graphs = function() as_glycan_structure(input$canonical_graphs),
    repeated_character = function() {
      as_glycan_structure(rep(input$strings[1:10], 100L))
    },
    repeated_graphs = function() {
      as_glycan_structure(rep(input$graphs[1:10], 100L))
    },
    recover_character = function() {
      as_glycan_structure(input$strings, on_failure = "na")
    },
    trusted = function() new_glycan_structure(input$iupac, input$lookup),
    all_missing = function() do.call(glycan_structure, rep(list(NA), 10000L))
  )
  signature <- function(x) {
    list(
      iupac = as.character(x),
      names = names(x),
      graphs = lapply(attr(x, "graphs"), function(graph) {
        list(
          vertices = igraph::as_data_frame(graph, "vertices"),
          edges = igraph::as_data_frame(graph, "edges"),
          attributes = igraph::graph_attr(graph)
        )
      })
    )
  }
  invisible(as_glycan_structure(input$strings[1:5]))
  timings <- signatures <- list()
  for (case in names(cases)) {
    gc()
    time <- system.time(result <- cases[[case]]())
    timings[[case]] <- data.frame(
      case = case,
      elapsed = unname(time["elapsed"]),
      cpu = unname(sum(time[c("user.self", "sys.self")]))
    )
    signatures[[case]] <- signature(result)
    cat(case, time["elapsed"], "s\n")
  }
  allocation_file <- tempfile()
  gc()
  Rprofmem(allocation_file)
  invisible(cases$all_missing())
  Rprofmem(NULL)
  allocations <- suppressWarnings(as.numeric(sub(
    " .*",
    "",
    readLines(allocation_file)
  )))
  unlink(allocation_file)
  saveRDS(
    list(
      timings = do.call(rbind, timings),
      signatures = signatures,
      missing_allocated_mib = sum(allocations, na.rm = TRUE) / 1024^2,
      session = sessionInfo()
    ),
    args[4]
  )
}
