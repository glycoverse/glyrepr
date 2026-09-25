# Run baseline before editing, then current in a fresh process:
# Rscript benchmarks/native-graph/audit.R baseline /tmp/glyrepr-native-baseline
# Rscript benchmarks/native-graph/audit.R current
args <- commandArgs(TRUE)
phase <- args[[1]]
if (length(args) > 1L) {
  library(glyrepr, lib.loc = args[[2]])
} else {
  pkgload::load_all(quiet = TRUE)
}
dir.create(
  "benchmarks/native-graph/results",
  recursive = TRUE,
  showWarnings = FALSE
)
root <- "benchmarks/native-graph/results"
if (phase == "baseline" && !file.exists(file.path(root, "corpus.rds"))) {
  strings <- character()
  collect <- function(x) {
    if (is.character(x)) {
      strings <<- c(strings, x)
    } else if (is.call(x) || is.expression(x) || is.pairlist(x)) {
      for (i in seq_along(x)) {
        try(collect(x[[i]]), silent = TRUE)
      }
    }
  }
  for (file in list.files("tests/testthat", "\\.R$", full.names = TRUE)) {
    collect(parse(file))
  }
  candidates <- unique(strings[
    grepl("\\([ab?][12?]-", strings) & nchar(strings) < 1000
  ])
  parsed <- suppressWarnings(as_glycan_structure(candidates, on_failure = "na"))
  corpus <- unique(unname(as.character(parsed[!is.na(parsed)])))
  saveRDS(corpus, file.path(root, "corpus.rds"))
} else {
  corpus <- readRDS(file.path(root, "corpus.rds"))
}
x <- as_glycan_structure(corpus)
graphs <- as.list(x)
# Attach attributes before permuting to expose incorrect ID remapping.
annotated <- lapply(graphs, function(g) {
  g <- igraph::set_vertex_attr(
    g,
    "source_id",
    value = seq_len(igraph::vcount(g))
  )
  g <- igraph::set_edge_attr(
    g,
    "weight",
    value = seq_len(igraph::ecount(g)) * 1.5
  )
  g$source <- list(format = "audit")
  if (has_floating_parts(g) || has_floating_substituents(g)) {
    g
  } else {
    igraph::permute(g, rev(seq_len(igraph::vcount(g))))
  }
})
describe_graph <- function(g) {
  if (is.null(g)) {
    return(NULL)
  }
  list(
    edges = igraph::as_edgelist(g, names = FALSE),
    vertices = igraph::vertex_attr(g),
    edge_attrs = igraph::edge_attr(g),
    graph_attrs = igraph::graph_attr(g)
  )
}
describe <- function(x) {
  if (inherits(x, "glyrepr_structure")) {
    return(list(
      values = as.character(x),
      graphs = lapply(attr(x, "graphs"), describe_graph)
    ))
  }
  if (inherits(x, "igraph")) {
    return(describe_graph(x))
  }
  if (is.list(x)) {
    return(list(attributes = attributes(x), data = lapply(x, describe)))
  }
  x
}
unknown <- remove_linkages(x)
floating <- as_glycan_structure(c(
  "{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-",
  "{3S|1,2}Gal(b1-4)Glc(?1-",
  "{Fuc(a1-2)|3,4}{Neu5Ac(a2-6)|3,4}Gal(b1-3)GalNAc(a1-"
))
tables <- list(
  nodes = structure_nodes(x),
  edges = structure_edges(x),
  floating_parts = structure_floating_parts(x),
  floating_substituents = structure_floating_substituents(x),
  anomers = get_anomer(x),
  alditols = get_alditol(x)
)
change <- function(g, ...) igraph::set_graph_attr(g, "anomer", "??")
cases <- list(
  generic = function() convert_to_generic(x),
  remove_linkages = function() remove_linkages(x),
  remove_substituents = function() remove_substituents(x),
  fill_anomer_pos = function() fill_anomer_pos(unknown),
  graph_constructor = function() as_glycan_structure(annotated),
  batch_canonicalize = function() canonicalize_glycan_graphs(annotated),
  canonicalize = function() lapply(annotated, canonicalize_glycan_graph),
  validate = function() lapply(annotated, validate_glycan_graph),
  serialize = function() lapply(graphs, graph_to_iupac),
  structure_to_iupac_graph = function() lapply(annotated, structure_to_iupac),
  tables = function() do.call(structure_from_tibbles, tables),
  smap = function() smap_structure(x, change),
  smap2 = function() smap2_structure(x, 1, change),
  spmap = function() spmap_structure(list(x, 1), change),
  simap = function() simap_structure(x, change),
  localize = function() {
    localize_floating_parts(
      floating[1],
      tibble::tibble(glycan_id = 1L, part_id = 1L, parent_node = 2L)
    )
  },
  enumerate = function() enumerate_floating_localizations(floating),
  enumerate_graph = function() {
    lapply(as.list(floating), enumerate_floating_graph_localizations)
  }
)
outputs <- lapply(cases, function(f) describe(f()))
exports <- sort(getNamespaceExports("glyrepr"))
outputs$api <- lapply(stats::setNames(exports, exports), function(name) {
  object <- getExportedValue("glyrepr", name)
  if (is.function(object)) formals(object) else class(object)
})
condition_data <- function(f) {
  warnings <- list()
  value <- withCallingHandlers(
    tryCatch(f(), error = identity),
    warning = function(cnd) {
      warnings[[length(warnings) + 1L]] <<- list(
        class = class(cnd),
        message = conditionMessage(cnd),
        positions = cnd$positions,
        reasons = cnd$reasons
      )
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(value, "error")) {
    value <- list(
      class = class(value),
      message = conditionMessage(value),
      position = value$position,
      input_name = value$input_name,
      reason = value$reason
    )
  }
  list(value = describe(value), warnings = warnings)
}
bad <- graphs[[1]]
igraph::V(bad)$mono[1] <- "unknown"
outputs$conditions <- list(
  validate = condition_data(function() validate_glycan_graph(bad)),
  strict = condition_data(function() {
    as_glycan_structure(list(graphs[[1]], bad))
  }),
  recovery = condition_data(function() {
    as_glycan_structure(
      list(ok = graphs[[1]], bad = bad, missing = NULL),
      on_failure = "na"
    )
  }),
  batch = condition_data(function() {
    canonicalize_glycan_graphs(
      list(ok = graphs[[1]], bad = bad, missing = NULL),
      on_failure = "na"
    )
  })
)
# Exercise every bounded floating example in the corpus, beyond the timing trio.
float_ids <- which(has_floating_parts(x) | has_floating_substituents(x))
outputs$floating_corpus <- lapply(float_ids, function(i) {
  condition_data(function() {
    enumerate_floating_localizations(x[i], max_variants = 64)
  })
})
saveRDS(outputs, file.path(root, paste0(phase, "-outputs.rds")))
if (phase != "baseline") {
  baseline <- readRDS(file.path(root, "baseline-outputs.rds"))
  diffs <- lapply(names(outputs), function(name) {
    diff <- waldo::compare(baseline[[name]], outputs[[name]], max_diffs = 10)
    if (length(diff)) {
      paste(name, paste(diff, collapse = "\n"), sep = "\n")
    } else {
      NULL
    }
  })
  names(diffs) <- names(outputs)
  writeLines(
    as.character(unlist(diffs)),
    file.path(root, "parity-differences.txt")
  )
  if (any(lengths(diffs))) {
    stop(
      "Parity differences: ",
      paste(names(diffs)[lengths(diffs) > 0], collapse = ", ")
    )
  }
}
if (identical(Sys.getenv("GLYREPR_AUDIT_PARITY_ONLY"), "1")) {
  cat(
    phase,
    "parity complete:",
    length(corpus),
    "structures;",
    length(exports),
    "exports;",
    length(float_ids),
    "floating examples\n"
  )
  quit(status = 0)
}
timing <- lapply(names(cases), function(name) {
  f <- cases[[name]]
  repetitions <- if (name %in% c("localize", "enumerate", "enumerate_graph")) {
    10L
  } else {
    1L
  }
  elapsed <- replicate(5, {
    gc()
    unname(system.time(
      for (i in seq_len(repetitions)) {
        f()
      }
    )[["elapsed"]]) /
      repetitions
  })
  data.frame(
    case = name,
    median_seconds = median(elapsed),
    min_seconds = min(elapsed),
    max_seconds = max(elapsed)
  )
})
write.csv(
  do.call(rbind, timing),
  file.path(root, paste0(phase, "-timings.csv")),
  row.names = FALSE
)
writeLines(
  c(capture.output(sessionInfo()), paste("corpus_size:", length(corpus))),
  file.path(root, paste0(phase, "-session.txt"))
)
cat(
  phase,
  "audit complete:",
  length(corpus),
  "structures;",
  length(cases),
  "API workloads\n"
)
