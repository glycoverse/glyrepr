test_that("the cached graph template is detached from native igraph state", {
  model <- .compact_graph_template()
  skip_if(is.null(model), "Installed igraph uses an unsupported graph format")
  expect_null(unclass(model)[[10L]][["igraph"]])

  restored <- unserialize(serialize(.compact_graph_template, NULL))
  restored_model <- restored()
  expect_null(unclass(restored_model)[[10L]][["igraph"]])
  arrays <- .compact_iupac_arrays(c("Gal", "Gal(b1-4)Glc"))
  actual <- .compact_build_graphs_native(arrays, restored_model)
  expected <- lapply(arrays, .compact_build_graph)
  expect_identical(
    lapply(actual, igraph::as_data_frame, what = "both"),
    lapply(expected, igraph::as_data_frame, what = "both")
  )
})

test_that("the graph template can be serialized by an exit finalizer", {
  # A new process cannot see pkgload's in-memory development namespace.
  skip_if_not(
    file.exists(file.path(
      getNamespaceInfo(asNamespace("glyrepr"), "path"),
      "Meta",
      "package.rds"
    )),
    "Requires the installed package; run R CMD check"
  )
  script <- tempfile(fileext = ".R")
  libpaths <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  log <- tempfile(fileext = ".log")
  on.exit(unlink(c(script, libpaths, output, log)), add = TRUE)
  saveRDS(.libPaths(), libpaths)
  # Exit finalizers run in reverse registration order. Register the saver
  # before creating the graph, reproducing covr's shutdown serialization.
  writeLines(
    c(
      "args <- commandArgs(trailingOnly = TRUE)",
      ".libPaths(readRDS(args[[1L]]))",
      "holder <- new.env(parent = emptyenv())",
      "holder$path <- args[[2L]]",
      "reg.finalizer(holder, function(e) saveRDS(e$fun, e$path), onexit = TRUE)",
      'holder$fun <- getFromNamespace(".compact_graph_template", "glyrepr")',
      "stopifnot(!isNamespace(environment(holder$fun)))",
      "invisible(holder$fun())",
      'q(save = "no", status = 0L, runLast = FALSE)'
    ),
    script
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), shQuote(libpaths), shQuote(output)),
    stdout = log,
    stderr = log
  )
  expect_identical(
    status,
    0L,
    info = paste(readLines(log, warn = FALSE), collapse = "\n")
  )
  expect_true(file.exists(output))
  restored <- readRDS(output)
  expect_type(restored, "closure")
  model <- environment(restored)$model
  if (!identical(model, FALSE)) {
    expect_null(unclass(model)[[10L]][["igraph"]])
  }
})

test_that("batch graph construction preserves igraph behavior and serialization", {
  x <- c(
    "Gal(b1-4)GlcNAc-ol",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc",
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-",
    "{6S|1,2}Gal(b1-4)Glc",
    "Gal"
  )
  arrays <- .compact_iupac_arrays(x)
  graphs <- .compact_build_graphs(arrays)
  reference <- lapply(arrays, .compact_build_graph)
  describe <- function(g) {
    list(
      edges = igraph::as_edgelist(g, names = FALSE),
      vertices = igraph::vertex_attr(g),
      edge_attrs = igraph::edge_attr(g),
      attrs = igraph::graph_attr(g)
    )
  }
  # Check an untouched serialized value as well as one with a materialized
  # igraph C index. Both must rebuild independent indexes after deserialization.
  restored <- unserialize(serialize(graphs, NULL))
  expect_identical(lapply(graphs, describe), lapply(reference, describe))
  expect_identical(lapply(restored, describe), lapply(reference, describe))
  restored <- unserialize(serialize(graphs, NULL))
  expect_identical(lapply(restored, describe), lapply(reference, describe))
  for (i in seq_along(graphs)) {
    expect_identical(
      igraph::degree(graphs[[i]]),
      igraph::degree(reference[[i]])
    )
    expect_identical(
      describe(igraph::permute(
        graphs[[i]],
        rev(seq_len(igraph::vcount(graphs[[i]])))
      )),
      describe(igraph::permute(
        reference[[i]],
        rev(seq_len(igraph::vcount(reference[[i]])))
      ))
    )
  }
  ids <- vapply(
    c(graphs, .compact_build_graphs(arrays)),
    igraph::graph_id,
    character(1)
  )
  expect_identical(anyDuplicated(ids), 0L)
  changed <- igraph::set_vertex_attr(graphs[[1]], "mono", 1, "Man")
  changed <- igraph::add_vertices(changed, 1)
  expect_identical(describe(graphs[[1]]), describe(reference[[1]]))
  expect_equal(igraph::vcount(changed), igraph::vcount(graphs[[1]]) + 1)
})

test_that("unknown igraph object formats use the public graph builder", {
  arrays <- .compact_iupac_arrays(c("Gal", "Gal(b1-4)Glc"))
  expect_null(.compact_build_graphs_native(arrays, NULL))
  expect_null(.compact_build_graphs_native(arrays, list()))
  local_mocked_bindings(.compact_graph_template = function() NULL)
  actual <- .compact_build_graphs(arrays)
  expected <- lapply(arrays, .compact_build_graph)
  expect_identical(
    lapply(actual, igraph::as_data_frame, what = "both"),
    lapply(expected, igraph::as_data_frame, what = "both")
  )
  expect_identical(
    .compact_bliss_labels(c(1L, 2L), c(1L, 1L)),
    as.integer(
      igraph::canonical_permutation(
        igraph::make_graph(c(1, 2), directed = FALSE),
        colors = c(1, 1)
      )$labeling
    )
  )
})

test_that("floating graph callbacks remain rooted during garbage collection", {
  input <- c(
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-",
    "{6S|1,2}Gal(b1-4)Glc",
    "{6S|1,2}{3P|1,2}Gal(b1-4)Glc"
  )
  reference <- .compact_iupac_arrays(input)
  previous <- gctorture2(100L)
  on.exit(gctorture2(previous), add = TRUE)
  actual <- lapply(seq_len(4L), function(i) .compact_iupac_arrays(input))
  gctorture2(previous)
  for (result in actual) {
    expect_identical(result, reference)
  }
})
