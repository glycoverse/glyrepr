test_that("batch graph canonicalization retains per-input attributes", {
  g <- .parse_iupac_condensed_single("Gal(b1-4)[Fuc(a1-3)]GlcNAc")
  g <- igraph::set_vertex_attr(
    g,
    "source_id",
    value = seq_len(igraph::vcount(g))
  )
  g <- igraph::set_edge_attr(g, "weight", value = c(11, 22))
  g$source <- "first"
  other <- g
  other$source <- "second"
  input <- list(first = g, missing = NULL, second = other)
  result <- canonicalize_glycan_graphs(input)
  expect_identical(
    result$status,
    c(first = "ok", missing = "missing", second = "ok")
  )
  expect_null(result$graphs$missing)
  expect_identical(result$graphs$first$source, "first")
  expect_identical(result$graphs$second$source, "second")
  for (i in c(1, 3)) {
    expected <- canonicalize_glycan_graph(input[[i]])
    expect_identical(result$iupac[[i]], graph_to_iupac(expected))
    expect_equal(
      igraph::as_data_frame(result$graphs[[i]], what = "both"),
      igraph::as_data_frame(expected, what = "both")
    )
  }
  trusted <- canonicalize_glycan_graphs(input, validate = FALSE)
  expect_identical(trusted$iupac, result$iupac)
  expect_equal(
    lapply(trusted$graphs[c(1, 3)], igraph::as_data_frame, what = "both"),
    lapply(result$graphs[c(1, 3)], igraph::as_data_frame, what = "both")
  )
})

test_that("batch graph canonicalization preserves floating graph attributes", {
  g <- .parse_iupac_condensed_single("{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc")
  g$source <- list(format = "WURCS")
  g <- igraph::set_vertex_attr(
    g,
    "source_id",
    value = seq_len(igraph::vcount(g))
  )
  result <- canonicalize_glycan_graphs(list(g))
  expected <- canonicalize_glycan_graph(validate_glycan_graph(g))
  expect_identical(result$iupac[[1]], graph_to_iupac(expected))
  expect_equal(
    igraph::as_data_frame(result$graphs[[1]], what = "both"),
    igraph::as_data_frame(expected, what = "both")
  )
  expect_identical(result$graphs[[1]]$source, g$source)
  expect_equal(result$graphs[[1]]$floating_parts, expected$floating_parts)
})

test_that("batch graph recovery distinguishes missing and invalid inputs", {
  expect_snapshot(
    result <- canonicalize_glycan_graphs(
      list(absent = NULL, bad = 1),
      on_failure = "na"
    )
  )
  expect_identical(result$status, c(absent = "missing", bad = "invalid"))
  expect_identical(result$iupac, c(absent = NA_character_, bad = NA_character_))
  expect_identical(names(result$graphs), c("absent", "bad"))
  expect_snapshot(error = TRUE, canonicalize_glycan_graphs(list(NULL, 1)))
  expect_identical(
    canonicalize_glycan_graphs(list()),
    list(
      graphs = list(),
      iupac = character(),
      status = character(),
      reason = character()
    )
  )
})

test_that("strict graph errors carry original position and input name", {
  cnd <- tryCatch(
    canonicalize_glycan_graphs(list(absent = NULL, bad = 1)),
    glyrepr_error_structure_failure = identity
  )
  expect_s3_class(cnd, "glyrepr_error_structure_failure")
  expect_identical(cnd$position, 2L)
  expect_identical(cnd$input_name, "bad")
  expect_type(cnd$reason, "character")
})
