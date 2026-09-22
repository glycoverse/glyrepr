test_that("array construction preserves names, missing values and duplicates", {
  a <- list(
    mono = c("Glc", "Gal"),
    sub = c("", ""),
    edges = c(1L, 2L),
    linkage = "b1-4",
    anomer = "?1"
  )
  result <- structure_from_arrays(list(a = a, missing = NULL, duplicate = a))
  expect_identical(
    as.character(result),
    c(
      a = "Gal(b1-4)Glc(?1-",
      missing = NA_character_,
      duplicate = "Gal(b1-4)Glc(?1-"
    )
  )
  expect_length(attr(result, "graphs"), 1)
  expect_identical(structure_from_arrays(list()), glycan_structure())
})

test_that("arrays preserve floating metadata and reducing end state", {
  a <- list(
    mono = c("Glc", "Gal", "Neu5Ac"),
    sub = rep("", 3),
    edges = c(1L, 2L),
    linkage = "b1-4",
    anomer = "?1",
    alditol = TRUE,
    floating_parts = list(list(
      root = 3L,
      nodes = 3L,
      linkage = "a2-6",
      parents = c(1L, 2L)
    )),
    floating_substituents = list(list(substituent = "3S", parents = c(1L, 2L)))
  )
  expected <- as_glycan_structure(
    "{3S|2,3}{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc-ol(?1-"
  )
  result <- structure_from_arrays(list(a))
  expect_identical(as.character(result), as.character(expected))
  expect_equal(
    structure_floating_parts(result),
    structure_floating_parts(expected)
  )
  expect_equal(
    structure_floating_substituents(result),
    structure_floating_substituents(expected)
  )
})

test_that("malformed array records fail before graph construction", {
  expect_snapshot(error = TRUE, structure_from_arrays(list(list(mono = "Gal"))))
  a <- list(
    mono = "Gal",
    sub = "",
    edges = c(1, 2),
    linkage = "b1-4",
    anomer = "?1"
  )
  expect_snapshot(error = TRUE, structure_from_arrays(list(a)))
  expect_snapshot(
    result <- structure_from_arrays(
      list(bad = a, absent = NULL),
      on_failure = "na"
    )
  )
  expect_identical(is.na(result), c(bad = TRUE, absent = TRUE))
})

test_that("native arrays agree with the graph reference for arbitrary node order", {
  strings <- c(
    "Gal(b1-4)[Fuc(a1-3)]GlcNAc",
    "Gal6S(b1-4)GlcNAc-ol",
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-",
    "{6S|1,2}Gal(b1-4)Glc",
    "{Neu5Ac(a2-6)|2}Gal",
    "{6S|1}Gal",
    "D-Fucf(a1-3)GlcNAc",
    "Hex(??-?)HexNAc"
  )
  graphs <- lapply(strings, .parse_iupac_condensed_single)
  records <- lapply(graphs, function(g) {
    list(
      mono = igraph::V(g)$mono,
      sub = igraph::V(g)$sub,
      edges = as.integer(t(igraph::as_edgelist(g, names = FALSE))),
      linkage = igraph::E(g)$linkage,
      anomer = g$anomer,
      alditol = isTRUE(g$alditol),
      floating_parts = g$floating_parts,
      floating_substituents = g$floating_substituents
    )
  })
  # Reverse all source node IDs, including cross-component parent domains.
  reversed <- lapply(records, function(a) {
    n <- length(a$mono)
    a$mono <- rev(a$mono)
    a$sub <- rev(a$sub)
    a$edges <- n + 1L - a$edges
    a$floating_parts <- lapply(a$floating_parts, function(p) {
      p$root <- n + 1L - p$root
      p$nodes <- n + 1L - p$nodes
      p$parents <- n + 1L - p$parents
      p
    })
    a$floating_substituents <- lapply(a$floating_substituents, function(p) {
      p$parents <- n + 1L - p$parents
      p
    })
    a
  })
  expected <- as_glycan_structure(graphs)
  # Successful native records must not need the reference graph pipeline.
  testthat::local_mocked_bindings(process_glycan_structure_element = function(
    ...
  ) {
    stop("unexpected graph fallback")
  })
  for (input in list(records, reversed)) {
    result <- structure_from_arrays(input)
    expect_identical(as.character(result), as.character(expected))
    for (accessor in list(
      structure_nodes,
      structure_edges,
      structure_floating_parts,
      structure_floating_substituents,
      get_alditol,
      get_anomer
    )) {
      expect_equal(accessor(result), accessor(expected))
    }
  }
})

test_that("native topology validation rejects cycles and incomplete components", {
  a <- list(
    mono = c("Gal", "Glc"),
    sub = c("", ""),
    edges = c(1L, 2L, 2L, 1L),
    linkage = c("b1-4", "b1-3"),
    anomer = "?1"
  )
  bad <- list(a)
  a$edges <- integer()
  a$linkage <- character()
  bad <- c(bad, list(a))
  a$floating_parts <- list(list(
    root = 2L,
    nodes = c(1L, 2L),
    linkage = "b1-4",
    parents = 1L
  ))
  bad <- c(bad, list(a))
  native <- .compact_arrays_native(
    lapply(bad, .validate_structure_arrays),
    base::order,
    .compact_bliss_labels
  )
  expect_identical(
    vapply(native, `[[`, character(1), "status"),
    rep("invalid", 3)
  )
  expect_snapshot(result <- structure_from_arrays(bad, on_failure = "na"))
  expect_identical(unname(is.na(result)), rep(TRUE, 3))
})

test_that("array failure conditions expose original positions and names", {
  good <- list(
    mono = "Gal",
    sub = "",
    edges = integer(),
    linkage = character(),
    anomer = "?1"
  )
  bad <- good
  bad$mono <- "not-a-residue"
  input <- list(
    good = good,
    bad = bad,
    absent = NULL,
    repeated = bad,
    last = good
  )
  cnd <- tryCatch(
    structure_from_arrays(input),
    glyrepr_error_structure_failure = identity
  )
  expect_s3_class(cnd, "glyrepr_error_structure_failure")
  expect_identical(cnd$position, 2L)
  expect_identical(cnd$input_name, "bad")
  expect_identical(cnd$reason, "Unknown monosaccharide.")
  warning <- NULL
  result <- withCallingHandlers(
    structure_from_arrays(input, on_failure = "na"),
    glyrepr_warning_structure_failure = function(cnd) {
      warning <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(warning$positions, c(2L, 4L))
  expect_identical(warning$reasons, rep("Unknown monosaccharide.", 2))
  expect_identical(
    as.character(result),
    c(
      good = "Gal(?1-",
      bad = NA_character_,
      absent = NA_character_,
      repeated = NA_character_,
      last = "Gal(?1-"
    )
  )
  expect_length(attr(result, "graphs"), 1)
})

test_that("native size guards use reference recovery without rejecting wide trees", {
  n <- 2049L
  a <- list(
    mono = c("Glc", rep("Gal", n - 1L)),
    sub = rep("", n),
    edges = as.integer(rbind(1L, 2:n)),
    linkage = rep("b1-?", n - 1L),
    anomer = "?1"
  )
  native <- .compact_arrays_native(
    list(.validate_structure_arrays(a)),
    base::order,
    .compact_bliss_labels
  )
  expect_identical(native[[1]]$status, "unsupported")
  result <- structure_from_arrays(list(a))
  expected <- as_glycan_structure(.compact_build_graph(a))
  expect_identical(as.character(result), as.character(expected))
  expect_equal(structure_nodes(result), structure_nodes(expected))
  expect_equal(structure_edges(result), structure_edges(expected))
})

test_that("chemical and metadata errors are isolated within array batches", {
  good <- list(
    mono = c("Glc", "Gal"),
    sub = c("", ""),
    edges = c(1L, 2L),
    linkage = "b1-4",
    anomer = "?1"
  )
  bad <- lapply(
    c("mono", "sub", "linkage", "anomer", "alditol", "edges"),
    function(field) {
      a <- good
      a[[field]] <- switch(
        field,
        mono = c("Glc", NA_character_),
        sub = c("", "6/4S"),
        linkage = "a9-4",
        anomer = "x1",
        alditol = NA,
        edges = c(1, 1.5)
      )
      a
    }
  )
  bad[[7]] <- good
  bad[[7]]$floating_substituents <- list(list(
    substituent = "6S",
    parents = c(1L, 1L)
  ))
  outcomes <- .structure_array_outcomes(c(list(good), bad, list(NULL, good)))
  expect_identical(
    vapply(outcomes, inherits, logical(1), "error"),
    c(FALSE, rep(TRUE, 7), FALSE, FALSE)
  )
  expect_identical(outcomes[[1]]$iupac, outcomes[[10]]$iupac)
})

test_that("empty and wholly missing array batches retain names without warnings", {
  expect_identical(
    structure_from_arrays(list(a = NULL, b = NULL)),
    as_glycan_structure(c(a = NA_character_, b = NA_character_))
  )
  empty <- stats::setNames(list(), character())
  expect_identical(names(structure_from_arrays(empty)), character())
})
