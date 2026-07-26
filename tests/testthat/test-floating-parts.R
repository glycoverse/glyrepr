test_that("floating parts parse and serialize for a bi-antennary N-glycan", {
  main <- paste0(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
    "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(a1-"
  )
  iupac <- paste0("{Neu5Ac(a2-3)}", main)

  glycan <- as_glycan_structure(iupac)
  graph <- get_structure_graphs(glycan, return_list = FALSE)

  expect_identical(structure_to_iupac(glycan), iupac)
  expect_equal(igraph::vcount(graph), 10)
  expect_equal(igraph::components(graph, mode = "weak")$no, 2)
  expect_equal(
    graph$floating_parts,
    list(list(root = 10L, linkage = "a2-3", parents = integer()))
  )

  ordinary <- as_glycan_structure(main)
  expect_identical(structure_to_iupac(ordinary), main)
  expect_null(
    igraph::graph_attr(
      get_structure_graphs(ordinary, return_list = FALSE),
      "floating_parts"
    )
  )
})


test_that("unrestricted and explicit candidate parents remain distinguishable", {
  main <- paste0(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
    "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(a1-"
  )
  unrestricted_iupac <- paste0("{Neu5Ac(a2-3)}", main)
  explicit_iupac <- paste0("{Neu5Ac(a2-3)|1,4}", main)

  unrestricted <- as_glycan_structure(unrestricted_iupac)
  explicit <- as_glycan_structure(explicit_iupac)
  unrestricted_graph <- get_structure_graphs(
    unrestricted,
    return_list = FALSE
  )
  explicit_graph <- get_structure_graphs(explicit, return_list = FALSE)

  expect_identical(
    unrestricted_graph$floating_parts[[1]]$parents,
    integer()
  )
  expect_identical(explicit_graph$floating_parts[[1]]$parents, c(1L, 4L))
  expect_identical(structure_to_iupac(unrestricted), unrestricted_iupac)
  expect_identical(structure_to_iupac(explicit), explicit_iupac)
  expect_false(unrestricted == explicit)
})

test_that("has_floating_parts is vectorized and preserves missing values", {
  glycans <- as_glycan_structure(c(
    floating = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
    ordinary = "Gal(a1-",
    missing = NA
  ))

  expect_identical(
    has_floating_parts(glycans),
    c(floating = TRUE, ordinary = FALSE, missing = NA)
  )
})


test_that("explicit candidate parents use canonical main-tree node indices", {
  main <- "Gal(b1-3)GalNAc(a1-"
  parsed <- as_glycan_structure(
    paste0("{Neu5Ac(a2-3)|2,1}", main)
  )
  parsed_graph <- get_structure_graphs(parsed, return_list = FALSE)

  expect_identical(
    structure_to_iupac(parsed),
    paste0("{Neu5Ac(a2-3)|1,2}", main)
  )
  expect_identical(parsed_graph$floating_parts[[1]]$parents, c(1L, 2L))
  expect_identical(igraph::V(parsed_graph)$name, c("1", "2", "3"))
  expect_identical(igraph::V(parsed_graph)$mono, c("Gal", "GalNAc", "Neu5Ac"))

  graph <- igraph::make_empty_graph(3, directed = TRUE)
  graph <- igraph::add_edges(graph, c(2, 3))
  igraph::V(graph)$name <- c("floating", "reducing", "gal")
  igraph::V(graph)$mono <- c("Neu5Ac", "GalNAc", "Gal")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-3"
  graph$anomer <- "a1"
  graph$floating_parts <- list(
    list(root = 1L, linkage = "a2-3", parents = 3L)
  )

  constructed <- glycan_structure(graph)
  canonical_graph <- get_structure_graphs(
    constructed,
    return_list = FALSE
  )

  expect_identical(
    structure_to_iupac(constructed),
    "{Neu5Ac(a2-3)|1}Gal(b1-3)GalNAc(a1-"
  )
  expect_equal(
    canonical_graph$floating_parts,
    list(list(root = 3L, linkage = "a2-3", parents = 1L))
  )
})


test_that("explicit candidate indices are remapped with the main structure", {
  noncanonical_main <- "Fuc(a1-6)[Gal(a1-3)]Man(a1-"

  expect_identical(
    structure_to_iupac(
      as_glycan_structure(
        paste0("{Neu5Ac(a2-3)}", noncanonical_main)
      )
    ),
    "{Neu5Ac(a2-3)}Gal(a1-3)[Fuc(a1-6)]Man(a1-"
  )
  expect_identical(
    structure_to_iupac(
      as_glycan_structure(
        paste0("{Neu5Ac(a2-3)|1}", noncanonical_main)
      )
    ),
    "{Neu5Ac(a2-3)|2}Gal(a1-3)[Fuc(a1-6)]Man(a1-"
  )
})


test_that("multi-residue floating substructures round-trip", {
  iupac <- "{Gal(b1-4)GlcNAc(b1-3)}GalNAc(a1-"

  glycan <- as_glycan_structure(iupac)
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  components <- igraph::components(graph, mode = "weak")

  expect_identical(structure_to_iupac(glycan), iupac)
  expect_equal(sort(components$csize), c(1L, 2L))
  expect_identical(igraph::V(graph)$mono, c("GalNAc", "Gal", "GlcNAc"))
  expect_equal(
    graph$floating_parts,
    list(list(root = 3L, linkage = "b1-3", parents = integer()))
  )
})


test_that("malformed floating-part IUPAC is rejected", {
  main <- "Gal(b1-3)GalNAc(a1-"
  invalid <- c(
    "{Neu5Ac(a2-3)}",
    paste0("{}", main),
    paste0("{Neu5Ac}", main),
    paste0("{Neu5Ac(a2-3)|}", main),
    paste0("{Neu5Ac(a2-3)|0}", main),
    paste0("{Neu5Ac(a2-3)|1,1}", main),
    paste0("{Neu5Ac(a2-3)|3}", main),
    paste0("{Neu5Ac(a2-3)|1|2}", main),
    paste0("{{Neu5Ac(a2-3)}}", main),
    paste0(main, "{Neu5Ac(a2-3)}")
  )

  errors <- purrr::map(invalid, function(iupac) {
    tryCatch(
      {
        as_glycan_structure(iupac)
        NULL
      },
      error = identity
    )
  })

  expect_identical(
    purrr::map_lgl(errors, inherits, what = "error"),
    rep(TRUE, length(invalid))
  )
  expect_snapshot(
    cat(
      paste0(
        invalid,
        "\n",
        purrr::map_chr(errors, conditionMessage),
        collapse = "\n\n"
      )
    )
  )
})


test_that("invalid floating-part graph annotations are rejected", {
  graph <- igraph::make_empty_graph(3, directed = TRUE)
  graph <- igraph::add_edges(graph, c(2, 3))
  igraph::V(graph)$name <- as.character(seq_len(3))
  igraph::V(graph)$mono <- c("Neu5Ac", "GalNAc", "Gal")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-3"
  graph$anomer <- "a1"
  graph$floating_parts <- list(
    list(root = 1L, linkage = "a2-3", parents = 3L)
  )

  cases <- list()
  cases$unannotated <- igraph::delete_graph_attr(graph, "floating_parts")

  cases$malformed_attr <- graph
  cases$malformed_attr$floating_parts <- "not a list"

  cases$nonroot <- graph
  cases$nonroot$floating_parts[[1]]$root <- 3L

  cases$floating_parent <- graph
  cases$floating_parent$floating_parts[[1]]$parents <- 1L

  cases$duplicate_parents <- graph
  cases$duplicate_parents$floating_parts[[1]]$parents <- c(3L, 3L)

  cases$duplicate_component <- graph
  cases$duplicate_component$floating_parts <- c(
    graph$floating_parts,
    graph$floating_parts
  )

  cases$undeclared_component <- igraph::add_vertices(
    graph,
    1,
    name = "4",
    mono = "Fuc",
    sub = ""
  )

  errors <- purrr::map(cases, function(case) {
    tryCatch(
      {
        glycan_structure(case)
        NULL
      },
      error = identity
    )
  })

  expect_identical(
    unname(purrr::map_lgl(errors, inherits, what = "error")),
    rep(TRUE, length(cases))
  )
  expect_snapshot(
    cat(
      paste0(
        names(errors),
        "\n",
        purrr::map_chr(errors, conditionMessage),
        collapse = "\n\n"
      )
    )
  )
})


test_that("floating-part identity supports equality and vctrs operations", {
  main <- "Gal(b1-3)GalNAc(a1-"
  unrestricted <- as_glycan_structure(
    paste0("{Neu5Ac(a2-3)}", main)
  )
  explicit <- as_glycan_structure(
    paste0("{Neu5Ac(a2-3)|1}", main)
  )

  combined <- vctrs::vec_c(unrestricted, unrestricted, explicit, NA)

  expect_identical(
    vctrs::vec_match(combined, vctrs::vec_c(unrestricted, explicit)),
    c(1L, 1L, 2L, NA_integer_)
  )
  expect_length(vctrs::vec_unique(combined), 3)
  expect_length(attr(combined, "graphs"), 2)
  expect_identical(
    structure_to_iupac(vctrs::vec_slice(combined, c(3L, 1L, 4L))),
    c(
      paste0("{Neu5Ac(a2-3)|1}", main),
      paste0("{Neu5Ac(a2-3)}", main),
      NA_character_
    )
  )
})


test_that("floating-part vectors preserve duplicates, missing values, and names", {
  unrestricted <- "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-"
  explicit <- "{Neu5Ac(a2-3)|1}Gal(b1-3)GalNAc(a1-"
  input <- c(
    first = unrestricted,
    missing = NA_character_,
    duplicate = unrestricted,
    restricted = explicit
  )

  result <- as_glycan_structure(input)

  expect_identical(names(result), names(input))
  expect_identical(structure_to_iupac(result), input)
  expect_identical(unname(is.na(result)), c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(
    unname(result[c(1L, 3L)]),
    unname(result[c(3L, 1L)])
  )
  expect_length(attr(result, "graphs"), 2)
})


test_that("floating components preserve multiplicity and canonical order", {
  input <- paste0(
    "{Neu5Ac(a2-3)}",
    "{Fuc(a1-3)}",
    "{Neu5Ac(a2-3)}",
    "Gal(a1-"
  )
  expected <- paste0(
    "{Fuc(a1-3)}",
    "{Neu5Ac(a2-3)}",
    "{Neu5Ac(a2-3)}",
    "Gal(a1-"
  )

  glycan <- as_glycan_structure(input)
  graph <- get_structure_graphs(glycan, return_list = FALSE)

  expect_identical(structure_to_iupac(glycan), expected)
  expect_length(graph$floating_parts, 3)
  expect_identical(
    purrr::map_int(graph$floating_parts, "root"),
    c(2L, 3L, 4L)
  )
  expect_equal(igraph::components(graph, mode = "weak")$no, 4)
  expect_identical(count_mono(glycan, "Neu5Ac"), 2L)
  expect_identical(count_mono(glycan, "Fuc"), 1L)
})
