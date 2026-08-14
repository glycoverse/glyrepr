test_that("floating substituents parse and serialize", {
  input <- c(
    unrestricted = "{6S}Gal(a1-3)Gal(a1-",
    restricted = "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-",
    unknown_position = "{?S}Gal(a1-3)Glc(a1-3)Man(a1-"
  )

  glycans <- as_glycan_structure(input)
  graphs <- as.list(glycans)

  expect_identical(structure_to_iupac(glycans), input)
  expect_equal(
    graphs[[1]]$floating_substituents,
    list(list(substituent = "6S", parents = integer()))
  )
  expect_equal(
    graphs[[2]]$floating_substituents,
    list(list(substituent = "6S", parents = c(1L, 2L)))
  )
  expect_equal(
    graphs[[3]]$floating_substituents,
    list(list(substituent = "?S", parents = integer()))
  )
  expect_identical(igraph::V(graphs[[1]])$sub, c("", ""))

  explicit_all <- as_glycan_structure("{6S|1,2}Gal(a1-3)Gal(a1-")
  expect_false(unname(glycans[[1]] == explicit_all))
})


test_that("singleton floating substituents resolve onto their residue", {
  explicit <- as_glycan_structure(
    "{6S|1}Gal(a1-3)Glc(a1-"
  )
  unrestricted <- as_glycan_structure("{?S}Gal(a1-")

  expect_identical(as.character(explicit), "Gal6S(a1-3)Glc(a1-")
  expect_identical(as.character(unrestricted), "Gal?S(a1-")
  expect_false(has_floating_substituents(explicit))
  expect_false(has_floating_substituents(unrestricted))
  expect_null(
    igraph::graph_attr(
      get_structure_graphs(explicit),
      "floating_substituents"
    )
  )
})


test_that("floating substituents canonicalize tokens, parents, and order", {
  input <- paste0(
    "{?S}",
    "{6/4Ac|2,1}",
    "{3Me|2,1}",
    "Gal(a1-?)Glc(a1-?)Man(a1-"
  )
  expected <- paste0(
    "{3Me|1,2}",
    "{4/6Ac|1,2}",
    "{?S}",
    "Gal(a1-?)Glc(a1-?)Man(a1-"
  )

  glycan <- as_glycan_structure(input)
  graph <- get_structure_graphs(glycan)

  expect_identical(as.character(glycan), expected)
  expect_identical(
    purrr::map_chr(graph$floating_substituents, "substituent"),
    c("3Me", "4/6Ac", "?S")
  )
  expect_identical(
    purrr::map(graph$floating_substituents, "parents"),
    list(c(1L, 2L), c(1L, 2L), integer())
  )
})


test_that("floating substituents preserve multiplicity", {
  glycan <- as_glycan_structure(
    "{6S}{6S}Gal(a1-3)Glc(a1-"
  )
  graph <- get_structure_graphs(glycan)

  expect_identical(
    as.character(glycan),
    "{6S}{6S}Gal(a1-3)Glc(a1-"
  )
  expect_length(graph$floating_substituents, 2)
  expect_identical(count_mono(glycan, "S"), 2L)
})


test_that("floating substituent assignments must be conflict-free", {
  valid <- as_glycan_structure(
    "{4/6S|1}{6Me|1}Gal(a1-"
  )
  expect_identical(as.character(valid), "Gal4/6S6Me(a1-")

  invalid <- c(
    "{6S|1}Gal6Me(a1-",
    "{6S|1,2}Gal6Me(a1-3)Glc(a1-",
    "{6S}{6Ac}Gal(a1-",
    "{4/6S|1}{4Ac|1}{6Me|1}Gal(a1-"
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


test_that("malformed floating substituent IUPAC is rejected", {
  main <- "Gal(a1-3)Glc(a1-"
  invalid <- c(
    paste0("{S}", main),
    paste0("{0S}", main),
    paste0("{6X}", main),
    paste0("{6S|}", main),
    paste0("{6S|0}", main),
    paste0("{6S|1,1}", main),
    paste0("{6S|3}", main),
    paste0("{6S|1|2}", main)
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


test_that("invalid floating substituent graph metadata is rejected", {
  graph <- get_structure_graphs(
    as_glycan_structure("Gal(a1-3)Glc(a1-")
  )
  graph$floating_substituents <- list(
    list(substituent = "6S", parents = c(1L, 2L))
  )

  cases <- list()
  cases$malformed_attr <- graph
  cases$malformed_attr$floating_substituents <- "not a list"

  cases$missing_field <- graph
  cases$missing_field$floating_substituents[[1]]$parents <- NULL

  cases$unknown_substituent <- graph
  cases$unknown_substituent$floating_substituents[[1]]$substituent <- "6X"

  cases$noncanonical_substituent <- graph
  cases$noncanonical_substituent$floating_substituents[[
    1
  ]]$substituent <- "6/4S"

  cases$duplicate_parents <- graph
  cases$duplicate_parents$floating_substituents[[1]]$parents <- c(1L, 1L)

  cases$invalid_parent <- graph
  cases$invalid_parent$floating_substituents[[1]]$parents <- 3L

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


test_that("floating substituents and floating parts coexist", {
  glycan <- as_glycan_structure(paste0(
    "{Neu5Ac(a2-6)|2,3}",
    "{8S|2,3}",
    "Gal(a1-3)Glc(a1-"
  ))
  graph <- get_structure_graphs(glycan)

  expect_identical(
    as.character(glycan),
    paste0(
      "{8S|2,3}",
      "{Neu5Ac(a2-6)|2,3}",
      "Gal(a1-3)Glc(a1-"
    )
  )
  expect_identical(
    graph$floating_substituents,
    list(list(substituent = "8S", parents = c(2L, 3L)))
  )
  expect_identical(
    graph$floating_parts[[1]]$parents,
    c(2L, 3L)
  )
})


test_that("resolving a floating part preserves substituent candidates", {
  glycan <- as_glycan_structure(paste0(
    "{6S}",
    "{Fuc(a1-2)|2}",
    "Gal(a1-3)Glc(a1-"
  ))

  expect_identical(
    as.character(glycan),
    "{6S}Fuc(a1-2)Gal(a1-3)Glc(a1-"
  )
})

test_that("floating substituents can target floating-part nodes", {
  glycan <- as_glycan_structure(paste0(
    "{6S|1}",
    "{Gal(b1-4)|2,3}",
    "Glc(a1-3)Man(a1-"
  ))

  expect_identical(
    as.character(glycan),
    "{Gal6S(b1-4)|2,3}Glc(a1-3)Man(a1-"
  )
  expect_null(get_structure_graphs(glycan)$floating_substituents)
})


test_that("floating substituents affect symmetric canonicalization", {
  make_graph <- function(edge_order, parents) {
    graph <- igraph::make_empty_graph(3, directed = TRUE)
    graph <- igraph::add_edges(graph, edge_order)
    igraph::V(graph)$name <- as.character(seq_len(3))
    igraph::V(graph)$mono <- c("Glc", "Gal", "Gal")
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- c("a1-?", "a1-?")
    graph$anomer <- "a1"
    graph$floating_substituents <- list(
      list(substituent = "6S", parents = parents)
    )
    graph
  }

  graph_a <- make_graph(c(1, 2, 1, 3), c(1L, 2L))
  graph_b <- make_graph(c(1, 3, 1, 2), c(1L, 2L))
  glycan_a <- glycan_structure(graph_a)
  glycan_b <- glycan_structure(graph_b)

  expect_identical(
    as.character(glycan_a),
    "{6S|1,3}Gal(a1-?)[Gal(a1-?)]Glc(a1-"
  )
  expect_identical(as.character(glycan_a), as.character(glycan_b))
  expect_true(unname(glycan_a == glycan_b))
})


test_that("floating substituents are included in downstream operations", {
  glycan <- as_glycan_structure(
    "{6S}{?Me}Gal3Ac(a1-3)Glc(a1-"
  )
  composition <- as_glycan_composition(glycan)

  expect_identical(count_mono(glycan, "S"), 1L)
  expect_identical(count_mono(glycan, "Me"), 1L)
  expect_identical(count_mono(glycan, "Ac"), 1L)
  expect_identical(count_mono(glycan, include_subs = TRUE), 5L)
  expect_identical(count_mono(composition, "S"), 1L)

  clean <- remove_substituents(glycan)
  clean_graph <- get_structure_graphs(clean)
  expect_identical(as.character(clean), "Gal(a1-3)Glc(a1-")
  expect_identical(igraph::V(clean_graph)$sub, c("", ""))
  expect_null(igraph::graph_attr(clean_graph, "floating_substituents"))
})


test_that("structure transformations preserve floating substituents", {
  glycan <- as_glycan_structure(
    "{?S|2,3}Fuc(a1-2)[Gal(a1-3)]Man(a1-"
  )
  strip_linkages <- function(graph) {
    graph <- igraph::set_edge_attr(graph, "linkage", value = "??-?")
    igraph::set_graph_attr(graph, "anomer", value = "??")
  }

  mapped <- smap_structure(glycan, strip_linkages)
  topological <- remove_linkages(glycan)
  generic <- convert_to_generic(glycan)

  expect_identical(
    as.character(mapped),
    "{?S|1,3}Gal(??-?)[Fuc(??-?)]Man(??-"
  )
  expect_identical(as.character(topological), as.character(mapped))
  expect_identical(
    get_structure_graphs(mapped)$floating_substituents[[1]]$parents,
    c(1L, 3L)
  )
  expect_identical(
    as.character(generic),
    "{?S|2,3}dHex(a1-2)[Hex(a1-3)]Hex(a1-"
  )
  expect_identical(count_mono(generic, "S"), 1L)
})


test_that("has_floating_substituents is vectorized", {
  glycans <- as_glycan_structure(c(
    floating = "{6S}Gal(a1-3)Glc(a1-",
    ordinary = "Gal6S(a1-",
    missing = NA
  ))

  expect_identical(
    has_floating_substituents(glycans),
    c(floating = TRUE, ordinary = FALSE, missing = NA)
  )
  expect_identical(
    has_floating_substituents(get_structure_graphs(glycans[[1]])),
    TRUE
  )
})
