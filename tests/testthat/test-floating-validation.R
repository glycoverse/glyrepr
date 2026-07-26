make_floating_validation_graph <- function(
  main_count,
  main_edges = integer(),
  main_linkages = character(),
  floating_linkages = character(),
  floating_parents = vector("list", length(floating_linkages))
) {
  graph <- igraph::make_empty_graph(
    main_count + length(floating_linkages),
    directed = TRUE
  )
  if (length(main_edges) > 0) {
    graph <- igraph::add_edges(graph, main_edges)
  }

  igraph::V(graph)$mono <- c(
    rep("Gal", main_count),
    rep("Neu5Ac", length(floating_linkages))
  )
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- main_linkages
  graph$anomer <- "a1"
  graph$floating_parts <- lapply(
    seq_along(floating_linkages),
    function(part_id) {
      list(
        root = as.integer(main_count + part_id),
        linkage = floating_linkages[[part_id]],
        parents = floating_parents[[part_id]]
      )
    }
  )
  graph
}


test_that("floating metadata rejects unsupported fields", {
  graph <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = "a2-3",
    floating_parents = list(1L)
  )
  graph$floating_parts[[1]]$probability <- 0.5

  expect_snapshot(
    error = TRUE,
    glycan_structure(graph)
  )
})

test_that("floating metadata rejects duplicated fields", {
  graph <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = "a2-3",
    floating_parents = list(1L)
  )
  names(graph$floating_parts[[1]])[[3]] <- "root"

  expect_snapshot(
    error = TRUE,
    glycan_structure(graph)
  )
})

test_that("floating metadata requires every core field", {
  graph <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = "a2-3",
    floating_parents = list(1L)
  )
  graph$floating_parts[[1]]$parents <- NULL

  expect_snapshot(
    error = TRUE,
    glycan_structure(graph)
  )
})

test_that("floating node metadata must describe exactly one component", {
  graph <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = "a2-3",
    floating_parents = list(integer())
  )

  duplicate_nodes <- graph
  duplicate_nodes$floating_parts[[1]]$nodes <- c(2L, 2L)
  expect_snapshot(
    error = TRUE,
    glycan_structure(duplicate_nodes)
  )

  missing_root <- graph
  missing_root$floating_parts[[1]]$nodes <- 1L
  expect_snapshot(
    error = TRUE,
    glycan_structure(missing_root)
  )

  wrong_component <- graph
  wrong_component$floating_parts[[1]]$nodes <- c(1L, 2L)
  expect_snapshot(
    error = TRUE,
    glycan_structure(wrong_component)
  )
})


test_that("explicit parents reject definitely occupied acceptor slots", {
  graph <- make_floating_validation_graph(
    main_count = 2,
    main_edges = c(1, 2),
    main_linkages = "b1-3",
    floating_linkages = "a2-3",
    floating_parents = list(1L)
  )

  expect_snapshot(
    error = TRUE,
    glycan_structure(graph)
  )

  alternative_position <- graph
  alternative_position$floating_parts[[1]]$linkage <- "a2-3/6"
  expect_s3_class(
    glycan_structure(alternative_position),
    "glyrepr_structure"
  )

  unknown_occupied_position <- graph
  igraph::E(unknown_occupied_position)$linkage <- "b1-?"
  expect_s3_class(
    glycan_structure(unknown_occupied_position),
    "glyrepr_structure"
  )
})


test_that("floating parts require a simultaneous attachment assignment", {
  conflict <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = c("a2-3", "a2-3"),
    floating_parents = list(1L, 1L)
  )

  expect_snapshot(
    error = TRUE,
    glycan_structure(conflict)
  )

  alternatives <- conflict
  alternatives$floating_parts[[1]]$linkage <- "a2-3/6"
  expect_s3_class(
    glycan_structure(alternatives),
    "glyrepr_structure"
  )
})


test_that("ambiguous main edges participate in slot matching", {
  one_floating <- make_floating_validation_graph(
    main_count = 2,
    main_edges = c(1, 2),
    main_linkages = "b1-3/6",
    floating_linkages = "a2-3",
    floating_parents = list(1L)
  )
  expect_s3_class(
    glycan_structure(one_floating),
    "glyrepr_structure"
  )

  conflict <- make_floating_validation_graph(
    main_count = 2,
    main_edges = c(1, 2),
    main_linkages = "b1-3/6",
    floating_linkages = c("a2-3", "a2-6"),
    floating_parents = list(1L, 1L)
  )
  expect_snapshot(
    error = TRUE,
    glycan_structure(conflict)
  )
})


test_that("unrestricted parents use all available main-tree slots", {
  graph <- make_floating_validation_graph(
    main_count = 2,
    main_edges = c(1, 2),
    main_linkages = "b1-4",
    floating_linkages = c("a2-3", "a2-3"),
    floating_parents = list(integer(), integer())
  )

  expect_s3_class(
    glycan_structure(graph),
    "glyrepr_structure"
  )

  conflict <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = c("a2-3", "a2-3"),
    floating_parents = list(integer(), integer())
  )
  expect_snapshot(
    error = TRUE,
    glycan_structure(conflict)
  )
})


test_that("unknown floating acceptor positions remain permissive", {
  graph <- make_floating_validation_graph(
    main_count = 1,
    floating_linkages = c("a2-?", "a2-?"),
    floating_parents = list(1L, 1L)
  )

  expect_s3_class(
    glycan_structure(graph),
    "glyrepr_structure"
  )
})
