test_that("localize_floating_parts attaches selected parts", {
  glycans <- as_glycan_structure(c(
    floating = paste0(
      "{Fuc(a1-2)|3,4}",
      "{Neu5Ac(a2-6)|3,4}",
      "Gal(b1-3)GalNAc(a1-"
    ),
    missing = NA,
    ordinary = "Gal(b1-4)GlcNAc(b1-"
  ))
  assignments <- tibble::tibble(
    glycan_id = 1L,
    part_id = 2L,
    parent_node = 3L
  )

  localized <- localize_floating_parts(glycans, assignments)

  expect_s3_class(localized, "glyrepr_structure")
  expect_identical(names(localized), names(glycans))
  expect_identical(is.na(localized), is.na(glycans))
  expect_identical(
    as.character(localized[[3]]),
    as.character(glycans[[3]])
  )
  expect_true(has_floating_parts(localized[[1]]))
  expect_identical(
    structure_floating_parts(localized[[1]])$linkage,
    "a1-2"
  )
})

test_that("localize_floating_parts localizes unrestricted parts", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-"
  )

  localized <- localize_floating_parts(
    glycan,
    tibble::tibble(
      glycan_id = 1L,
      part_id = 1L,
      parent_node = 2L
    )
  )

  expect_identical(
    as.character(localized),
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  )
  expect_false(has_floating_parts(localized))
})

test_that("localize_floating_parts remaps remaining parent indices", {
  glycan <- as_glycan_structure(
    paste0(
      "{Fuc(a1-2)|3,4}",
      "{Neu5Ac(a2-6)|3,4}",
      "Gal(b1-3)GalNAc(a1-"
    )
  )

  localized <- localize_floating_parts(
    glycan,
    tibble::tibble(
      glycan_id = 1L,
      part_id = 2L,
      parent_node = 3L
    )
  )

  expect_identical(
    as.character(localized),
    "{Fuc(a1-2)|3,4}Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"
  )
  expect_identical(
    structure_floating_parts(localized)$parents,
    list(c(3L, 4L))
  )
})

test_that("localize_floating_parts preserves graph vertex IDs", {
  structure <- as_glycan_structure(
    paste0(
      "{Fuc(a1-2)|3,4}",
      "{Neu5Ac(a2-6)|3,4}",
      "Gal(b1-3)GalNAc(a1-"
    )
  )
  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  edges_before <- igraph::as_edgelist(graph, names = FALSE)

  result <- localize_floating_parts(
    graph,
    tibble::tibble(
      glycan_id = 1L,
      part_id = 2L,
      parent_node = 3L
    )
  )

  expect_s3_class(result, "igraph")
  expect_identical(igraph::V(result)$name, names_before)
  expect_identical(
    igraph::as_edgelist(result, names = FALSE),
    rbind(edges_before, c(3, 2))
  )
  expect_identical(
    result$floating_parts,
    list(list(
      root = 1L,
      nodes = 1L,
      linkage = "a1-2",
      parents = c(3L, 4L)
    ))
  )
})

test_that("graph localization keeps unrestricted whole-graph domains", {
  structure <- as_glycan_structure(
    paste0(
      "{Fuc(a1-2)|2,3}",
      "{Neu5Ac(a2-6)}",
      "Gal(b1-3)GalNAc(a1-"
    )
  )
  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  edges_before <- igraph::as_edgelist(graph, names = FALSE)

  result <- localize_floating_parts(
    graph,
    tibble::tibble(
      glycan_id = 1L,
      part_id = 1L,
      parent_node = 3L
    )
  )

  expect_identical(igraph::V(result)$name, names_before)
  expect_identical(
    igraph::as_edgelist(result, names = FALSE),
    rbind(edges_before, c(3, 1))
  )
  expect_identical(
    result$floating_parts,
    list(list(
      root = 2L,
      nodes = 2L,
      linkage = "a2-6",
      parents = integer()
    ))
  )
})

test_that("part localization keeps unrestricted substituent domains", {
  structure <- as_glycan_structure(
    paste0(
      "{6S}",
      "{Fuc(a1-2)|2,3}",
      "Gal(a1-3)Glc(a1-"
    )
  )
  assignments <- tibble::tibble(
    glycan_id = 1L,
    part_id = 1L,
    parent_node = 2L
  )

  localized <- localize_floating_parts(structure, assignments)

  expect_identical(
    as.character(localized),
    "{6S}Fuc(a1-2)Gal(a1-3)Glc(a1-"
  )
  expect_identical(
    structure_floating_substituents(localized)$parents,
    list(integer())
  )

  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  graph_localized <- localize_floating_parts(graph, assignments)

  expect_identical(igraph::V(graph_localized)$name, names_before)
  expect_identical(
    graph_localized$floating_substituents,
    list(list(substituent = "6S", parents = integer()))
  )
  expect_identical(
    structure_floating_candidates(graph_localized)$parent_node,
    c(1L, 2L, 3L)
  )
})

test_that("localize_floating_parts returns an unchanged graph for no assignments", {
  graph <- get_structure_graphs(as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  ))
  assignments <- tibble::tibble(
    glycan_id = integer(),
    part_id = integer(),
    parent_node = integer()
  )

  expect_identical(localize_floating_parts(graph, assignments), graph)
})

test_that("localize_floating_parts returns x for no assignments", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  assignments <- tibble::tibble(
    glycan_id = integer(),
    part_id = integer(),
    parent_node = integer()
  )

  expect_identical(
    localize_floating_parts(glycan, assignments),
    glycan
  )
})

test_that("localize_floating_parts validates assignment tables", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )

  expect_snapshot(
    localize_floating_parts(glycan, tibble::tibble()),
    error = TRUE
  )
  expect_snapshot(
    localize_floating_parts(
      glycan,
      tibble::tibble(
        glycan_id = c(1L, 1L),
        part_id = c(1L, 1L),
        parent_node = c(1L, 2L)
      )
    ),
    error = TRUE
  )
  expect_snapshot(
    localize_floating_parts(
      glycan,
      tibble::tibble(
        glycan_id = 2L,
        part_id = 1L,
        parent_node = 1L
      )
    ),
    error = TRUE
  )
})

test_that("localize_floating_parts rejects invalid targets", {
  glycans <- as_glycan_structure(c(
    restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
    missing = NA,
    ordinary = "Gal(a1-"
  ))

  expect_snapshot(
    localize_floating_parts(
      glycans,
      tibble::tibble(
        glycan_id = 1L,
        part_id = 1L,
        parent_node = 1L
      )
    ),
    error = TRUE
  )
  expect_snapshot(
    localize_floating_parts(
      glycans,
      tibble::tibble(
        glycan_id = 2L,
        part_id = 1L,
        parent_node = 1L
      )
    ),
    error = TRUE
  )
  expect_snapshot(
    localize_floating_parts(
      glycans,
      tibble::tibble(
        glycan_id = 3L,
        part_id = 1L,
        parent_node = 1L
      )
    ),
    error = TRUE
  )
})

test_that("localize_floating_parts validates simultaneous slot conflicts", {
  glycan <- as_glycan_structure(
    paste0(
      "{Fuc(a1-3)|3,4}",
      "{Neu5Ac(a2-3)|3,4}",
      "Gal(b1-4)GalNAc(a1-"
    )
  )

  expect_snapshot(
    localize_floating_parts(
      glycan,
      tibble::tibble(
        glycan_id = c(1L, 1L),
        part_id = c(1L, 2L),
        parent_node = c(3L, 3L)
      )
    ),
    error = TRUE
  )
})

test_that("localize_floating_parts checks occupied main-tree slots", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-"
  )

  expect_snapshot(
    localize_floating_parts(
      glycan,
      tibble::tibble(
        glycan_id = 1L,
        part_id = 1L,
        parent_node = 3L
      )
    ),
    error = TRUE
  )
})

test_that("enumerate_floating_localizations returns every valid variant", {
  glycan <- as_glycan_structure(
    paste0(
      "{Fuc(a1-3)|3,4}",
      "{Neu5Ac(a2-3)|3,4}",
      "Gal(b1-4)GalNAc(a1-"
    )
  )

  variants <- enumerate_floating_localizations(glycan)

  expect_s3_class(variants, "tbl_df")
  expect_named(
    variants,
    c("input_id", "variant_id", "structure", "assignments")
  )
  expect_identical(variants$input_id, c(1L, 1L))
  expect_identical(variants$variant_id, c(1L, 2L))
  expect_s3_class(variants$structure, "glyrepr_structure")
  expect_false(any(has_floating_parts(variants$structure)))
  expect_identical(
    purrr::map(variants$assignments, "glycan_id"),
    list(c(1L, 1L), c(1L, 1L))
  )
  expect_true(all(
    purrr::map_lgl(
      variants$assignments,
      ~ length(unique(.x$parent_node)) == 2
    )
  ))
})

test_that("enumeration permits parents on other floating parts", {
  glycan <- as_glycan_structure(
    "{Man(a1-?)}{Man(a1-?)}Man(a1-"
  )
  graph <- get_structure_graphs(glycan)

  graph_variants <- enumerate_floating_graph_localizations(graph)
  all_variants <- enumerate_floating_localizations(
    glycan,
    deduplicate = FALSE
  )
  unique_variants <- enumerate_floating_localizations(glycan)

  expect_identical(nrow(graph_variants), 3L)
  expect_identical(nrow(all_variants), 3L)
  expect_identical(
    sort(purrr::map_chr(
      graph_variants$assignments,
      ~ paste(.x$parent_node, collapse = ",")
    )),
    c("2,3", "3,1", "3,3")
  )
  expect_identical(
    as.character(unique_variants$structure),
    c(
      "Man(a1-?)Man(a1-?)Man(a1-",
      "Man(a1-?)[Man(a1-?)]Man(a1-"
    )
  )
})

test_that("partial localization merges floating components", {
  glycan <- as_glycan_structure(
    "{Man(a1-?)}{Man(a1-?)}Man(a1-"
  )
  assignments <- tibble::tibble(
    glycan_id = 1L,
    part_id = 1L,
    parent_node = 2L
  )

  localized <- localize_floating_parts(glycan, assignments)
  graph <- localize_floating_parts(
    get_structure_graphs(glycan),
    assignments
  )

  expect_identical(
    as.character(localized),
    "Man(a1-?)Man(a1-?)Man(a1-"
  )
  expect_identical(
    graph$floating_parts,
    list(list(
      root = 2L,
      nodes = c(1L, 2L),
      linkage = "a1-?",
      parents = integer()
    ))
  )
})

test_that("enumerate_floating_localizations handles ambiguous linkages", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-3/6)|2,3}Gal(b1-4)GalNAc(a1-"
  )

  variants <- enumerate_floating_localizations(glycan)

  expect_equal(nrow(variants), 2)
  expect_true(all(
    purrr::map_lgl(
      as.list(variants$structure),
      ~ "a2-3/6" %in% igraph::edge_attr(.x, "linkage")
    )
  ))
})

test_that("enumerate_floating_localizations filters occupied slots", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-"
  )

  variants <- enumerate_floating_localizations(glycan)

  expect_equal(nrow(variants), 1)
  expect_identical(
    as.character(variants$structure),
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  )
  expect_identical(variants$assignments[[1]]$parent_node, 2L)
})

test_that("enumerate_floating_localizations deduplicates canonical variants", {
  glycan <- as_glycan_structure(
    paste0(
      "{Neu5Ac(a2-3)|2,3}",
      "Gal(??-?)[Gal(??-?)]GlcNAc(??-"
    )
  )

  variants <- enumerate_floating_localizations(glycan)

  expect_equal(nrow(variants), 1)
  expect_identical(variants$variant_id, 1L)
  expect_identical(variants$assignments[[1]]$parent_node, 2L)
})

test_that("enumerate_floating_localizations can retain assignment provenance", {
  glycan <- as_glycan_structure(
    paste0(
      "{Neu5Ac(a2-3)|2,3}",
      "Gal(??-?)[Gal(??-?)]GlcNAc(??-"
    )
  )

  variants <- enumerate_floating_localizations(
    glycan,
    deduplicate = FALSE
  )

  expect_identical(variants$variant_id, c(1L, 2L))
  expect_identical(
    as.character(variants$structure),
    rep(as.character(variants$structure[[1]]), 2L)
  )
  expect_identical(
    purrr::map_int(variants$assignments, ~ .x$parent_node),
    c(2L, 3L)
  )
})

test_that("enumerate_floating_localizations localizes substituents", {
  glycan <- as_glycan_structure(
    "{6S}Gal(a1-3)Gal(a1-"
  )

  variants <- enumerate_floating_localizations(
    glycan,
    deduplicate = FALSE
  )

  expect_identical(
    as.character(variants$structure),
    c("Gal6S(a1-3)Gal(a1-", "Gal(a1-3)Gal6S(a1-")
  )
  expect_false(any(has_floating_substituents(variants$structure)))
  expect_identical(
    purrr::map_int(variants$assignments, ~ .x$parent_node),
    c(1L, 2L)
  )
  expect_true(all(
    purrr::map_lgl(variants$assignments, ~ is.na(.x$part_id))
  ))
  expect_identical(
    purrr::map_int(variants$assignments, ~ .x$substituent_id),
    c(1L, 1L)
  )
})

test_that("enumerate_floating_localizations validates mixed metadata", {
  glycan <- as_glycan_structure(
    paste0(
      "{6S|2,3}",
      "{Fuc(a1-6)|2,3}",
      "Gal(a1-3)Glc(a1-"
    )
  )

  variants <- enumerate_floating_localizations(
    glycan,
    deduplicate = FALSE
  )

  expect_identical(nrow(variants), 2L)
  expect_false(any(has_floating_parts(variants$structure)))
  expect_false(any(has_floating_substituents(variants$structure)))
  expect_true(all(
    purrr::map_lgl(
      variants$assignments,
      ~ length(unique(.x$parent_node)) == 2L
    )
  ))
  expect_true(all(
    purrr::map_lgl(
      variants$assignments,
      ~ identical(.x$part_id, c(1L, NA_integer_)) &&
        identical(.x$substituent_id, c(NA_integer_, 1L))
    )
  ))
})

test_that("graph localization materializes floating substituents", {
  graph <- get_structure_graphs(
    as_glycan_structure("{?S|1,2}Gal(a1-3)Glc(a1-"),
    return_list = FALSE
  )

  variants <- enumerate_floating_graph_localizations(graph)

  expect_identical(variants$variant_id, c(1L, 2L))
  expect_identical(
    purrr::map_int(variants$assignments, ~ .x$parent_node),
    c(1L, 2L)
  )
  expect_true(all(
    purrr::map_lgl(
      variants$graph,
      ~ identical(igraph::V(.x)$name, igraph::V(graph)$name)
    )
  ))
  expect_true(all(
    purrr::map_lgl(
      variants$graph,
      ~ !("floating_substituents" %in% igraph::graph_attr_names(.x))
    )
  ))
  expect_identical(
    purrr::map(variants$graph, ~ igraph::V(.x)$sub),
    list(c("?S", ""), c("", "?S"))
  )
})

test_that("graph localizations preserve original vertex IDs", {
  glycan <- as_glycan_structure(
    paste0(
      "{Neu5Ac(a2-3)|2,3}",
      "Gal(??-?)[Gal(??-?)]GlcNAc(??-"
    )
  )
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  original_edges <- igraph::as_edgelist(graph, names = FALSE)

  variants <- enumerate_floating_graph_localizations(graph)

  expect_s3_class(variants, "tbl_df")
  expect_named(variants, c("variant_id", "graph", "assignments"))
  expect_identical(variants$variant_id, c(1L, 2L))
  expect_identical(
    purrr::map_int(variants$assignments, ~ .x$parent_node),
    c(2L, 3L)
  )
  expect_identical(
    purrr::map_lgl(variants$graph, igraph::is_igraph),
    c(TRUE, TRUE)
  )
  expect_identical(
    purrr::map_lgl(
      variants$graph,
      ~ identical(igraph::V(.x)$name, igraph::V(graph)$name)
    ),
    c(TRUE, TRUE)
  )
  expect_identical(
    purrr::map_lgl(
      variants$graph,
      ~ identical(
        igraph::vertex_attr(.x),
        igraph::vertex_attr(graph)
      )
    ),
    c(TRUE, TRUE)
  )
  expect_identical(
    purrr::map_lgl(
      variants$graph,
      ~ identical(
        igraph::as_edgelist(.x, names = FALSE)[
          seq_len(nrow(original_edges)),
          ,
          drop = FALSE
        ],
        original_edges
      )
    ),
    c(TRUE, TRUE)
  )
  expect_identical(
    purrr::map2(
      variants$graph,
      variants$assignments,
      ~ unname(igraph::as_edgelist(.x, names = FALSE)[
        igraph::ecount(.x),
        ,
        drop = TRUE
      ])
    ),
    purrr::map(
      variants$assignments,
      ~ as.numeric(c(.x$parent_node, 1L))
    )
  )
  expect_identical(
    purrr::map_lgl(
      variants$graph,
      ~ !("floating_parts" %in% igraph::graph_attr_names(.x))
    ),
    c(TRUE, TRUE)
  )
})

test_that("graph localization returns an ordinary graph unchanged", {
  graph <- get_structure_graphs(
    as_glycan_structure("Gal(b1-4)GlcNAc(b1-"),
    return_list = FALSE
  )

  variants <- enumerate_floating_graph_localizations(graph)

  expect_identical(variants$variant_id, 1L)
  expect_identical(variants$graph, list(graph))
  expect_identical(variants$assignments, list(empty_floating_assignments()))
})

test_that("graph localization validates inputs and its conservative bound", {
  graph <- get_structure_graphs(
    as_glycan_structure(
      paste0(
        "{Fuc(a1-3)|3,4}",
        "{Neu5Ac(a2-3)|3,4}",
        "Gal(b1-4)GalNAc(a1-"
      )
    ),
    return_list = FALSE
  )

  expect_snapshot(
    enumerate_floating_graph_localizations(graph, max_variants = 3),
    error = TRUE
  )
  expect_error(
    enumerate_floating_graph_localizations("not a graph"),
    class = "error"
  )
})

test_that("enumerate_floating_localizations retains every input position", {
  glycans <- as_glycan_structure(c(
    missing = NA,
    ordinary = "Gal(a1-",
    floating = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  ))

  variants <- enumerate_floating_localizations(glycans)

  expect_identical(variants$input_id, c(1L, 2L, 3L, 3L))
  expect_identical(variants$variant_id, c(1L, 1L, 1L, 2L))
  expect_identical(
    unname(is.na(variants$structure)),
    c(TRUE, FALSE, FALSE, FALSE)
  )
  expect_identical(
    purrr::map_int(variants$assignments, nrow),
    c(0L, 0L, 1L, 1L)
  )
  expect_identical(
    names(variants$structure),
    c("missing", "ordinary", "floating", "floating")
  )
})

test_that("enumerate_floating_localizations handles empty vectors", {
  variants <- enumerate_floating_localizations(glycan_structure())

  expect_named(
    variants,
    c("input_id", "variant_id", "structure", "assignments")
  )
  expect_equal(nrow(variants), 0)
  expect_s3_class(variants$structure, "glyrepr_structure")
})

test_that("enumerate_floating_localizations enforces a conservative bound", {
  glycan <- as_glycan_structure(
    paste0(
      "{Fuc(a1-3)|3,4}",
      "{Neu5Ac(a2-3)|3,4}",
      "Gal(b1-4)GalNAc(a1-"
    )
  )

  expect_snapshot(
    enumerate_floating_localizations(glycan, max_variants = 3),
    error = TRUE
  )
  expect_error(
    enumerate_floating_localizations(glycan, max_variants = 0),
    class = "error"
  )
  expect_error(
    enumerate_floating_localizations(glycan, deduplicate = NA),
    class = "error"
  )
})
