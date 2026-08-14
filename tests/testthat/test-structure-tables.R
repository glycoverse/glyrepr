test_that("structure_nodes expands duplicated structures", {
  glycans <- c(o_glycan_core_1(), o_glycan_core_1())

  nodes <- structure_nodes(glycans)

  expect_s3_class(nodes, "tbl_df")
  expect_named(nodes, c("glycan_id", "node_id", "mono", "sub"))
  expect_equal(nodes$glycan_id, c(1L, 1L, 2L, 2L))
  expect_equal(nodes$node_id, c(1L, 2L, 1L, 2L))
  expect_equal(nodes$mono, c("Gal", "GalNAc", "Gal", "GalNAc"))
  expect_equal(nodes$sub, c("", "", "", ""))
})

test_that("structure_edges expands duplicated structures", {
  glycans <- c(o_glycan_core_1(), o_glycan_core_1())

  edges <- structure_edges(glycans)

  expect_s3_class(edges, "tbl_df")
  expect_named(
    edges,
    c(
      "glycan_id",
      "edge_id",
      "from_node",
      "to_node",
      "linkage"
    )
  )
  expect_equal(edges$glycan_id, c(1L, 2L))
  expect_equal(edges$edge_id, c(1L, 1L))
  expect_equal(edges$from_node, c(2L, 2L))
  expect_equal(edges$to_node, c(1L, 1L))
  expect_equal(edges$linkage, c("b1-3", "b1-3"))
})

test_that("structure_nodes includes glycan_name for named structures", {
  glycans <- c(
    first = o_glycan_core_1(),
    missing = NA,
    second = o_glycan_core_1()
  )

  nodes <- structure_nodes(glycans)

  expect_named(nodes, c("glycan_id", "glycan_name", "node_id", "mono", "sub"))
  expect_equal(nodes$glycan_id, c(1L, 1L, 3L, 3L))
  expect_equal(nodes$glycan_name, c("first", "first", "second", "second"))
  expect_equal(nodes$node_id, c(1L, 2L, 1L, 2L))
})

test_that("structure_edges includes glycan_name for named structures", {
  glycans <- c(
    first = o_glycan_core_1(),
    missing = NA,
    second = o_glycan_core_1()
  )

  edges <- structure_edges(glycans)

  expect_named(
    edges,
    c(
      "glycan_id",
      "glycan_name",
      "edge_id",
      "from_node",
      "to_node",
      "linkage"
    )
  )
  expect_equal(edges$glycan_id, c(1L, 3L))
  expect_equal(edges$glycan_name, c("first", "second"))
})

test_that("structure_floating_parts returns normalized attachment metadata", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
    restricted = "{Fuc(a1-2)|3,4}{Neu5Ac(a2-6)|3,4}Gal(b1-3)GalNAc(a1-",
    ordinary = "Gal(a1-",
    missing = NA
  ))

  parts <- structure_floating_parts(glycans)

  expect_named(
    parts,
    c(
      "glycan_id",
      "glycan_name",
      "part_id",
      "root_node",
      "nodes",
      "linkage",
      "parents"
    )
  )
  expect_equal(parts$glycan_id, c(1L, 2L, 2L))
  expect_equal(parts$glycan_name, c("unrestricted", "restricted", "restricted"))
  expect_equal(parts$part_id, c(1L, 1L, 2L))
  expect_equal(parts$nodes, list(1L, 1L, 2L))
  expect_equal(parts$linkage, c("a2-3", "a1-2", "a2-6"))
  expect_equal(
    parts$parents,
    list(integer(), c(3L, 4L), c(3L, 4L))
  )

  empty <- structure_floating_parts(glycan_structure())
  expect_named(
    empty,
    c(
      "glycan_id",
      "part_id",
      "root_node",
      "nodes",
      "linkage",
      "parents"
    )
  )
  expect_identical(empty$nodes, list())
})

test_that("structure_floating_substituents returns normalized metadata", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{6S}Gal(a1-3)Glc(a1-",
    restricted = "{?Me|1,2}Gal(a1-3)Glc(a1-",
    ordinary = "Gal6S(a1-",
    missing = NA
  ))

  substituents <- structure_floating_substituents(glycans)

  expect_named(
    substituents,
    c(
      "glycan_id",
      "glycan_name",
      "substituent_id",
      "substituent",
      "parents"
    )
  )
  expect_identical(substituents$glycan_id, c(1L, 2L))
  expect_identical(
    substituents$glycan_name,
    c("unrestricted", "restricted")
  )
  expect_identical(substituents$substituent_id, c(1L, 1L))
  expect_identical(substituents$substituent, c("6S", "?Me"))
  expect_identical(substituents$parents, list(integer(), c(1L, 2L)))

  empty <- structure_floating_substituents(glycan_structure())
  expect_named(
    empty,
    c("glycan_id", "substituent_id", "substituent", "parents")
  )
  expect_identical(empty$parents, list())
})

test_that("structure_floating_candidates expands every attachment candidate", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
    restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
    ordinary = "Gal(a1-",
    missing = NA
  ))

  candidates <- structure_floating_candidates(glycans)

  expect_named(
    candidates,
    c(
      "glycan_id",
      "glycan_name",
      "part_id",
      "root_node",
      "parent_node",
      "linkage",
      "scope",
      "substituent_id",
      "substituent"
    )
  )
  expect_equal(candidates$glycan_id, c(1L, 1L, 2L, 2L))
  expect_equal(
    candidates$glycan_name,
    c("unrestricted", "unrestricted", "restricted", "restricted")
  )
  expect_equal(candidates$part_id, rep(1L, 4))
  expect_equal(candidates$root_node, rep(1L, 4))
  expect_equal(candidates$parent_node, c(2L, 3L, 2L, 3L))
  expect_equal(candidates$linkage, c("a2-3", "a2-3", "a2-6", "a2-6"))
  expect_equal(candidates$scope, c("all", "all", "explicit", "explicit"))
  expect_true(all(is.na(candidates$substituent_id)))
  expect_true(all(is.na(candidates$substituent)))
})

test_that("structure_floating_candidates expands floating substituents", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{6S}Gal(a1-3)Gal(a1-",
    restricted = "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-",
    unknown_position = "{?S}Gal(a1-3)Glc(a1-3)Man(a1-"
  ))

  candidates <- structure_floating_candidates(glycans)

  expect_named(
    candidates,
    c(
      "glycan_id",
      "glycan_name",
      "part_id",
      "root_node",
      "parent_node",
      "linkage",
      "scope",
      "substituent_id",
      "substituent"
    )
  )
  expect_identical(candidates$glycan_id, c(1L, 1L, 2L, 2L, 3L, 3L, 3L))
  expect_identical(
    candidates$glycan_name,
    c(
      "unrestricted",
      "unrestricted",
      "restricted",
      "restricted",
      rep("unknown_position", 3L)
    )
  )
  expect_true(all(is.na(candidates$part_id)))
  expect_true(all(is.na(candidates$root_node)))
  expect_identical(candidates$parent_node, c(1L, 2L, 1L, 2L, 1L, 2L, 3L))
  expect_true(all(is.na(candidates$linkage)))
  expect_identical(
    candidates$scope,
    c("all", "all", "explicit", "explicit", "all", "all", "all")
  )
  expect_identical(candidates$substituent_id, rep(1L, 7L))
  expect_identical(candidates$substituent, c(rep("6S", 4L), rep("?S", 3L)))
})

test_that("structure_floating_candidates returns typed empty tables", {
  unnamed <- structure_floating_candidates(glycan_structure())
  named <- structure_floating_candidates(
    c(ordinary = as_glycan_structure("Gal(a1-"), missing = NA)
  )

  expect_named(
    unnamed,
    c(
      "glycan_id",
      "part_id",
      "root_node",
      "parent_node",
      "linkage",
      "scope",
      "substituent_id",
      "substituent"
    )
  )
  expect_equal(nrow(unnamed), 0)
  expect_named(
    named,
    c(
      "glycan_id",
      "glycan_name",
      "part_id",
      "root_node",
      "parent_node",
      "linkage",
      "scope",
      "substituent_id",
      "substituent"
    )
  )
  expect_equal(nrow(named), 0)
})

test_that("structure_component_membership identifies every graph component", {
  glycans <- as_glycan_structure(c(
    floating = paste0(
      "{Fuc(a1-2)|3,4}",
      "{Neu5Ac(a2-6)|3,4}",
      "Gal(b1-3)GalNAc(a1-"
    ),
    ordinary = "Gal(a1-",
    missing = NA
  ))

  membership <- structure_component_membership(glycans)

  expect_named(
    membership,
    c(
      "glycan_id",
      "glycan_name",
      "node_id",
      "component_type",
      "part_id"
    )
  )
  expect_equal(membership$glycan_id, c(1L, 1L, 1L, 1L, 2L))
  expect_equal(
    membership$glycan_name,
    c(rep("floating", 4), "ordinary")
  )
  expect_equal(membership$node_id, c(1L, 2L, 3L, 4L, 1L))
  expect_equal(
    membership$component_type,
    c("floating", "floating", "main", "main", "main")
  )
  expect_equal(membership$part_id, c(1L, 2L, NA, NA, NA))
})

test_that("structure_component_membership returns typed empty tables", {
  unnamed <- structure_component_membership(glycan_structure())
  named <- structure_component_membership(
    c(missing = glycan_structure(NA))
  )

  expect_named(
    unnamed,
    c("glycan_id", "node_id", "component_type", "part_id")
  )
  expect_equal(nrow(unnamed), 0)
  expect_named(
    named,
    c(
      "glycan_id",
      "glycan_name",
      "node_id",
      "component_type",
      "part_id"
    )
  )
  expect_equal(nrow(named), 0)
})

test_that("structure_candidate_edges exposes potential virtual edges", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
    restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
    ordinary = "Gal(a1-",
    missing = NA
  ))

  edges <- structure_candidate_edges(glycans)

  expect_named(
    edges,
    c(
      "glycan_id",
      "glycan_name",
      "part_id",
      "from_node",
      "to_node",
      "linkage",
      "scope"
    )
  )
  expect_equal(edges$glycan_id, c(1L, 1L, 2L, 2L))
  expect_equal(edges$part_id, rep(1L, 4))
  expect_equal(edges$from_node, c(2L, 3L, 2L, 3L))
  expect_equal(edges$to_node, rep(1L, 4))
  expect_equal(edges$linkage, c("a2-3", "a2-3", "a2-6", "a2-6"))
  expect_equal(edges$scope, c("all", "all", "explicit", "explicit"))
})

test_that("structure_candidate_edges returns a typed empty table", {
  edges <- structure_candidate_edges(glycan_structure())

  expect_named(
    edges,
    c(
      "glycan_id",
      "part_id",
      "from_node",
      "to_node",
      "linkage",
      "scope"
    )
  )
  expect_equal(nrow(edges), 0)
})

test_that("structure table accessors accept one glycan graph", {
  structure <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  graph <- get_structure_graphs(structure)
  accessors <- list(
    structure_nodes,
    structure_edges,
    structure_floating_parts,
    structure_floating_substituents,
    structure_floating_candidates,
    structure_component_membership,
    structure_candidate_edges
  )

  for (accessor in accessors) {
    expect_identical(accessor(graph), accessor(structure))
  }
})

test_that("structure tables retain current graph vertex positions", {
  graph <- get_structure_graphs(o_glycan_core_1())
  graph <- igraph::permute(graph, c(2L, 1L))

  nodes <- structure_nodes(graph)
  edges <- structure_edges(graph)
  edge_ends <- igraph::as_edgelist(graph, names = FALSE)

  expect_named(nodes, c("glycan_id", "node_id", "mono", "sub"))
  expect_identical(nodes$glycan_id, rep(1L, igraph::vcount(graph)))
  expect_identical(nodes$node_id, seq_len(igraph::vcount(graph)))
  expect_identical(nodes$mono, igraph::V(graph)$mono)
  expect_identical(edges$from_node, as.integer(edge_ends[, 1]))
  expect_identical(edges$to_node, as.integer(edge_ends[, 2]))
})

test_that("structure_from_tibbles recreates structure vectors", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)
  anomers <- get_anomer(glycans)

  rebuilt <- structure_from_tibbles(nodes, edges, anomers)

  expect_s3_class(rebuilt, "glyrepr_structure")
  expect_equal(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
})

test_that("structure tables preserve cross-component parent IDs", {
  glycan <- as_glycan_structure(
    "{Fuc(a1-2)|2,3}{Man(a1-3)|1,3}Glc(a1-"
  )
  parts <- structure_floating_parts(glycan)

  expect_identical(parts$parents, list(c(2L, 3L), c(1L, 3L)))
  expect_identical(
    structure_candidate_edges(glycan)$from_node,
    c(2L, 3L, 1L, 3L)
  )

  rebuilt <- structure_from_tibbles(
    structure_nodes(glycan),
    structure_edges(glycan),
    get_anomer(glycan),
    parts
  )

  expect_identical(as.character(rebuilt), as.character(glycan))
  expect_true(unname(rebuilt == glycan))
})

test_that("structure tables round-trip alditol status", {
  glycans <- as_glycan_structure(c(
    reduced = "Gal(b1-4)GlcNAc-ol(a1-",
    missing = NA,
    ordinary = "Gal(b1-4)GlcNAc(a1-"
  ))

  rebuilt <- structure_from_tibbles(
    structure_nodes(glycans),
    structure_edges(glycans),
    get_anomer(glycans),
    alditols = get_alditol(glycans)
  )

  expect_identical(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
  expect_identical(get_alditol(rebuilt), get_alditol(glycans))
  expect_identical(names(rebuilt), names(glycans))

  ordinary <- structure_from_tibbles(
    structure_nodes(glycans),
    structure_edges(glycans),
    get_anomer(glycans)
  )
  expect_identical(
    unname(get_alditol(ordinary)),
    c(FALSE, NA, FALSE)
  )

  without_missing <- glycans[c(1, 3)]
  all_reduced <- structure_from_tibbles(
    structure_nodes(without_missing),
    structure_edges(without_missing),
    get_anomer(without_missing),
    alditols = TRUE
  )
  expect_identical(unname(get_alditol(all_reduced)), c(TRUE, TRUE))
})

test_that("structure_from_tibbles validates alditol status", {
  glycan <- o_glycan_core_1()
  nodes <- structure_nodes(glycan)
  edges <- structure_edges(glycan)

  expect_snapshot(
    structure_from_tibbles(nodes, edges, "a1", alditols = logical()),
    error = TRUE
  )
  expect_snapshot(
    structure_from_tibbles(nodes, edges, "a1", alditols = NA),
    error = TRUE
  )
  expect_snapshot(
    structure_from_tibbles(nodes, edges, "a1", alditols = "FALSE"),
    error = TRUE
  )
})

test_that("structure tables round-trip floating parts", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
    missing = NA,
    restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
    duplicate = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-"
  ))

  rebuilt <- structure_from_tibbles(
    structure_nodes(glycans),
    structure_edges(glycans),
    get_anomer(glycans),
    structure_floating_parts(glycans)
  )

  expect_identical(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
  expect_identical(names(rebuilt), names(glycans))
  expect_identical(is.na(rebuilt), is.na(glycans))
  expect_identical(
    structure_floating_parts(rebuilt)$nodes,
    structure_floating_parts(glycans)$nodes
  )
})

test_that("structure tables round-trip floating substituents", {
  glycans <- as_glycan_structure(c(
    unrestricted = "{6S}Gal(a1-3)Glc(a1-",
    missing = NA,
    restricted = "{?Me|1,2}Gal(a1-3)Glc(a1-",
    duplicate = "{6S}Gal(a1-3)Glc(a1-"
  ))

  rebuilt <- structure_from_tibbles(
    structure_nodes(glycans),
    structure_edges(glycans),
    get_anomer(glycans),
    structure_floating_parts(glycans),
    structure_floating_substituents(glycans)
  )

  expect_identical(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
  expect_identical(names(rebuilt), names(glycans))
  expect_identical(is.na(rebuilt), is.na(glycans))
  expect_identical(
    structure_floating_substituents(rebuilt),
    structure_floating_substituents(glycans)
  )
})

test_that("structure tables resolve a single substituent candidate", {
  nodes <- tibble::tibble(
    glycan_id = 1L,
    node_id = 1L,
    mono = "Gal",
    sub = ""
  )
  edges <- empty_structure_edges()
  floating_substituents <- tibble::tibble(
    glycan_id = 1L,
    substituent_id = 1L,
    substituent = "6S",
    parents = list(1L)
  )

  result <- structure_from_tibbles(
    nodes,
    edges,
    "a1",
    floating_substituents = floating_substituents
  )

  expect_identical(as.character(result), "Gal6S(a1-")
  expect_false(has_floating_substituents(result))
  expect_equal(nrow(structure_floating_substituents(result)), 0)
})

test_that("structure tables resolve a single candidate parent", {
  nodes <- tibble::tibble(
    glycan_id = c(1L, 1L),
    node_id = c(1L, 2L),
    mono = c("Gal", "Neu5Ac"),
    sub = c("", "")
  )
  edges <- tibble::tibble(
    glycan_id = integer(),
    edge_id = integer(),
    from_node = integer(),
    to_node = integer(),
    linkage = character()
  )
  floating_parts <- tibble::tibble(
    glycan_id = 1L,
    part_id = 1L,
    root_node = 2L,
    linkage = "a2-3",
    parents = list(1L)
  )

  result <- structure_from_tibbles(
    nodes,
    edges,
    "a1",
    floating_parts
  )

  expect_identical(as.character(result), "Neu5Ac(a2-3)Gal(a1-")
  expect_false(has_floating_parts(result))
  expect_equal(nrow(structure_floating_parts(result)), 0)
  expect_identical(structure_edges(result)$linkage, "a2-3")
})

test_that("structure_from_tibbles restores names from glycan_name", {
  glycans <- c(first = o_glycan_core_1(), second = n_glycan_core())
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)
  anomers <- unname(get_anomer(glycans))

  rebuilt <- structure_from_tibbles(nodes, edges, anomers)

  expect_equal(names(rebuilt), c("first", "second"))
  expect_equal(
    unname(structure_to_iupac(rebuilt)),
    unname(structure_to_iupac(glycans))
  )
})

test_that("structure_from_tibbles preserves named missing positions", {
  glycans <- c(
    first = o_glycan_core_1(),
    missing = NA,
    second = o_glycan_core_1()
  )
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)
  anomers <- get_anomer(glycans)

  rebuilt <- structure_from_tibbles(nodes, edges, anomers)

  expect_equal(names(rebuilt), names(glycans))
  expect_equal(is.na(rebuilt), is.na(glycans))
})

test_that("structure_from_tibbles handles single-node and reordered rows", {
  glycan <- as_glycan_structure("Glc3S(a1-")
  nodes <- structure_nodes(glycan)
  edges <- structure_edges(glycan)

  rebuilt <- structure_from_tibbles(
    nodes[rev(seq_len(nrow(nodes))), ],
    edges[rev(seq_len(nrow(edges))), ],
    get_anomer(glycan)
  )

  expect_equal(structure_to_iupac(rebuilt), structure_to_iupac(glycan))
})

test_that("structure_from_tibbles is insensitive to node and edge row order", {
  glycans <- c(
    first = n_glycan_core(),
    second = o_glycan_core_1(),
    third = n_glycan_core()
  )
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)

  nodes <- nodes[order(nodes$node_id, nodes$glycan_id, decreasing = TRUE), ]
  edges <- edges[order(edges$edge_id, edges$glycan_id, decreasing = TRUE), ]

  rebuilt <- structure_from_tibbles(nodes, edges, get_anomer(glycans))

  expect_equal(names(rebuilt), names(glycans))
  expect_equal(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
})

test_that("structure table round trip preserves missing positions", {
  glycans <- c(o_glycan_core_1(), glycan_structure(NA), o_glycan_core_1())
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)
  anomers <- get_anomer(glycans)

  rebuilt <- structure_from_tibbles(nodes, edges, anomers)

  expect_equal(is.na(rebuilt), is.na(glycans))
  expect_equal(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
})

test_that("structure table helpers reject invalid inputs", {
  nodes <- structure_nodes(o_glycan_core_1())
  edges <- structure_edges(o_glycan_core_1())

  expect_error(
    structure_from_tibbles(nodes[-1], edges, "a1"),
    "must contain"
  )
  expect_error(
    structure_from_tibbles(nodes, edges, character()),
    "outside"
  )
  expect_error(
    structure_from_tibbles(nodes[0, ], edges, "a1"),
    "without nodes"
  )
})
