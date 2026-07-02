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

test_that("structure_from_tibbles recreates structure vectors", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  nodes <- structure_nodes(glycans)
  edges <- structure_edges(glycans)
  anomers <- get_anomer(glycans)

  rebuilt <- structure_from_tibbles(nodes, edges, anomers)

  expect_s3_class(rebuilt, "glyrepr_structure")
  expect_equal(structure_to_iupac(rebuilt), structure_to_iupac(glycans))
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
