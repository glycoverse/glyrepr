test_that("get_structure_level works for each glycan separately", {
  glycan1 <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  glycan2 <- as_glycan_structure("Gal(b1-?)GalNAc(a1-")
  glycan3 <- as_glycan_structure("Gal(??-?)GalNAc(??-")
  glycan4 <- as_glycan_structure("Hex(??-?)HexNAc(??-")
  glycan5 <- as_glycan_structure("Hex(b1-3)HexNAc(a1-")

  expect_equal(get_structure_level(glycan1), "intact")
  expect_equal(get_structure_level(glycan2), "partial")
  expect_equal(get_structure_level(glycan3), "topological")
  expect_equal(get_structure_level(glycan4), "basic")
  expect_equal(get_structure_level(glycan5), "basic")
})

test_that("get_structure_level works for multiple glycans", {
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "GalNAc(a1-"
  ))
  expect_equal(get_structure_level(glycans), c("intact", "intact"))
})

test_that("reduce_structure_level works for each glycan separately", {
  glycan_intact <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  glycan_partial <- as_glycan_structure("Gal(b1-?)GalNAc(a1-")
  glycan_topological <- as_glycan_structure("Gal(??-?)GalNAc(??-")

  # intact to topological
  expect_equal(
    as.character(reduce_structure_level(glycan_intact, to_level = "topological")),
    "Gal(??-?)GalNAc(??-"
  )
  # partial to topological
  expect_equal(
    as.character(reduce_structure_level(glycan_partial, to_level = "topological")),
    "Gal(??-?)GalNAc(??-"
  )
  # topological to topological (identity)
  expect_equal(
    as.character(reduce_structure_level(glycan_topological, to_level = "topological")),
    "Gal(??-?)GalNAc(??-"
  )

  # intact to basic
  expect_equal(
    as.character(reduce_structure_level(glycan_intact, to_level = "basic")),
    "Hex(??-?)HexNAc(??-"
  )
  # partial to basic
  expect_equal(
    as.character(reduce_structure_level(glycan_partial, to_level = "basic")),
    "Hex(??-?)HexNAc(??-"
  )
  # topological to basic
  expect_equal(
    as.character(reduce_structure_level(glycan_topological, to_level = "basic")),
    "Hex(??-?)HexNAc(??-"
  )
})

test_that("reduce_structure_level rejects higher level", {
  glycan <- as_glycan_structure("Hex(??-?)HexNAc(??-")
  expect_snapshot(reduce_structure_level(glycan, to_level = "topological"), error = TRUE)
})

test_that("reduce_structure_level works for multiple glycans", {
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-?)GalNAc(a1-",
    "Gal(??-?)GalNAc(??-"
  ))
  expect_equal(
    as.character(reduce_structure_level(glycans, to_level = "topological")),
    c("Gal(??-?)GalNAc(??-", "Gal(??-?)GalNAc(??-", "Gal(??-?)GalNAc(??-")
  )
})