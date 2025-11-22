test_that("get_structure_level works", {
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-?)GalNAc(a1-",
    "Gal(??-?)GalNAc(??-",
    "Hex(??-?)HexNAc(??-",
    "Hex(b1-3)HexNAc(a1-"
  ))
  expect_equal(get_structure_level(glycans), c("intact", "partial", "topological", "basic", "basic"))
})

test_that("reduce_structure_level works", {
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",  # intact
    "Gal(b1-?)GalNAc(a1-",  # partial
    "Gal(??-?)GalNAc(??-"  # topological
  ))
  expect_equal(
    as.character(reduce_structure_level(glycans, to_level = "topological")),
    c(
      "Gal(??-?)GalNAc(??-",
      "Gal(??-?)GalNAc(??-",
      "Gal(??-?)GalNAc(??-"
    )
  )
  expect_equal(
    as.character(reduce_structure_level(glycans, to_level = "basic")),
    c(
      "Hex(??-?)HexNAc(??-",
      "Hex(??-?)HexNAc(??-",
      "Hex(??-?)HexNAc(??-"
    )
  )
})

test_that("reduce_structure_level rejects higher level", {
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",  # intact
    "Gal(b1-?)GalNAc(a1-",  # partial
    "Gal(??-?)GalNAc(??-",  # topological
    "Hex(??-?)HexNAc(??-"  # basic
  ))
  expect_snapshot(reduce_structure_level(glycans, to_level = "topological"), error = TRUE)
})