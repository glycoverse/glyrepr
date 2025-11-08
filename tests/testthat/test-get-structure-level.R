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