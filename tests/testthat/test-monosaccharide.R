patrick::with_parameters_test_that(
  "known monosaccharides",
  {
    expect_true(is_known_monosaccharide(mono))
  },
  mono = c("Gal", "Hex")
)


patrick::with_parameters_test_that(
  "unknown monosaccharides",
  {
    expect_false(is_known_monosaccharide(mono))
  },
  mono = c("X", "Hx", "Nac")
)


test_that("get all generic monosaccharides", {
  res <- available_monosaccharides("generic")
  expect_contains(res, c("Hex", "HexNAc", "dHex"))
  expect_true(length(unique(res)) == length(res))
  expect_identical(res, unique(monosaccharide_definitions$generic))
})


test_that("get all concrete monosaccharides", {
  res <- available_monosaccharides("concrete")
  expect_contains(res, c("Gal", "Man", "GlcNAc", "Galf", "GlcfNAc"))
  expect_true(length(unique(res)) == length(res))
})


test_that("every concrete monosaccharide has a furanose form", {
  ringless <- names(furanose_monosaccharides)
  furanose <- unname(furanose_monosaccharides)

  expect_length(furanose, length(ringless))
  expect_setequal(ringless, monosaccharide_definitions$concrete)
  expect_setequal(
    available_monosaccharides("concrete"),
    c(ringless, furanose)
  )
  expect_identical(
    is_known_monosaccharide(furanose),
    rep(TRUE, length(furanose))
  )
  expect_identical(
    get_mono_type(furanose),
    rep("concrete", length(furanose))
  )
  expect_identical(
    infer_anomer_pos(furanose),
    infer_anomer_pos(ringless)
  )
  expect_identical(
    convert_to_generic(furanose),
    convert_to_generic(ringless)
  )
})


test_that("unusual configurations are appended to the concrete vocabulary", {
  natural <- natural_monosaccharide_definitions$concrete
  natural_furanose <- unname(natural_furanose_monosaccharides)
  unusual <- unname(
    unusual_configuration_monosaccharides[configuration_monos]
  )
  unusual_furanose <- unname(
    unusual_configuration_furanose_monosaccharides[unusual]
  )

  expect_identical(
    available_monosaccharides("concrete"),
    c(natural, natural_furanose, unusual, unusual_furanose)
  )
  expect_contains(unusual, c("DFuc", "LGul", "LNeu5Ac", "LKdn"))
  expect_contains(unusual_furanose, c("DFucf", "LGulf", "LNeuf5Ac"))
})


test_that("every unusual configuration has concrete residue behavior", {
  natural <- names(unusual_configuration_monosaccharides)
  unusual <- unname(unusual_configuration_monosaccharides)

  expect_identical(
    is_known_monosaccharide(unusual),
    rep(TRUE, length(unusual))
  )
  expect_identical(
    get_mono_type(unusual),
    rep("concrete", length(unusual))
  )
  expect_identical(
    unname(infer_anomer_pos(unusual)),
    unname(infer_anomer_pos(natural))
  )
  expect_identical(
    unname(convert_to_generic(unusual)),
    unname(convert_to_generic(natural))
  )
})


test_that("natural configurations remain unprefixed", {
  explicit_natural <- paste0(
    natural_monosaccharide_configurations,
    names(natural_monosaccharide_configurations)
  )

  expect_identical(
    is_known_monosaccharide(explicit_natural),
    rep(FALSE, length(explicit_natural))
  )
  expect_identical(
    is_known_monosaccharide(c("LFuc", "DGul", "DNeu5Ac")),
    rep(FALSE, 3)
  )
})


test_that("configuration-unspecified residues have no unusual forms", {
  prefixed <- c(
    paste0("D", configuration_unspecified_monosaccharides),
    paste0("L", configuration_unspecified_monosaccharides)
  )

  expect_identical(
    is_known_monosaccharide(prefixed),
    rep(FALSE, length(prefixed))
  )
})


test_that("get all monosaccharides", {
  res <- available_monosaccharides()
  expect_contains(res, c("Hex", "Man"))
  expect_true(length(unique(res)) == length(res))
})


test_that("infer_anomer_pos returns anomer positions for concrete monosaccharides", {
  monos <- c("Gal", "GalNAc", "Neu5Ac", "Kdo", "Fru")

  expect_equal(infer_anomer_pos(monos), c(1L, 1L, 2L, 2L, 2L))
})


test_that("infer_anomer_pos returns anomer positions for generic monosaccharides", {
  monos <- c("Hex", "HexNAc", "NeuAc", "gKdo")

  expect_equal(infer_anomer_pos(monos), c(1L, 1L, 2L, 2L))
})


test_that("infer_anomer_pos supports mixed generic and concrete monosaccharides", {
  monos <- c("Gal", "Hex", "Neu5Ac", "NeuAc")

  expect_equal(infer_anomer_pos(monos), c(1L, 1L, 2L, 2L))
})


test_that("infer_anomer_pos preserves names", {
  monos <- c(a = "Gal", b = "Neu5Ac")

  expect_equal(names(infer_anomer_pos(monos)), names(monos))
})


test_that("infer_anomer_pos rejects unknown monosaccharides", {
  expect_error(
    infer_anomer_pos("X"),
    regexp = "known monosaccharide"
  )
})


test_that("get_anomer_pos remains an alias for infer_anomer_pos", {
  monos <- c(a = "Hex", b = "Neu5Ac")

  expect_equal(get_anomer_pos(monos), infer_anomer_pos(monos))
})
