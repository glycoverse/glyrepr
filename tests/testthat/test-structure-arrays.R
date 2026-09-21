test_that("array construction preserves names, missing values and duplicates", {
  a <- list(
    mono = c("Glc", "Gal"),
    sub = c("", ""),
    edges = c(1L, 2L),
    linkage = "b1-4",
    anomer = "?1"
  )
  result <- structure_from_arrays(list(a = a, missing = NULL, duplicate = a))
  expect_identical(
    as.character(result),
    c(
      a = "Gal(b1-4)Glc(?1-",
      missing = NA_character_,
      duplicate = "Gal(b1-4)Glc(?1-"
    )
  )
  expect_length(attr(result, "graphs"), 1)
  expect_identical(structure_from_arrays(list()), glycan_structure())
})

test_that("arrays preserve floating metadata and reducing end state", {
  a <- list(
    mono = c("Glc", "Gal", "Neu5Ac"),
    sub = rep("", 3),
    edges = c(1L, 2L),
    linkage = "b1-4",
    anomer = "?1",
    alditol = TRUE,
    floating_parts = list(list(
      root = 3L,
      nodes = 3L,
      linkage = "a2-6",
      parents = c(1L, 2L)
    )),
    floating_substituents = list(list(substituent = "3S", parents = c(1L, 2L)))
  )
  expected <- as_glycan_structure(
    "{3S|2,3}{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc-ol(?1-"
  )
  result <- structure_from_arrays(list(a))
  expect_identical(as.character(result), as.character(expected))
  expect_equal(
    structure_floating_parts(result),
    structure_floating_parts(expected)
  )
  expect_equal(
    structure_floating_substituents(result),
    structure_floating_substituents(expected)
  )
})

test_that("malformed array records fail before graph construction", {
  expect_snapshot(error = TRUE, structure_from_arrays(list(list(mono = "Gal"))))
  a <- list(
    mono = "Gal",
    sub = "",
    edges = c(1, 2),
    linkage = "b1-4",
    anomer = "?1"
  )
  expect_snapshot(error = TRUE, structure_from_arrays(list(a)))
  expect_snapshot(
    result <- structure_from_arrays(
      list(bad = a, absent = NULL),
      on_failure = "na"
    )
  )
  expect_identical(is.na(result), c(bad = TRUE, absent = TRUE))
})
