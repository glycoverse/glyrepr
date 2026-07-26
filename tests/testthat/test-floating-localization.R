test_that("localize_floating_parts attaches selected parts", {
  glycans <- as_glycan_structure(c(
    floating = paste0(
      "{Fuc(a1-2)|1,2}",
      "{Neu5Ac(a2-6)|1,2}",
      "Gal(b1-3)GalNAc(a1-"
    ),
    missing = NA,
    ordinary = "Gal(b1-4)GlcNAc(b1-"
  ))
  assignments <- tibble::tibble(
    glycan_id = 1L,
    part_id = 2L,
    parent_node = 1L
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
      parent_node = 1L
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
      "{Fuc(a1-2)|1,2}",
      "{Neu5Ac(a2-6)|1,2}",
      "Gal(b1-3)GalNAc(a1-"
    )
  )

  localized <- localize_floating_parts(
    glycan,
    tibble::tibble(
      glycan_id = 1L,
      part_id = 2L,
      parent_node = 1L
    )
  )

  expect_identical(
    as.character(localized),
    "{Fuc(a1-2)|2,3}Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"
  )
  expect_identical(
    structure_floating_parts(localized)$parents,
    list(c(2L, 3L))
  )
})

test_that("localize_floating_parts returns x for no assignments", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
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
    "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
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
    restricted = "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-",
    missing = NA,
    ordinary = "Gal(a1-"
  ))

  expect_snapshot(
    localize_floating_parts(
      glycans,
      tibble::tibble(
        glycan_id = 1L,
        part_id = 1L,
        parent_node = 3L
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
      "{Fuc(a1-3)|1,2}",
      "{Neu5Ac(a2-3)|1,2}",
      "Gal(b1-4)GalNAc(a1-"
    )
  )

  expect_snapshot(
    localize_floating_parts(
      glycan,
      tibble::tibble(
        glycan_id = c(1L, 1L),
        part_id = c(1L, 2L),
        parent_node = c(1L, 1L)
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
        parent_node = 2L
      )
    ),
    error = TRUE
  )
})
