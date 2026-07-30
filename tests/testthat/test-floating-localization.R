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
      parent_node = 3L
    )
  )

  expect_identical(
    as.character(localized),
    "{Fuc(a1-2)|2,3}Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"
  )
  expect_identical(
    structure_floating_parts(localized)$parents,
    list(c(3L, 4L))
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
      "{Fuc(a1-3)|1,2}",
      "{Neu5Ac(a2-3)|1,2}",
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

test_that("enumerate_floating_localizations handles ambiguous linkages", {
  glycan <- as_glycan_structure(
    "{Neu5Ac(a2-3/6)|1,2}Gal(b1-4)GalNAc(a1-"
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
      "{Neu5Ac(a2-3)|1,2}",
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
      "{Neu5Ac(a2-3)|1,2}",
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

test_that("enumerate_floating_localizations retains every input position", {
  glycans <- as_glycan_structure(c(
    missing = NA,
    ordinary = "Gal(a1-",
    floating = "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
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
      "{Fuc(a1-3)|1,2}",
      "{Neu5Ac(a2-3)|1,2}",
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
