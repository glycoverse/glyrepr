test_that("monosaccharides with substituents are colored correctly", {
  expect_equal(get_mono_color("Neu5Ac9Ac"), "#A54399") # Purple for Neu5Ac
  expect_equal(get_mono_color("Neu5Gc9Ac"), "#8FCCE9") # Light blue for Neu5Gc
  expect_equal(get_mono_color("Gal6S"), "#FFD400") # Yellow for Gal
  expect_equal(get_mono_color("Gal4/6S"), "#FFD400") # Yellow for Gal
  expect_equal(get_mono_color("Man3Me"), "#00A651") # Green for Man
  expect_equal(get_mono_color("Glc3Pyr"), "#0072BC") # Blue for Glc
  expect_equal(get_mono_color("Glc"), "#0072BC") # Blue for Glc
})

test_that("furanose forms use the same colors as their ringless forms", {
  ringless <- names(furanose_monosaccharides)
  furanose <- unname(furanose_monosaccharides)

  expect_identical(
    purrr::map_chr(furanose, get_mono_color),
    purrr::map_chr(ringless, get_mono_color)
  )
  expect_identical(get_mono_color("Galf3Me"), get_mono_color("Gal3Me"))
  expect_identical(
    get_mono_color("Neuf5Ac9Ac"),
    get_mono_color("Neu5Ac9Ac")
  )
})

test_that("unusual configurations use their natural counterparts' colors", {
  natural <- names(unusual_configuration_monosaccharides)
  unusual <- unname(unusual_configuration_monosaccharides)

  expect_identical(
    purrr::map_chr(unusual, get_mono_color),
    purrr::map_chr(natural, get_mono_color)
  )
  expect_identical(get_mono_color("D-Fuc3S"), get_mono_color("Fuc3S"))
  expect_identical(
    get_mono_color("L-Neu5Ac9Ac"),
    get_mono_color("Neu5Ac9Ac")
  )
})

test_that("special sialic acid cases are handled correctly", {
  # Test color assignment for special cases
  expect_equal(get_mono_color("Neu4Ac5Ac"), "#A54399") # Purple for Neu5Ac
  expect_equal(get_mono_color("Neu4Ac5Gc"), "#8FCCE9") # Light blue for Neu5Gc
})

test_that("all valid linkage forms can be colored", {
  old_options <- options(cli.num_colors = 256)
  on.exit(options(old_options), add = TRUE)
  annotations <- c("a1-3", "??-?", "a2-3/6", "a?-?", "?1-4")

  purrr::walk(annotations, function(annotation) {
    text <- paste0("(", annotation, ")")
    colored <- add_gray_linkages(text)
    expect_true(cli::ansi_has_any(colored))
    expect_identical(cli::ansi_strip(colored), text)
  })

  incomplete <- add_gray_linkages("(??-")
  expect_true(cli::ansi_has_any(incomplete))
  expect_identical(cli::ansi_strip(incomplete), "(??-")
})
