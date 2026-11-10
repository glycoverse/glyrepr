test_that("installed native code loads and parses under Chinese collation", {
  # The DLL must not have been loaded before selecting the locale. A separate
  # installed-package process also puts a deadline on native initialization.
  package_path <- find.package("glyrepr")
  skip_if_not(file.exists(file.path(package_path, "Meta", "package.rds")))
  script <- tempfile(fileext = ".R")
  output <- tempfile(fileext = ".log")
  withr::defer(unlink(c(script, output)))
  libraries <- unique(c(dirname(package_path), .libPaths()))
  writeLines(
    c(
      paste0(".libPaths(", paste(deparse(libraries), collapse = ""), ")"),
      'locales <- if (.Platform$OS.type == "windows") {',
      '  c("Chinese (Simplified)_China.utf8", "Chinese_China.utf8", "Chinese_China.936")',
      '} else { c("zh_CN.UTF-8", "zh_CN.utf8") }',
      'selected <- ""',
      'for (locale in locales) {',
      '  selected <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale))',
      '  if (nzchar(selected)) break',
      '}',
      'if (!nzchar(selected)) quit(status = 42L)',
      'before <- Sys.getlocale()',
      'library(glyrepr)',
      'stopifnot(identical(Sys.getlocale(), before))',
      'x <- c("Gal6S", "Gal3/6S", "{6S|1,2}Gal(b1-4)Glc",',
      '       "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-")',
      'arrays <- glyrepr:::.compact_iupac_arrays(x)',
      'stopifnot(all(vapply(arrays, function(a) identical(a$status, "ok"), logical(1))))',
      'actual <- as_glycan_structure(x)',
      'reference <- as_glycan_structure(lapply(x, glyrepr:::.parse_iupac_condensed_single))',
      'stopifnot(identical(as.character(actual), as.character(reference)))',
      'stopifnot(identical(Sys.getlocale(), before))'
    ),
    script
  )
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = output,
    stderr = output,
    timeout = 30
  ))
  if (identical(status, 42L)) {
    skip("Chinese collation is not available on this system")
  }
  expect_identical(status, 0L, info = paste(readLines(output), collapse = "\n"))
})

test_that("native ASCII token scanning preserves parser validation", {
  valid <- c(
    "Gal6S",
    "Gal3/6S",
    "Gal?S",
    "Gal3S6P",
    "Gal6PPEtn",
    "Gal6PEtn",
    "Gal6PC",
    "Gal6NAc",
    "Gal6NGc",
    "Gal6Pyr",
    "Gal6Me",
    "Gal6Ac",
    "Gal6Gc",
    "Gal6N",
    "{6S|1,2}Gal(b1-4)Glc"
  )
  arrays <- .compact_iupac_arrays(valid)
  expect_true(all(vapply(
    arrays,
    function(a) identical(a$status, "ok"),
    logical(1)
  )))
  reference <- as_glycan_structure(lapply(valid, .parse_iupac_condensed_single))
  expect_identical(
    as.character(as_glycan_structure(valid)),
    as.character(reference)
  )
  invalid <- paste0(
    "{6S|",
    c("0", "01", "1,", ",1", "1,,2", "1,01", "-1", "1,1", "2147483648"),
    "}Gal(b1-4)Glc"
  )
  invalid <- c(invalid, "{0S}Gal", "{1//2S}Gal", "{1SX}Gal")
  invalid_arrays <- .compact_iupac_arrays(invalid)
  expect_true(all(vapply(
    invalid_arrays,
    function(a) identical(a$status, "error"),
    logical(1)
  )))
})
