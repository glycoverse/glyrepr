skip_on_old_win <- function() {
  if (getRversion() < "4.4.0" && .Platform$OS.type == "windows") {
    testthat::skip("Skipping test on Windows with R version < 4.4.0")
  }
}
