#' @export
print.glycan_graph <- function(x, ...) {
  if (inherits(x, "ne_glycan_graph")) {
    cli::cat_line("Glycan Graph (NE)")
  } else {
    cli::cat_line("Glycan Graph (DN)")
  }
  composition <- get_composition(x)
  comp_str <- paste(names(composition), composition, sep = ": ", collapse = ", ")
  cli::cat_line(comp_str)
}
