#' Get Available Substituents
#'
#' Get the available substituents for monosaccharides.
#'
#' @return A character vector.
#' @export
available_substituents <- function() {
  c("Me", "Ac", "NAc", "P", "S", "Pyr", "PC", "PPEtn", "PEtn", "N")
}
