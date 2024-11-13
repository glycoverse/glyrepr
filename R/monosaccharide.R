monosaccharides <- tibble::tribble(
  ~simple, ~generic, ~concrete,
  # Hexoses
  "H", "Hex", "Glc",
  "H", "Hex", "Man",
  "H", "Hex", "Gal",
  NA, "Hex", "Gul",
  NA, "Hex", "Alt",
  NA, "Hex", "All",
  NA, "Hex", "Tal",
  NA, "Hex", "Ido",
  # HexNAcs
  "N", "HexNAc", "GlcNAc",
  "N", "HexNAc", "GalNAc",
  NA, "HexNAc", "ManNAc",
  NA, "HexNAc", "GulNAc",
  NA, "HexNAc", "AltNAc",
  NA, "HexNAc", "AllNAc",
  NA, "HexNAc", "TalNAc",
  NA, "HexNAc", "IdoNAc",
  # Hexosamines
  NA, "HexN", "GlcN",
  NA, "HexN", "ManN",
  NA, "HexN", "GalN",
  NA, "HexN", "GulN",
  NA, "HexN", "AltN",
  NA, "HexN", "AllN",
  NA, "HexN", "TalN",
  NA, "HexN", "IdoN",
  # Hexuronates
  NA, "HexA", "GlcA",
  NA, "HexA", "ManA",
  NA, "HexA", "GalA",
  NA, "HexA", "GulA",
  NA, "HexA", "AltA",
  NA, "HexA", "AllA",
  NA, "HexA", "TalA",
  NA, "HexA", "IdoA",
  # Deoxyhexoses
  "F", "dHex", "Fuc",
  NA, "dHex", "Qui",
  NA, "dHex", "Rha",
  NA, "dHex", "6dGul",
  NA, "dHex", "6dAlt",
  NA, "dHex", "6dTal",
  # DeoxyhexNAcs
  NA, "dHexNAc", "QuiNAc",
  NA, "dHexNAc", "RhaNAc",
  NA, "dHexNAc", "6dAltNAc",
  NA, "dHexNAc", "6dTalNAc",
  NA, "dHexNAc", "FucNAc",
  # Dideoxyhexoses
  NA, "ddHex", "Oli",
  NA, "ddHex", "Tyv",
  NA, "ddHex", "Abe",
  NA, "ddHex", "Par",
  NA, "ddHex", "Dig",
  NA, "ddHex", "Col",
  # Pentoses
  NA, "Pent", "Ara",
  NA, "Pent", "Lyx",
  NA, "Pent", "Xyl",
  NA, "Pent", "Rib",
  # Sialic acids
  "A", "NeuAc", "Neu5Ac",
  "G", "NeuGc", "Neu5Gc",
  "S", "Sia", "Sia",
  NA, "Kdn", "Kdn"
)


#' Get Available Monosaacharides
#'
#' This function returns a tibble of available monosaccharides.
#' The tibble has three columns: `simple`, `generic`, and `concrete`.
#' "simple" is the simple monosaccharide name, e.g. "H", "N", "F", etc.
#' "generic" is the generic monosaccharide name, e.g. "Hex", "HexNAc", "dHex", etc.
#' "concrete" is the concrete monosaccharide name, e.g. "Glc", "Gal", "GlcNAc", etc.
#'
#' @return A tibble of available monosaccharides.
#' @export
available_monosaccharides <- function() {
  monosaccharides
}


#' Check if a Monosaccharide is Known
#'
#' This function checks if a vector of monosaccharide names are known.
#'
#' @param mono A character vector of monosaccharide names.
#'
#' @return A logical vector.
#' @export
is_known_monosaccharide <- function(mono) {
  (
    mono %in% monosaccharides$simple |
    mono %in% monosaccharides$generic |
    mono %in% monosaccharides$concrete
  )
}
