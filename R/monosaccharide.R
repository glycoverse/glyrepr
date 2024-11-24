# This table referred to https://www.ncbi.nlm.nih.gov/glycans/snfg.html
monosaccharides <- tibble::tribble(
  ~simple, ~generic, ~concrete,
  # Hexose
  "H", "Hex", "Glc",
  "H", "Hex", "Man",
  "H", "Hex", "Gal",
  NA, "Hex", "Gul",
  NA, "Hex", "Alt",
  NA, "Hex", "All",
  NA, "Hex", "Tal",
  NA, "Hex", "Ido",
  # HexNAc
  "N", "HexNAc", "GlcNAc",
  "N", "HexNAc", "GalNAc",
  NA, "HexNAc", "ManNAc",
  NA, "HexNAc", "GulNAc",
  NA, "HexNAc", "AltNAc",
  NA, "HexNAc", "AllNAc",
  NA, "HexNAc", "TalNAc",
  NA, "HexNAc", "IdoNAc",
  # Hexosamine
  NA, "HexN", "GlcN",
  NA, "HexN", "ManN",
  NA, "HexN", "GalN",
  NA, "HexN", "GulN",
  NA, "HexN", "AltN",
  NA, "HexN", "AllN",
  NA, "HexN", "TalN",
  NA, "HexN", "IdoN",
  # Hexuronate
  NA, "HexA", "GlcA",
  NA, "HexA", "ManA",
  NA, "HexA", "GalA",
  NA, "HexA", "GulA",
  NA, "HexA", "AltA",
  NA, "HexA", "AllA",
  NA, "HexA", "TalA",
  NA, "HexA", "IdoA",
  # Deoxyhexose
  "F", "dHex", "Fuc",
  NA, "dHex", "Qui",
  NA, "dHex", "Rha",
  NA, "dHex", "6dGul",
  NA, "dHex", "6dAlt",
  NA, "dHex", "6dTal",
  # DeoxyhexNAc
  NA, "dHexNAc", "QuiNAc",
  NA, "dHexNAc", "RhaNAc",
  NA, "dHexNAc", "6dAltNAc",
  NA, "dHexNAc", "6dTalNAc",
  NA, "dHexNAc", "FucNAc",
  # Di-deoxyhexose
  NA, "ddHex", "Oli",
  NA, "ddHex", "Tyv",
  NA, "ddHex", "Abe",
  NA, "ddHex", "Par",
  NA, "ddHex", "Dig",
  NA, "ddHex", "Col",
  # Pentose
  NA, "Pent", "Ara",
  NA, "Pent", "Lyx",
  NA, "Pent", "Xyl",
  NA, "Pent", "Rib",
  # 3-deoxy-nonulosonic acids
  "A", "NeuAc", "Neu5Ac",
  "G", "NeuGc", "Neu5Gc",
  "S", NA, "Sia",
  NA, NA, "Neu",
  NA, NA, "Kdn",
  # 3,9-dideoxy-nonulosonic acids
  NA, NA, "Pse",
  NA, NA, "Leg",
  NA, NA, "Aci",
  NA, NA, "4eLeg",
  # Unknown
  NA, NA, "Bac",
  NA, NA, "LDmanHep",
  NA, NA, "Kdo",
  NA, NA, "Dha",
  NA, NA, "DDmanHep",
  NA, NA, "MurNAc",
  NA, NA, "MurNGc",
  NA, NA, "Mur",
  # Assigned
  NA, NA, "Api",
  NA, NA, "Fru",
  NA, NA, "Tag",
  NA, NA, "Sor",
  NA, NA, "Psi"
)


#' Get Available Monosaacharides
#'
#' This function returns a character vector of monosaccharide names of
#' the given type. See [decide_mono_type()] for monosaacharide types.
#'
#' @param mono_type A character string specifying the type of monosaccharides.
#'  Can be "all", "simple", "generic", or "concrete". Default is "all".
#'
#' @return A character vector of monosaccharide names.
#' @export
available_monosaccharides <- function(mono_type = "all") {
  checkmate::assert_choice(mono_type, c("all", "simple", "generic", "concrete"))
  if (mono_type == "all") {
    unique(purrr::discard(unlist(monosaccharides, use.names = FALSE), is.na))
  } else {
    unique(purrr::discard(monosaccharides[[mono_type]], is.na))
  }
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
  checkmate::assert_character(mono)
  (
    mono %in% monosaccharides$simple |
    mono %in% monosaccharides$generic |
    mono %in% monosaccharides$concrete
  )
}
