# This table referred to https://www.ncbi.nlm.nih.gov/glycans/snfg.html
monosaccharides <- tibble::tribble(
  ~generic, ~concrete,
  # Hexose
  "Hex", "Glc",
  "Hex", "Man",
  "Hex", "Gal",
  "Hex", "Gul",
  "Hex", "Alt",
  "Hex", "All",
  "Hex", "Tal",
  "Hex", "Ido",
  # HexNAc
  "HexNAc", "GlcNAc",
  "HexNAc", "GalNAc",
  "HexNAc", "ManNAc",
  "HexNAc", "GulNAc",
  "HexNAc", "AltNAc",
  "HexNAc", "AllNAc",
  "HexNAc", "TalNAc",
  "HexNAc", "IdoNAc",
  # Hexosamine
  "HexN", "GlcN",
  "HexN", "ManN",
  "HexN", "GalN",
  "HexN", "GulN",
  "HexN", "AltN",
  "HexN", "AllN",
  "HexN", "TalN",
  "HexN", "IdoN",
  # Hexuronate
  "HexA", "GlcA",
  "HexA", "ManA",
  "HexA", "GalA",
  "HexA", "GulA",
  "HexA", "AltA",
  "HexA", "AllA",
  "HexA", "TalA",
  "HexA", "IdoA",
  # Deoxyhexose
  "dHex", "Fuc",
  "dHex", "Qui",
  "dHex", "Rha",
  "dHex", "6dGul",
  "dHex", "6dAlt",
  "dHex", "6dTal",
  # DeoxyhexNAc
  "dHexNAc", "QuiNAc",
  "dHexNAc", "RhaNAc",
  "dHexNAc", "6dAltNAc",
  "dHexNAc", "6dTalNAc",
  "dHexNAc", "FucNAc",
  # Di-deoxyhexose
  "ddHex", "Oli",
  "ddHex", "Tyv",
  "ddHex", "Abe",
  "ddHex", "Par",
  "ddHex", "Dig",
  "ddHex", "Col",
  # Pentose
  "Pent", "Ara",
  "Pent", "Lyx",
  "Pent", "Xyl",
  "Pent", "Rib",
  # 3-deoxy-nonulosonic acids
  "NeuAc", "Neu5Ac",
  "NeuGc", "Neu5Gc",
  NA, "Sia",
  NA, "Neu",
  NA, "Kdn",
  # 3,9-dideoxy-nonulosonic acids
  NA, "Pse",
  NA, "Leg",
  NA, "Aci",
  NA, "4eLeg",
  # Unknown
  NA, "Bac",
  NA, "LDmanHep",
  NA, "Kdo",
  NA, "Dha",
  NA, "DDmanHep",
  NA, "MurNAc",
  NA, "MurNGc",
  NA, "Mur",
  # Assigned
  NA, "Api",
  NA, "Fru",
  NA, "Tag",
  NA, "Sor",
  NA, "Psi"
)


#' Get Available Monosaacharides
#'
#' This function returns a character vector of monosaccharide names of
#' the given type. See [get_mono_type()] for monosaacharide types.
#'
#' @param mono_type A character string specifying the type of monosaccharides.
#'  Can be "all", "generic", or "concrete". Default is "all".
#'
#' @returns A character vector of monosaccharide names.
#'
#' @examples
#' available_monosaccharides()
#' 
#' @export
available_monosaccharides <- function(mono_type = "all") {
  checkmate::assert_choice(mono_type, c("all", "generic", "concrete"))
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
#' @returns A logical vector.
#'
#' @examples
#' is_known_monosaccharide(c("Gal", "Hex"))
#' is_known_monosaccharide(c("X", "Hx", "Nac"))
#'
#' @export
is_known_monosaccharide <- function(mono) {
  checkmate::assert_character(mono)
  (
    mono %in% monosaccharides$generic |
    mono %in% monosaccharides$concrete
  )
}
