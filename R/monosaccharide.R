# This table referred to https://www.ncbi.nlm.nih.gov/glycans/snfg.html
monosaccharides <- tibble::tribble(
  ~generic  , ~concrete  , ~anomer_pos ,
  # Hexose
  "Hex"     , "Glc"      , 1L          ,
  "Hex"     , "Man"      , 1L          ,
  "Hex"     , "Gal"      , 1L          ,
  "Hex"     , "Gul"      , 1L          ,
  "Hex"     , "Alt"      , 1L          ,
  "Hex"     , "All"      , 1L          ,
  "Hex"     , "Tal"      , 1L          ,
  "Hex"     , "Ido"      , 1L          ,
  # HexNAc
  "HexNAc"  , "GlcNAc"   , 1L          ,
  "HexNAc"  , "GalNAc"   , 1L          ,
  "HexNAc"  , "ManNAc"   , 1L          ,
  "HexNAc"  , "GulNAc"   , 1L          ,
  "HexNAc"  , "AltNAc"   , 1L          ,
  "HexNAc"  , "AllNAc"   , 1L          ,
  "HexNAc"  , "TalNAc"   , 1L          ,
  "HexNAc"  , "IdoNAc"   , 1L          ,
  # Hexosamine
  "HexN"    , "GlcN"     , 1L          ,
  "HexN"    , "ManN"     , 1L          ,
  "HexN"    , "GalN"     , 1L          ,
  "HexN"    , "GulN"     , 1L          ,
  "HexN"    , "AltN"     , 1L          ,
  "HexN"    , "AllN"     , 1L          ,
  "HexN"    , "TalN"     , 1L          ,
  "HexN"    , "IdoN"     , 1L          ,
  # Hexuronate
  "HexA"    , "GlcA"     , 1L          ,
  "HexA"    , "ManA"     , 1L          ,
  "HexA"    , "GalA"     , 1L          ,
  "HexA"    , "GulA"     , 1L          ,
  "HexA"    , "AltA"     , 1L          ,
  "HexA"    , "AllA"     , 1L          ,
  "HexA"    , "TalA"     , 1L          ,
  "HexA"    , "IdoA"     , 1L          ,
  # Deoxyhexose
  "dHex"    , "Fuc"      , 1L          ,
  "dHex"    , "Qui"      , 1L          ,
  "dHex"    , "Rha"      , 1L          ,
  "dHex"    , "6dGul"    , 1L          ,
  "dHex"    , "6dAlt"    , 1L          ,
  "dHex"    , "6dTal"    , 1L          ,
  # DeoxyhexNAc
  "dHexNAc" , "QuiNAc"   , 1L          ,
  "dHexNAc" , "RhaNAc"   , 1L          ,
  "dHexNAc" , "6dAltNAc" , 1L          ,
  "dHexNAc" , "6dTalNAc" , 1L          ,
  "dHexNAc" , "FucNAc"   , 1L          ,
  # Di-deoxyhexose
  "ddHex"   , "Oli"      , 1L          ,
  "ddHex"   , "Tyv"      , 1L          ,
  "ddHex"   , "Abe"      , 1L          ,
  "ddHex"   , "Par"      , 1L          ,
  "ddHex"   , "Dig"      , 1L          ,
  "ddHex"   , "Col"      , 1L          ,
  # Pentose
  "Pen"     , "Ara"      , 1L          ,
  "Pen"     , "Lyx"      , 1L          ,
  "Pen"     , "Xyl"      , 1L          ,
  "Pen"     , "Rib"      , 1L          ,
  # 3-deoxy-nonulosonic acids
  "NeuAc"   , "Neu5Ac"   , 2L          ,
  "NeuGc"   , "Neu5Gc"   , 2L          ,
  "gNeu"    , "Neu"      , 2L          ,
  "gKdn"    , "Kdn"      , 2L          ,
  # 3,9-dideoxy-nonulosonic acids
  "gPse"    , "Pse"      , 2L          ,
  "gLeg"    , "Leg"      , 2L          ,
  "gAci"    , "Aci"      , 2L          ,
  "g4eLeg"  , "4eLeg"    , 2L          ,
  # Unknown
  "gBac"    , "Bac"      , 1L          ,
  "Hep"     , "LDmanHep" , 1L          ,
  "gKdo"    , "Kdo"      , 2L          ,
  "HepA"    , "Dha"      , 2L          ,
  "Hep"     , "DDmanHep" , 1L          ,
  "MurAc"   , "MurNAc"   , 1L          ,
  "MurGc"   , "MurNGc"   , 1L          ,
  "gMur"    , "Mur"      , 1L          ,
  # Assigned
  "Pen"     , "Api"      , 1L          ,
  "Hex"     , "Fru"      , 2L          ,
  "Hex"     , "Tag"      , 2L          ,
  "Hex"     , "Sor"      , 2L          ,
  "Hex"     , "Psi"      , 2L
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
    monos <- monosaccharides[c("generic", "concrete")]
    unique(purrr::discard(unlist(monos, use.names = FALSE), is.na))
  } else {
    unique(purrr::discard(monosaccharides[[mono_type]], is.na))
  }
}


#' Get Anomer Positions
#'
#' This function returns the anomer position for concrete monosaccharide names.
#'
#' @param mono A character vector of concrete monosaccharide names.
#'
#' @returns An integer vector of anomer positions.
#'
#' @examples
#' get_anomer_pos(c("Gal", "Neu5Ac"))
#'
#' @export
get_anomer_pos <- function(mono) {
  checkmate::assert_character(mono, any.missing = FALSE)

  anomer_pos <- monosaccharides$anomer_pos[match(
    mono,
    monosaccharides$concrete
  )]
  unknown_monos <- mono[is.na(anomer_pos)]

  if (length(unknown_monos) > 0) {
    cli::cli_abort(c(
      "{.arg mono} must contain only concrete monosaccharide names.",
      "x" = "Invalid value{?s}: {.val {unique(unknown_monos)}}.",
      "i" = "Call {.fun available_monosaccharides} with {.val concrete} to see supported names."
    ))
  }

  names(anomer_pos) <- names(mono)
  anomer_pos
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
  (mono %in% monosaccharides$generic | mono %in% monosaccharides$concrete)
}
