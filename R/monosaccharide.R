#' @importFrom tidyr unnest
monosaccharides <- unnest(cols = concrete, tibble::enframe(name = "generic", value = "concrete", list(
  "Hex" = c("Glc", "Man", "Gal", "Gul", "Alt", "All", "Tal", "Ido"),
  "HexNAc" = c("GlcNAc", "ManNAc", "GalNAc", "GulNAc",  "AltNAc", "AllNAc", "TalNAc", "IdoNAc"),
  "HexN" = c("GlcN", "ManN", "GalN", "GulN", "AltN", "AllN", "TalN", "IdoN"),
  "HexA" = c("GlcA", "ManA", "GalA", "GulA", "AltA", "AllA", "TalA", "IdoA"),
  "dHex" = c("Qui", "Rha", "6dGul", "6dAlt", "6dTal", "Fuc"),
  "dHexNAc" = c("QuiNAc", "RhaNAc", "6dAltNAc", "6dTalNAc", "FucNAc"),
  "ddHex" = c("Oli", "Tyv", "Abe", "Par", "Dig", "Col"),
  "Pent" = c("Ara", "Lyx", "Xyl", "Rib"),
  "NeuAc" = "Neu5Ac",
  "NeuGc" = "Neu5Gc"
)))
