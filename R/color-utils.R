#' Get Color for Concrete Monosaccharides
#'
#' @param mono A monosaccharide name (character), potentially with substituents
#' @return A color code (character)
#' @keywords internal
get_mono_color <- function(mono) {
  # Extract base monosaccharide name (without substituents)
  base_mono <- .extract_base_mono(mono)

  dplyr::case_match(
    base_mono,
    c("Glc", "GlcNAc", "GlcN", "GlcA", "Qui", "QuiNAc", "Oli", "Bac", "Api") ~ "#0072BC",
    c("Man", "ManNAc", "ManN", "ManA", "Rha", "RhaNAc", "Tyv", "Ara", "Kdn", "Pse", "LDmanHep", "Fru") ~ "#00A651",
    c("Gal", "GalNAc", "GalN", "GalA", "Lyx", "Leg", "Kdo", "Tag") ~ "#FFD400",
    c("Gul", "GulNAc", "GulN", "GulA", "6dGul", "Abe", "Xyl", "Dha", "Sor") ~ "#F47920",
    c("Alt", "AltNAc", "AltN", "AltA", "6dAlt", "6dAltNAc", "Par", "Rib", "Aci", "DDmanHep", "Psi") ~ "#F69EA1",
    c("All", "AllNAc", "AllN", "AllA", "Dig", "Neu5Ac", "MurNAc") ~ "#A54399",
    c("Tal", "TalNAc", "TalN", "TalA", "6dTal", "6dTalNAc", "Col", "Neu5Gc", "4eLeg", "MurNGc") ~ "#8FCCE9",
    c("Ido", "IdoNAc", "IdoN", "IdoA", "Neu", "Mur") ~ "#A17A4D",
    c("Fuc", "FucNAc", "Sia") ~ "#ED1C24",
    .default = "black"
  )
}

#' Add Colors to Monosaccharides
#'
#' @param monos A character vector of monosaccharide names
#' @param colored A logical value indicating whether to add colors
#' @return A character vector with ANSI color codes
#' @keywords internal
add_colors <- function(monos, colored = TRUE) {
  if (!colored) {
    return(monos)
  }
  
  colors <- purrr::map_chr(monos, get_mono_color)
  color_styles <- purrr::map(colors, cli::make_ansi_style)
  purrr::map2_chr(monos, color_styles, ~ .y(.x))
}



#' Replace Monosaccharides in String with Colored Versions
#'
#' @param text Character string containing monosaccharide names
#' @param mono_names Character vector of monosaccharide names to replace
#' @return Character string with monosaccharides replaced by colored versions
#' @keywords internal
replace_monos_with_colored <- function(text, mono_names) {
  # Get unique monosaccharides and their colors
  unique_monos <- unique(mono_names)
  colored_monos <- add_colors(unique_monos, colored = TRUE)
  names(colored_monos) <- unique_monos

  # Replace each monosaccharide in the text with its colored version
  result <- text
  for (i in seq_along(unique_monos)) {
    mono <- unique_monos[i]
    colored_mono <- colored_monos[i]

    # Special handling for monosaccharides that might have substituents
    # For example, "Neu5Ac" should match in "Neu5Ac9Ac"
    if (mono %in% c("Neu5Ac", "Neu5Gc", "Neu")) {
      # For sialic acids, match the base name followed by optional substituents
      pattern <- paste0("\\b", mono, "(?=[0-9]|\\(|$)")
      result <- stringr::str_replace_all(result, pattern, colored_mono)
    } else {
      # For other monosaccharides, use more flexible matching
      # Match the mono name followed by optional substituents (digits + letters)
      pattern <- paste0("\\b", mono, "(?=[0-9]|\\(|$)")
      result <- stringr::str_replace_all(result, pattern, colored_mono)
    }
  }

  result
}

#' Add Gray Color to Linkages in IUPAC String
#'
#' @param iupac_text Character string of IUPAC notation
#' @return Character string with linkages colored gray
#' @keywords internal
add_gray_linkages <- function(iupac_text) {
  # Gray color style
  gray_style <- cli::make_ansi_style("gray")
  
  result <- iupac_text
  
  # Pattern 1: Complete linkages like (b1-3), (a1-6)
  complete_pattern <- "\\(([ab?]\\d*-\\d*)\\)"
  result <- stringr::str_replace_all(result, complete_pattern, function(match) {
    linkage <- stringr::str_sub(match, 2, -2)  # Remove parentheses
    paste0("(", gray_style(linkage), ")")
  })
  
  # Pattern 2: Incomplete linkages at end like (a1-, (?1-
  incomplete_pattern <- "\\(([ab?]\\d*-)$"
  result <- stringr::str_replace_all(result, incomplete_pattern, function(match) {
    linkage <- stringr::str_sub(match, 2, -1)  # Remove opening parenthesis only
    paste0("(", gray_style(linkage))
  })
  
  result
}

#' Apply Colors to IUPAC String (Monosaccharides + Gray Linkages)
#'
#' @param iupac_text Character string of IUPAC notation
#' @param mono_names Character vector of monosaccharide names to color
#' @return Character string with colored monosaccharides and gray linkages
#' @keywords internal
colorize_iupac_string <- function(iupac_text, mono_names) {
  # First, color the monosaccharides
  result <- replace_monos_with_colored(iupac_text, mono_names)
  
  # Then, make linkages gray
  result <- add_gray_linkages(result)
  
  result
}

#' Extract Base Monosaccharide Name (Without Substituents)
#'
#' @param mono A monosaccharide name (character), potentially with substituents
#' @return The base monosaccharide name without substituents
#' @keywords internal
.extract_base_mono <- function(mono) {
  # Use a simplified version of the substituent extraction logic
  # This avoids circular dependency with iupac-to-structure.R

  # Handle special cases first
  if (mono == "Neu5Ac") {
    return("Neu5Ac")
  } else if (stringr::str_starts(mono, "Neu5Ac") && nchar(mono) > 6) {
    return("Neu5Ac")
  } else if (mono == "Neu4Ac5Ac") {
    return("Neu5Ac")
  } else if (mono == "Neu4Ac5Gc") {
    return("Neu5Gc")
  }

  # For other monosaccharides, remove substituents
  subs_pattern <- stringr::str_c(available_substituents(), collapse = "|")
  single_sub_pattern <- stringr::str_glue("[1-9\\?]({subs_pattern})")

  # Find all substituents
  all_subs <- stringr::str_extract_all(mono, single_sub_pattern)[[1]]

  if (length(all_subs) > 0) {
    # Remove all substituents to get base monosaccharide
    clean_mono <- mono
    for (sub in all_subs) {
      clean_mono <- stringr::str_remove(clean_mono, stringr::fixed(sub))
    }
    return(clean_mono)
  } else {
    return(mono)
  }
}