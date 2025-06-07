#' Print method for Glycan Graphs
#'
#' This function prints information about a glycan graph,
#' including its composition and structure.
#'
#' If all linkages are unknown ("??-?"),
#' the linkage information will be omitted.
#'
#' @param x A glycan graph.
#' @param ... Ignored.
#' @param verbose A logical value. If `TRUE`, the structure of the glycan graph will be printed.
#' @param colored A logical value. If `TRUE`, monosaccharides will be colored in structure printing.
#' If `TRUE`, the structure of the glycan graph will be printed.
#' Default is `TRUE`.
#'
#' @export
print.glycan_graph <- function(x, ..., verbose = TRUE, colored = TRUE) {
  checkmate::assert_flag(verbose)
  checkmate::assert_flag(colored)
  cli::cat_line("Glycan Graph")
  print_composition(x)
  if (verbose) {
    cli::cat_line(stringr::str_dup("-", 18))
    label_getter <- function(graph) {
      if (igraph::vcount(graph) == 1) {
        label <- igraph::V(graph)$mono
        if (colored) {
          label <- add_colors(label)
        }
        if (graph$alditol) {
          label <- paste0(label, "-ol")
        }
        if (igraph::V(graph)$sub != "") {
          label <- paste0(label, "-", igraph::V(graph)$sub)
        }
        return(label)
      }
      root <- igraph::V(graph)[igraph::degree(graph, mode = "in") == 0]
      monos <- igraph::V(graph)$mono
      if (colored) {
        monos <- add_colors(monos)
      }
      if (graph$alditol) {
        monos[[root]] <- paste0(monos[[root]], "-ol")
      }
      subs <- igraph::V(graph)$sub
      monos <- dplyr::if_else(subs == "", monos, stringr::str_c(monos, subs, sep = "-"))
      if (!has_linkages(x)) return(monos)
      linkages <- purrr::map_chr(
        igraph::V(graph)[2:igraph::vcount(graph)],
        ~ igraph::incident(graph, .x, mode = "in")$linkage
      )
      linkages_str <- stringr::str_c(" (", linkages, ")", sep = "")
      linkages_str <- dplyr::if_else(is.na(linkages_str), "", linkages_str)
      linkages_str <- c("", linkages_str)
      return(stringr::str_c(monos, linkages_str))
    }
    print_structure(x, label_getter)
  }
}


print_composition <- function(graph) {
  composition <- get_composition(graph)
  comp_str <- paste(names(composition), composition, sep = ": ", collapse = ", ")
  cli::cat_line(comp_str)
}


print_structure <- function(graph, label_getter) {
  root <- igraph::V(graph)[igraph::degree(graph, mode = "in") == 0]
  labels <- label_getter(graph)

  line_last <- "\u2514\u2500"  # "└─"
  line_mid <- "\u251c\u2500"  # "├─"
  line_vertical <- "\u2502 "  # "│ "
  space <- "  "

  print_node <- function(graph, node, prefix, is_last) {
    label <- labels[node]

    if (prefix == "") {
      if (graph$anomer != "??") {
        label <- paste0(label, " (", graph$anomer, "-)")
      }
      cli::cat_line(label)
    } else {
      connector <- if (is_last) line_last else line_mid
      cli::cat_line(stringr::str_sub(prefix, 3), connector, label)
    }

    new_prefix <- if (is_last) paste0(prefix, space) else paste0(prefix, line_vertical)

    children <- igraph::neighbors(graph, node, mode = "out")
    num_children <- length(children)

    for (i in seq_along(children)) {
      child <- children[i]
      is_last_child <- i == num_children
      Recall(graph, child, new_prefix, is_last_child)
    }
  }

  print_node(graph, root, prefix = "", is_last = TRUE)
}


get_mono_color <- function(mono) {
  dplyr::case_match(
    mono,
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


add_colors <- function(monos) {
  colors <- purrr::map_chr(monos, get_mono_color)
  color_styles <- purrr::map(colors, cli::make_ansi_style)
  purrr::map2_chr(monos, color_styles, ~ .y(.x))
}
