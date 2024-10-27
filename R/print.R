#' Print method for Glycan Graphs
#'
#' This function prints information about a glycan graph,
#' including its composition and structure.
#'
#' For NE glycan graphs, if all linkages are unknown ("??-?"),
#' the linkage information will be omitted.
#'
#' @param x A glycan graph.
#' @param ... Ignored.
#' @param verbose A logical value.
#' If `TRUE`, the structure of the glycan graph will be printed.
#' Default is `TRUE`.
#'
#' @export
print.ne_glycan_graph <- function(x, ..., verbose = TRUE) {
  cli::cat_line("Glycan Graph (NE)")
  print_composition(x)
  if (verbose) {
    cli::cat_line(stringr::str_dup("-", 18))
    label_getter <- function(graph) {
      if (igraph::vcount(graph) == 1) return(igraph::V(graph)$mono)
      monos <- igraph::V(graph)$mono
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


#' @rdname print.ne_glycan_graph
#' @export
print.dn_glycan_graph <- function(x, ..., verbose = TRUE) {
  cli::cat_line("Glycan Graph (DN)")
  print_composition(x)
  if (verbose) {
    cli::cat_line(stringr::str_dup("-", 18))
    label_getter <- function(graph) {
      dplyr::if_else(
        igraph::V(graph)$type == "mono",
        igraph::V(graph)$mono,
        igraph::V(graph)$linkage
      )
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
