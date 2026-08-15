# Is the graph directed?
is_directed_graph <- function(graph) {
  igraph::is_directed(graph)
}


# Is the graph an outward tree?
is_out_tree <- function(graph) {
  igraph::is_tree(graph, mode = "out")
}


# Does the graph have these vertex attributes?
has_vertex_attrs <- function(graph, attrs) {
  all(attrs %in% igraph::vertex_attr_names(graph))
}


# Does the graph have these edge attributes?
has_edge_attrs <- function(graph, attrs) {
  if (
    getRversion() < "4.4.0" &&
      .Platform$OS.type == "windows" &&
      igraph::vcount(graph) == 1
  ) {
    # It seems like when the graph has no edges,
    # setting edge attributes will not work on Windows with R < 4.4.0.
    # So we skip this check to circumvent R CMD check.
    return(TRUE)
  }
  all(attrs %in% igraph::edge_attr_names(graph))
}


.unique_no_na <- function(x) unique(x[!is.na(x)])


# Are all monosaaccharides known?
is_known_mono <- function(monos) {
  known_monos <- c(
    .unique_no_na(monosaccharides$generic),
    monosaccharides$concrete
  )
  monos %in% known_monos
}


# Is a valid subtituent?
valid_substituent <- function(sub) {
  # Apply to each element if input is a vector
  purrr::map_lgl(sub, function(single_sub) {
    # Empty substituent is always valid
    if (single_sub == "") {
      return(TRUE)
    }

    # Split by commas to handle multiple substituents
    individual_subs <- stringr::str_split(single_sub, ",")[[1]]

    # Check if each individual substituent is valid
    pattern <- substituent_token_pattern(anchored = TRUE)

    individual_valid <- purrr::map_lgl(
      individual_subs,
      ~ stringr::str_detect(.x, pattern)
    )

    # All individual substituents must be valid
    if (!all(individual_valid)) {
      return(FALSE)
    }

    is_canonical <- individual_subs ==
      purrr::map_chr(individual_subs, normalize_substituent_token)
    if (!all(is_canonical)) {
      return(FALSE)
    }

    positions <- substituent_position_tokens(individual_subs)
    numeric_positions <- substituent_position_values(individual_subs)

    # Check if positions are sorted in ascending order
    is_sorted <- all(numeric_positions == sort(numeric_positions))

    has_assignment <- has_conflict_free_assignment(
      substituent_position_domains(individual_subs)
    )

    is_sorted && has_assignment
  })
}


# Is a valid anomer?
valid_anomer <- function(anomer) {
  stringr::str_detect(anomer, "^[ab\\?][\\d\\?]$")
}


#' Build a Linkage Regex Pattern
#'
#' Creates the shared regex pattern for glycosidic linkages.
#'
#' @param anchored Whether to anchor the pattern at the start and end.
#'
#' @returns A regex pattern string.
#'
#' @noRd
linkage_pattern <- function(anchored = TRUE) {
  checkmate::assert_flag(anchored)

  anomer_p <- "[ab\\?]"
  pos1_p <- "([12]|\\?)"
  pos2_p <- "([1-9](/[1-9])*|\\?)"
  pattern <- stringr::str_glue("{anomer_p}{pos1_p}-{pos2_p}")

  if (anchored) {
    pattern <- stringr::str_glue("^{pattern}$")
  }

  pattern
}


#' Check if Linkages are Valid
#'
#' Valid linkages are in the form of "a1-2", "b1-4", "a?-1", etc.
#' Specifically, the pattern is `xy-z`:
#' - `x`: the anomer, either "a", "b", or "?".
#' - `y`: the first position, either "1", "2" or "?".
#' - `z`: the second position, either a 1-9 digit or "?".
#' Can also be multiple positions separated by "/", e.g. "1/2/3".
#' "?" could not be used with "/".
#'
#' @param linkages A character vector of linkages.
#'
#' @returns A logical vector.
#'
#' @examples
#' # Valid linkages
#' valid_linkages(c("a1-2", "?1-4", "a?-1", "b?-?", "??-?", "a1/2-3"))
#'
#' # Invalid linkages
#' valid_linkages(c("a1-2/?", "1-4", "a/b1-2", "c1-2", "a9-1"))
#'
#' @export
valid_linkages <- function(linkages) {
  checkmate::assert_character(linkages)
  stringr::str_detect(linkages, linkage_pattern())
}


# Check if any duplicated linkage positions exist.
# The same position of one residue cannot connect to multiple other residues.
any_dup_linkage_pos <- function(
  glycan,
  linkages = igraph::edge_attr(glycan, "linkage")
) {
  if (length(linkages) < 2L) {
    return(FALSE)
  }

  endpoints <- igraph::as_edgelist(glycan, names = FALSE)
  dash <- regexpr("-", linkages, fixed = TRUE)
  positions <- substring(linkages, dash + 1L)
  check <- positions != "?" & !grepl("/", positions, fixed = TRUE)
  positions <- positions[check]

  if (length(positions) < 2L) {
    return(FALSE)
  }

  parents <- endpoints[check, 1L]
  anyDuplicated(paste(parents, positions, sep = "\r")) > 0L
}
