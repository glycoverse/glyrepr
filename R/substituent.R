#' Available Substituents
#'
#' Get the available substituents for monosaccharides.
#'
#' @returns A character vector.
#'
#' @examples
#' available_substituents()
#'
#' @export
available_substituents <- function() {
  c("Me", "Ac", "NAc", "NGc", "P", "S", "Pyr", "PC", "PPEtn", "PEtn", "N", "Gc")
}

#' Build a Substituent Name Regex Pattern
#'
#' @param longest_first Whether to sort longer substituent names before shorter
#'   names in the regex alternation.
#'
#' @returns A regex alternation pattern for known substituent names.
#'
#' @noRd
substituent_name_pattern <- function(longest_first = FALSE) {
  checkmate::assert_flag(longest_first)

  substituents <- available_substituents()
  if (longest_first) {
    substituents <- substituents[order(nchar(substituents), decreasing = TRUE)]
  }

  stringr::str_c(stringr::str_escape(substituents), collapse = "|")
}

#' Build a Substituent Position Regex Pattern
#'
#' @returns A regex alternation pattern for one substituent position, ambiguous
#'   slash-separated positions, or an unknown position.
#'
#' @noRd
substituent_position_pattern <- function() {
  "(?:[1-9](?:/[1-9])*|\\?)"
}

#' Build a Substituent Token Regex Pattern
#'
#' @param longest_first Whether to sort longer substituent names before shorter
#'   names in the regex alternation.
#' @param anchored Whether to anchor the pattern at the start and end.
#'
#' @returns A regex pattern for one complete substituent token.
#'
#' @noRd
substituent_token_pattern <- function(
  longest_first = FALSE,
  anchored = FALSE
) {
  checkmate::assert_flag(longest_first)
  checkmate::assert_flag(anchored)

  position <- substituent_position_pattern()
  name <- substituent_name_pattern(longest_first = longest_first)
  pattern <- stringr::str_glue("{position}(?:{name})")

  if (anchored) {
    pattern <- stringr::str_glue("^{pattern}$")
  }

  pattern
}

#' Normalize Substituent String
#'
#' Takes a substituent string (potentially with multiple substituents) and
#' returns a normalized string with substituents sorted by position.
#'
#' @param sub A character string representing substituents, e.g., "4Ac,3Me" or "6S"
#' @returns A character string with substituents sorted by position, e.g., "3Me,4Ac"
#'
#' @examples
#' normalize_substituents("4Ac,3Me")  # Returns "3Me,4Ac"
#' normalize_substituents("6S")       # Returns "6S"
#' normalize_substituents("")         # Returns ""
#' @noRd
normalize_substituents <- function(sub) {
  checkmate::assert_character(sub, len = 1)

  if (sub == "") {
    return("")
  }

  # Split by commas
  individual_subs <- stringr::str_split(sub, ",")[[1]]

  # Remove any empty strings (in case of double commas)
  individual_subs <- individual_subs[individual_subs != ""]

  if (length(individual_subs) == 0) {
    return("")
  }

  collapse_substituent_tokens(individual_subs)
}

#' Get Substituent Position Tokens
#'
#' Extracts the leading substituent positions from substituent tokens.
#'
#' @param subs A character vector of substituent tokens.
#'
#' @returns A character vector of position tokens.
#'
#' @noRd
substituent_position_tokens <- function(subs) {
  pattern <- stringr::str_glue("^{substituent_position_pattern()}")

  purrr::map_chr(
    subs,
    ~ stringr::str_extract(.x, pattern)
  )
}

#' Normalize One Substituent Token
#'
#' Sorts and deduplicates slash-separated candidate positions while preserving
#' the substituent name.
#'
#' @param sub One valid substituent token.
#'
#' @returns A canonical substituent token.
#'
#' @noRd
normalize_substituent_token <- function(sub) {
  position <- substituent_position_tokens(sub)
  if (position == "?") {
    return(sub)
  }

  candidates <- stringr::str_split(position, stringr::fixed("/"))[[1]]
  candidates <- sort(unique(as.integer(candidates)))
  name <- stringr::str_sub(sub, stringr::str_length(position) + 1L)

  paste0(stringr::str_c(candidates, collapse = "/"), name)
}

#' Get Candidate Position Sets for Substituent Tokens
#'
#' @param subs A character vector of valid substituent tokens.
#'
#' @returns A list of character vectors. Unknown-position substituents are
#'   omitted because they do not constrain known position assignments.
#'
#' @noRd
substituent_position_domains <- function(subs) {
  positions <- substituent_position_tokens(subs)
  positions <- positions[positions != "?"]
  lapply(positions, stringr::str_split_1, pattern = stringr::fixed("/"))
}

#' Check for a Conflict-Free Assignment
#'
#' @param domains A list of candidate value vectors.
#'
#' @returns `TRUE` when each domain can be assigned a distinct value.
#'
#' @noRd
has_conflict_free_assignment <- function(domains) {
  if (length(domains) == 0) {
    return(TRUE)
  }
  if (any(lengths(domains) == 0)) {
    return(FALSE)
  }

  domains <- domains[order(lengths(domains))]
  assign_value <- function(index, used) {
    if (index > length(domains)) {
      return(TRUE)
    }

    available <- setdiff(domains[[index]], used)
    for (value in available) {
      if (assign_value(index + 1L, c(used, value))) {
        return(TRUE)
      }
    }
    FALSE
  }

  assign_value(1L, character())
}

#' Get Numeric Substituent Position Values
#'
#' Converts substituent position tokens to numeric values, using `Inf` for
#' unknown positions so they sort after known positions.
#'
#' @param subs A character vector of substituent tokens.
#'
#' @returns A numeric vector of position values.
#'
#' @noRd
substituent_position_values <- function(subs) {
  positions <- substituent_position_tokens(subs)

  purrr::map_dbl(positions, function(pos) {
    if (pos == "?") {
      Inf
    } else {
      as.numeric(stringr::str_extract(pos, "^[1-9]"))
    }
  })
}

#' Sort Substituent Tokens by Position
#'
#' Sorts substituent tokens by their leading position, placing unknown positions
#' after known positions.
#'
#' @param subs A character vector of substituent tokens.
#'
#' @returns A character vector sorted by position.
#'
#' @noRd
sort_substituent_tokens <- function(subs) {
  if (length(subs) == 0) {
    return(character(0))
  }

  subs[order(substituent_position_values(subs))]
}

#' Collapse Sorted Substituent Tokens
#'
#' Sorts substituent tokens by position and collapses them to the comma-separated
#' representation used in glycan structures.
#'
#' @param subs A character vector of substituent tokens.
#'
#' @returns A comma-separated substituent string.
#'
#' @noRd
collapse_substituent_tokens <- function(subs) {
  if (length(subs) == 0) {
    return("")
  }

  subs <- purrr::map_chr(subs, normalize_substituent_token)
  stringr::str_c(sort_substituent_tokens(subs), collapse = ",")
}

#' Remove All Substituents from a Glycan
#'
#' This function replaces all vertex substituents in a glycan structure with
#' empty strings and removes unresolved floating substituents.
#'
#' @param glycan A glyrepr_structure vector or a glycan `igraph`.
#'
#' @returns An object of the same representation as `glycan` with all
#'   substituents removed. Graph input retains its vertex IDs and order.
#'
#' @examples
#' (glycan <- o_glycan_core_1())
#' remove_substituents(glycan)
#'
#' @export
remove_substituents <- function(glycan) {
  if (inherits(glycan, "igraph")) {
    return(.remove_substituents_single(glycan))
  }

  if (!is_glycan_structure(glycan)) {
    cli::cli_abort(c(
      "Input must be a glyrepr_structure vector.",
      "i" = "Use `glycan_structure()` to create a glyrepr_structure from igraph objects."
    ))
  }

  smap_structure(glycan, .remove_substituents_single)
}

# Internal function to remove substituents from a single igraph
.remove_substituents_single <- function(glycan) {
  glycan <- igraph::set_vertex_attr(glycan, "sub", value = "")
  delete_floating_substituents_attr(glycan)
}
