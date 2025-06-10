#' Parse IUPAC-condensed string to glycan structure
#'
#' Internal functions for parsing IUPAC-condensed strings into igraph objects.
#' This supports the as_glycan_structure.character method.
#'
#' @param x A single IUPAC-condensed string
#' @return An igraph object representing the glycan structure
#' @keywords internal
.parse_iupac_condensed_single <- function(x) {
  if (is.na(x) || nchar(x) == 0 || stringr::str_detect(x, "^\\s*$")) {
    rlang::abort("Cannot parse empty or NA IUPAC-condensed string.")
  }
  
  tryCatch({
    # Validate input string - no leading/trailing whitespace
    if (stringr::str_detect(x, "^\\s+|\\s+$")) {
      rlang::abort("IUPAC-condensed string cannot have leading or trailing whitespace")
    }
    
    # Validate no internal whitespace
    if (stringr::str_detect(x, "\\s")) {
      rlang::abort("IUPAC-condensed string cannot contain whitespace")
    }
    
    # Validate proper bracket matching
    if (!.validate_brackets(x)) {
      rlang::abort("Malformed brackets in IUPAC-condensed string")
    }
    anomer <- .extract_anomer(x)  # may be NA
    if (!is.na(anomer)) {
      x <- stringr::str_sub(x, 1, -stringr::str_length(anomer)-3)
    }
    
    alditol <- .extract_alditol(x)
    if (alditol) {
      x <- stringr::str_replace(x, "-ol", "")
    }
    
    tokens <- .tokenize_iupac(x)
    
    # Deal with missing anomer
    first_mono_sub_res <- .extract_substituent(tokens[[1]])
    if (is.na(anomer)) {
      anomer <- local({
        anomer_pos <- .decide_anomer_pos(first_mono_sub_res[["mono"]])
        paste0("?", anomer_pos)
      })
    }
    
    # Create a new graph and add the first node
    graph <- igraph::make_empty_graph()
    graph <- igraph::add_vertices(
      graph, 1, name = "1",
      mono = first_mono_sub_res[["mono"]],
      sub = first_mono_sub_res[["sub"]]
    )
    
    if (length(tokens) == 1) {
      graph <- igraph::set_edge_attr(graph, "linkage", value = character(0))
      graph$anomer <- anomer
      graph$alditol <- alditol
      return(graph)
    }
    
    # Iterate over the tokens
    node_stack <- rstackdeque::rstack()
    node_stack <- rstackdeque::insert_top(node_stack, 1)
    current_node_id <- 1
    for (token in tokens[2:length(tokens)]) {
      if (token == "[") {
        node_stack <- rstackdeque::insert_top(node_stack, current_node_id)
      } else if (token == "]") {
        current_node_id <- rstackdeque::peek_top(node_stack)
        node_stack <- rstackdeque::without_top(node_stack)
      } else {
        parsed_token <- .parse_token(token)
        new_node_id <- igraph::vcount(graph) + 1
        graph <- igraph::add_vertices(
          graph, 1,
          name = as.character(new_node_id),
          mono = parsed_token[["mono"]],
          sub = parsed_token[["sub"]]
        )
        graph <- igraph::add_edges(
          graph,
          c(current_node_id, new_node_id),
          linkage = parsed_token[["linkage"]]
        )
        current_node_id <- new_node_id
      }
    }
    
    graph$anomer <- anomer
    graph$alditol <- alditol
    return(graph)
  }, error = function(e) {
    rlang::abort(c(
      "Could not parse IUPAC-condensed string: {.val {x}}",
      "i" = conditionMessage(e)
    ))
  })
}

# Validate bracket matching
.validate_brackets <- function(x) {
  # Check for balanced brackets
  open_count <- stringr::str_count(x, "\\[")
  close_count <- stringr::str_count(x, "\\]")
  
  if (open_count != close_count) {
    return(FALSE)
  }
  
  # Check for proper nesting
  depth <- 0
  chars <- stringr::str_split(x, "")[[1]]
  for (char in chars) {
    if (char == "[") {
      depth <- depth + 1
    } else if (char == "]") {
      depth <- depth - 1
      if (depth < 0) {
        return(FALSE)
      }
    }
  }
  
  return(depth == 0)
}

# Extract anomer from IUPAC condensed string
.extract_anomer <- function(iupac) {
  # e.g. "Neu5Ac(a2-" -> "a2"  (ending anomer specification)
  # e.g. "Neu5Ac" -> NA
  p <- "\\(([ab\\?][12\\?])-$"
  if (stringr::str_detect(iupac, p)) {
    # Check if it's a complete IUPAC string ending with anomer specification
    # This should be allowed for standalone monosaccharides like "Neu5Ac(a2-"
    stringr::str_extract(iupac, p, group = 1)
  } else {
    NA_character_
  }
}

# Extract alditol from IUPAC condensed string
.extract_alditol <- function(iupac) {
  stringr::str_detect(iupac, "-ol")
}

# Tokenize IUPAC condensed string
.tokenize_iupac <- function(iupac) {
  # anomer is either "a", "b", or "?".
  anomer_p <- "[ab\\?]"
  
  # pos1 is either "1", "2", or "?".
  pos1_p <- "([12\\?])"
  
  # pos2 is "1"-"9" (followed by any number of "/1-9"), or "?"
  pos2_p <- "([1-9](/[1-9])*|\\?)"
  
  # Linkage pattern
  linkage_pattern <- stringr::str_glue("{anomer_p}{pos1_p}-{pos2_p}")
  
  # Monosaccharide name pattern (including potential substituents)
  # Allow letters, digits, and ? for substituents like "Man?S", "Gal6S", etc.
  mono_pattern <- "[A-Za-z][A-Za-z0-9\\?]*"
  mono_linkage_pattern <- stringr::str_glue("{mono_pattern}(\\({linkage_pattern}\\))?")
  
  # The pattern is either a monosaccharide name or a bracket
  pattern <- paste(mono_linkage_pattern, "\\[", "\\]", sep = "|")
  
  tokens <- stringr::str_extract_all(iupac, pattern)[[1]]
  
  # Check if we extracted the full string
  extracted_string <- paste(tokens, collapse = "")
  if (extracted_string != iupac) {
    rlang::abort("Invalid characters or format in IUPAC-condensed string")
  }
  
  tokens <- stringr::str_replace(tokens, "\\[", "TEMP_LEFT")
  tokens <- stringr::str_replace(tokens, "\\]", "TEMP_RIGHT")
  tokens <- stringr::str_replace(tokens, "TEMP_LEFT", "\\]")
  tokens <- stringr::str_replace(tokens, "TEMP_RIGHT", "\\[")
  
  # Reverse the tokens to make the first monosaccharide the reducing end
  rev(tokens)
}

# Parse a single token
.parse_token <- function(token) {
  left_bracket_pos <- stringr::str_locate(token, "\\(")[1]
  
  if (is.na(left_bracket_pos)) {
    rlang::abort("Missing linkage information in token")
  }
  
  mono <- stringr::str_sub(token, 1, left_bracket_pos - 1)
  mono_sub_res <- .extract_substituent(mono)
  mono <- mono_sub_res[["mono"]]
  sub <- mono_sub_res[["sub"]]
  linkage <- stringr::str_sub(token, left_bracket_pos + 1, -2)
  
  # Validate linkage format
  if (!.validate_linkage(linkage)) {
    rlang::abort(paste0("Invalid linkage format: ", linkage))
  }
  
  c(mono = mono, sub = sub, linkage = linkage)
}

# Validate linkage format
.validate_linkage <- function(linkage) {
  # Empty linkage
  if (nchar(linkage) == 0) {
    return(FALSE)
  }
  
  # Pattern for valid linkage: anomer + pos1 + - + pos2
  # anomer: a, b, or ?
  # pos1: 1, 2, or ?
  # pos2: 1-9 (optionally followed by /1-9 patterns), or ?
  valid_pattern <- "^[ab\\?][12\\?]-([1-9](/[1-9])*|\\?)$"
  
  return(stringr::str_detect(linkage, valid_pattern))
}

# Extract substituent from monosaccharide name
.extract_substituent <- function(mono) {
  subs_pattern <- stringr::str_c(available_substituents(), collapse = "|")
  subs_pattern <- stringr::str_glue("[1-9\\?]({subs_pattern})$")  # Changed to exclude 0
  
  result <- if (mono == "Neu5Ac") {
    # "Neu5Ac" is special that it satisfies the regex below,
    # but should not be split.
    c(mono = mono, sub = "")
  } else if (mono == "Neu4Ac5Ac") {
    c(mono = "Neu5Ac", sub = "4Ac")
  } else if (mono == "Neu4Ac5Gc") {
    c(mono = "Neu5Gc", sub = "4Ac")
  } else if (stringr::str_detect(mono, subs_pattern)) {
    sub_loc <- stringr::str_locate(mono, subs_pattern)[1]
    sub <- stringr::str_sub(mono, sub_loc, -1)
    mono <- stringr::str_sub(mono, 1, sub_loc - 1)
    c(mono = mono, sub = sub)
  } else {
    c(mono = mono, sub = "")
  }
  
  # Validate that the monosaccharide is known
  if (!is_known_monosaccharide(result[["mono"]])) {
    rlang::abort(paste0("Unknown monosaccharide: ", result[["mono"]]))
  }
  
  result
}

# Decide anomer position for known monosaccharides
.decide_anomer_pos <- function(mono) {
  anomer_on_pos2 <- c(
    "Neu5Ac", "Neu5Gc", "Neu", "Kdn", "Pse", "Leg", "Aci",
    "4eLeg", "Kdo", "Dha", "Fru", "Tag", "Sor", "Psi"
  )
  dplyr::if_else(mono %in% anomer_on_pos2, "2", "1")
} 