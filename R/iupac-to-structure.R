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
    cli::cli_abort("Cannot parse empty or NA IUPAC-condensed string.")
  }
  
  tryCatch({
    # Validate input string - no leading/trailing whitespace
    if (stringr::str_detect(x, "^\\s+|\\s+$")) {
      cli::cli_abort("IUPAC-condensed string cannot have leading or trailing whitespace")
    }
    
    # Validate no internal whitespace
    if (stringr::str_detect(x, "\\s")) {
      cli::cli_abort("IUPAC-condensed string cannot contain whitespace")
    }
    
    # Validate proper bracket matching
    if (!.validate_brackets(x)) {
      cli::cli_abort("Malformed brackets in IUPAC-condensed string")
    }
    anomer <- .extract_anomer(x)
    x <- stringr::str_sub(x, 1, -stringr::str_length(anomer)-3)

    alditol <- .extract_alditol(x)
    if (alditol) {
      x <- stringr::str_replace(x, "-ol", "")
    }

    tokens <- .tokenize_iupac(x)

    # Require anomer information - no longer auto-supplement
    first_mono_sub_res <- .extract_substituent(tokens[[1]])
    
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
    cli::cli_abort(c(
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
  p <- "\\(([ab\\?][12\\?])-$"
  if (stringr::str_detect(iupac, p)) {
    # Check if it's a complete IUPAC string ending with anomer specification
    # This should be allowed for standalone monosaccharides like "Neu5Ac(a2-"
    stringr::str_extract(iupac, p, group = 1)
  } else {
    cli::cli_abort(c(
      "Can't extract anomer information.",
      "i" = "Anomer information is required for the reducing-end monosaccharide.",
      "i" = "For example, use 'Man(a1-' instead of 'Man'."
    ))
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
  # Allow letters, digits, and ? for substituents like "Man?S", "Glc3Me6S", etc.
  # Substituents are directly concatenated in IUPAC format, no commas
  mono_pattern <- "[A-Za-z][A-Za-z0-9\\?]*"
  mono_linkage_pattern <- stringr::str_glue("{mono_pattern}(\\({linkage_pattern}\\))?")
  
  # The pattern is either a monosaccharide name or a bracket
  pattern <- paste(mono_linkage_pattern, "\\[", "\\]", sep = "|")
  
  tokens <- stringr::str_extract_all(iupac, pattern)[[1]]
  
  # Check if we extracted the full string
  extracted_string <- paste(tokens, collapse = "")
  if (extracted_string != iupac) {
    cli::cli_abort("Invalid characters or format in IUPAC-condensed string")
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
    cli::cli_abort("Missing linkage information in token")
  }
  
  mono <- stringr::str_sub(token, 1, left_bracket_pos - 1)
  mono_sub_res <- .extract_substituent(mono)
  mono <- mono_sub_res[["mono"]]
  sub <- mono_sub_res[["sub"]]
  linkage <- stringr::str_sub(token, left_bracket_pos + 1, -2)
  
  # Validate linkage format
  if (!.validate_linkage(linkage)) {
    cli::cli_abort(paste0("Invalid linkage format: ", linkage))
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
  single_sub_pattern <- stringr::str_glue("[1-9\\?]({subs_pattern})")  # Pattern for a single substituent
  
  result <- if (mono == "Neu5Ac") {
    # "Neu5Ac" is special that it satisfies the regex below,
    # but should not be split.
    c(mono = mono, sub = "")
  } else if (stringr::str_starts(mono, "Neu5Ac") && nchar(mono) > 6) {
    # Handle "Neu5Ac" with additional substituents, e.g., "Neu5Ac9Ac"
    # Extract everything after "Neu5Ac" as substituents
    sub_part <- stringr::str_sub(mono, 7, -1)  # Everything after "Neu5Ac"
    
    # Parse the substituent part
    all_subs <- stringr::str_extract_all(sub_part, single_sub_pattern)[[1]]
    
    if (length(all_subs) > 0) {
      # Sort substituents by position
      positions <- purrr::map_chr(all_subs, ~ stringr::str_extract(.x, "^[1-9\\?]"))
      numeric_positions <- purrr::map_dbl(positions, function(pos) {
        if (pos == "?") Inf else as.numeric(pos)
      })
      
      sorted_indices <- order(numeric_positions)
      sorted_subs <- all_subs[sorted_indices]
      
      sub_string <- stringr::str_c(sorted_subs, collapse = ",")
      c(mono = "Neu5Ac", sub = sub_string)
    } else {
      c(mono = "Neu5Ac", sub = "")
    }
  } else if (mono == "Neu4Ac5Ac") {
    # Special case: Neu4Ac5Ac -> Neu5Ac + 4Ac (5Ac is absorbed into the mono name)
    c(mono = "Neu5Ac", sub = "4Ac")
  } else if (mono == "Neu4Ac5Gc") {
    # Special case: Neu4Ac5Gc -> Neu5Gc + 4Ac (5Gc is absorbed into the mono name)
    c(mono = "Neu5Gc", sub = "4Ac")
  } else {
    # Try to find all substituents in the monosaccharide name
    # Substituents are directly concatenated in IUPAC, e.g., "Glc3Me6S"
    all_subs <- stringr::str_extract_all(mono, single_sub_pattern)[[1]]
    
    if (length(all_subs) > 0) {
      # Remove all substituents from the mono name to get the base monosaccharide
      clean_mono <- mono
      for (sub in all_subs) {
        clean_mono <- stringr::str_remove(clean_mono, stringr::fixed(sub))
      }
      
      # Sort substituents by position
      # Extract positions for sorting
      positions <- purrr::map_chr(all_subs, ~ stringr::str_extract(.x, "^[1-9\\?]"))
      numeric_positions <- purrr::map_dbl(positions, function(pos) {
        if (pos == "?") Inf else as.numeric(pos)
      })
      
      # Sort substituents by position
      sorted_indices <- order(numeric_positions)
      sorted_subs <- all_subs[sorted_indices]
      
      # Combine into comma-separated string for internal storage
      sub_string <- if (length(sorted_subs) > 0) {
        stringr::str_c(sorted_subs, collapse = ",")
      } else {
        ""
      }
      
      c(mono = clean_mono, sub = sub_string)
    } else {
      c(mono = mono, sub = "")
    }
  }
  
  # Validate that the monosaccharide is known
  if (!is_known_monosaccharide(result[["mono"]])) {
    cli::cli_abort(paste0("Unknown monosaccharide: ", result[["mono"]]))
  }
  
  result
}

