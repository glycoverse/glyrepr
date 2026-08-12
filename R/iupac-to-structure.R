#' Parse IUPAC-condensed string to glycan structure
#'
#' Internal functions for parsing IUPAC-condensed strings into igraph objects.
#' This supports the as_glycan_structure.character method.
#'
#' @param x A single IUPAC-condensed string
#' @return An igraph object representing the glycan structure
#' @keywords internal
.parse_iupac_condensed_single <- function(x) {
  x <- stringr::str_replace_all(
    x,
    stringr::fixed("(?-?)"),
    "(??-?)"
  )
  x <- stringr::str_replace_all(
    x,
    "\\(([ab\\?][12\\?])-(?:[1-9]/)*\\?(?:/[1-9])*\\)",
    "(\\1-?)"
  )

  if (!isTRUE(startsWith(x, "{"))) {
    return(.parse_iupac_tree_single(x))
  }

  floating <- split_floating_iupac(x)

  main <- .parse_iupac_tree_single(floating$main)
  igraph::V(main)$source_index <- rev(seq_len(igraph::vcount(main)))
  main <- validate_glycan_graph(main)
  main <- canonicalize_glycan_graph(main)
  source_to_canonical <- match(
    seq_len(igraph::vcount(main)),
    igraph::V(main)$source_index
  )
  floating$parts <- purrr::map(floating$parts, function(part) {
    if (
      length(part$parents) > 0 &&
        any(part$parents > length(source_to_canonical))
    ) {
      cli::cli_abort(c(
        "Floating part parent index is outside the main glycan.",
        "i" = "The main glycan has {length(source_to_canonical)} node{?s}."
      ))
    }
    part$parents <- source_to_canonical[part$parents]
    part
  })
  main <- igraph::delete_vertex_attr(main, "source_index")

  parsed_parts <- purrr::map(
    floating$parts,
    parse_floating_iupac_part,
    main_size = igraph::vcount(main)
  )
  combine_floating_iupac_graphs(main, parsed_parts)
}

split_floating_iupac <- function(x) {
  if (is.na(x) || nchar(x) == 0 || stringr::str_detect(x, "^\\s*$")) {
    cli::cli_abort("Cannot parse empty or NA IUPAC-condensed string.")
  }
  if (!stringr::str_starts(x, stringr::fixed("{"))) {
    return(list(parts = list(), main = x))
  }
  if (stringr::str_detect(x, "^\\s+|\\s+$")) {
    cli::cli_abort(
      "IUPAC-condensed string cannot have leading or trailing whitespace"
    )
  }
  if (stringr::str_detect(x, "\\s")) {
    cli::cli_abort("IUPAC-condensed string cannot contain whitespace")
  }

  parts <- list()
  remainder <- x
  while (stringr::str_starts(remainder, stringr::fixed("{"))) {
    chars <- stringr::str_split(remainder, "")[[1]]
    closing_positions <- which(chars == "}")
    if (length(closing_positions) == 0) {
      cli::cli_abort("Malformed floating part in IUPAC-condensed string.")
    }
    closing <- closing_positions[[1]]
    if (closing == 1) {
      cli::cli_abort("Malformed floating part in IUPAC-condensed string.")
    }
    if (any(chars[seq_len(closing)] == "{" & seq_len(closing) != 1)) {
      cli::cli_abort("Floating parts cannot be nested.")
    }

    content <- stringr::str_sub(remainder, 2, closing - 1)
    parts[[length(parts) + 1]] <- split_floating_iupac_part(content)
    remainder <- stringr::str_sub(remainder, closing + 1)
  }

  if (!nzchar(remainder)) {
    cli::cli_abort("A floating glycan structure must have a main glycan.")
  }
  if (stringr::str_detect(remainder, "[{}]")) {
    cli::cli_abort(
      "Floating parts must precede the main IUPAC-condensed structure."
    )
  }

  list(parts = parts, main = remainder)
}

split_floating_iupac_part <- function(content) {
  pipe_count <- stringr::str_count(content, stringr::fixed("|"))
  if (pipe_count > 1) {
    cli::cli_abort(
      "A floating part can contain at most one parent-index separator."
    )
  }

  fields <- stringr::str_split(content, stringr::fixed("|"))[[1]]
  sequence <- fields[[1]]
  if (!nzchar(sequence)) {
    cli::cli_abort("A floating part cannot be empty.")
  }

  parents <- integer()
  if (length(fields) == 2) {
    parent_text <- fields[[2]]
    if (!stringr::str_detect(parent_text, "^[1-9][0-9]*(,[1-9][0-9]*)*$")) {
      cli::cli_abort(
        "Floating part parents must be comma-separated positive node indices."
      )
    }
    parent_tokens <- stringr::str_split(parent_text, ",")[[1]]
    parent_values <- suppressWarnings(as.double(parent_tokens))
    if (
      any(!is.finite(parent_values)) ||
        any(parent_values > .Machine$integer.max)
    ) {
      cli::cli_abort(
        "Floating part parent indices exceed the supported integer range."
      )
    }
    parents <- as.integer(parent_values)
    if (anyDuplicated(parents) > 0) {
      cli::cli_abort("Floating part parent indices must be unique.")
    }
  }

  list(sequence = sequence, parents = parents)
}

parse_floating_iupac_part <- function(part, main_size) {
  linkage_match <- stringr::str_match(
    part$sequence,
    paste0("\\((", linkage_pattern(anchored = FALSE), ")\\)$")
  )
  linkage <- linkage_match[[1, 2]]
  if (is.na(linkage)) {
    cli::cli_abort(
      "A floating part must end with its linkage to the main glycan."
    )
  }
  if (length(part$parents) > 0 && any(part$parents > main_size)) {
    cli::cli_abort(c(
      "Floating part parent index is outside the main glycan.",
      "i" = "The main glycan has {main_size} node{?s}."
    ))
  }

  sequence <- stringr::str_remove(
    part$sequence,
    paste0("\\(", linkage_pattern(anchored = FALSE), "\\)$")
  )
  donor <- stringr::str_sub(linkage, 1, 2)
  graph <- .parse_iupac_tree_single(paste0(sequence, "(", donor, "-"))

  list(
    graph = graph,
    linkage = linkage,
    parents = sort(part$parents)
  )
}

combine_floating_iupac_graphs <- function(main, parts) {
  graphs <- c(list(main), purrr::map(parts, "graph"))
  sizes <- purrr::map_int(graphs, igraph::vcount)
  offsets <- cumsum(c(0L, sizes[-length(sizes)]))

  graph <- igraph::make_empty_graph(sum(sizes), directed = TRUE)
  igraph::V(graph)$name <- as.character(seq_len(igraph::vcount(graph)))
  igraph::V(graph)$mono <- unlist(
    purrr::map(graphs, ~ igraph::V(.x)$mono),
    use.names = FALSE
  )
  igraph::V(graph)$sub <- unlist(
    purrr::map(graphs, ~ igraph::V(.x)$sub),
    use.names = FALSE
  )

  edge_vectors <- purrr::map2(
    graphs,
    offsets,
    function(component, offset) {
      edges <- igraph::as_edgelist(component, names = FALSE)
      if (length(edges) == 0) {
        return(integer())
      }
      as.integer(t(edges + offset))
    }
  )
  edge_vector <- unlist(edge_vectors, use.names = FALSE)
  if (length(edge_vector) > 0) {
    graph <- igraph::add_edges(graph, edge_vector)
  }
  igraph::E(graph)$linkage <- unlist(
    purrr::map(graphs, ~ igraph::E(.x)$linkage),
    use.names = FALSE
  )
  graph$anomer <- main$anomer

  metadata <- purrr::map2(
    parts,
    offsets[-1],
    function(part, offset) {
      cache <- build_seq_cache(part$graph)
      list(
        root = as.integer(offset + cache$root),
        nodes = as.integer(offset + seq_len(igraph::vcount(part$graph))),
        linkage = part$linkage,
        parents = part$parents
      )
    }
  )
  set_floating_parts_attr(graph, metadata)
}

.parse_iupac_tree_single <- function(x) {
  if (is.na(x) || nchar(x) == 0 || stringr::str_detect(x, "^\\s*$")) {
    cli::cli_abort("Cannot parse empty or NA IUPAC-condensed string.")
  }

  tryCatch(
    {
      # Validate input string - no leading/trailing whitespace
      if (stringr::str_detect(x, "^\\s+|\\s+$")) {
        cli::cli_abort(
          "IUPAC-condensed string cannot have leading or trailing whitespace"
        )
      }

      # Validate no internal whitespace
      if (stringr::str_detect(x, "\\s")) {
        cli::cli_abort("IUPAC-condensed string cannot contain whitespace")
      }

      # Validate proper bracket matching
      if (!.validate_brackets(x)) {
        cli::cli_abort("Malformed brackets in IUPAC-condensed string")
      }
      x <- .infer_reducing_end_anomer(x)
      anomer <- .extract_anomer(x)
      x <- stringr::str_sub(x, 1, -stringr::str_length(anomer) - 3)

      tokens <- .tokenize_iupac(x)

      # Require anomer information - no longer auto-supplement
      first_mono_sub_res <- .extract_substituent(tokens[[1]])

      # Create a new graph and add the first node
      graph <- igraph::make_empty_graph()
      graph <- igraph::add_vertices(
        graph,
        1,
        name = "1",
        mono = first_mono_sub_res[["mono"]],
        sub = first_mono_sub_res[["sub"]]
      )

      if (length(tokens) == 1) {
        graph <- igraph::set_edge_attr(graph, "linkage", value = character(0))
        graph$anomer <- anomer
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
            graph,
            1,
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
      return(graph)
    },
    error = function(e) {
      cli::cli_abort(c(
        "Could not parse IUPAC-condensed string: {.val {x}}",
        "i" = conditionMessage(e)
      ))
    }
  )
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

# Infer a missing reducing-end anomer position
.infer_reducing_end_anomer <- function(iupac) {
  if (stringr::str_detect(iupac, "\\(([ab\\?][12\\?])-$")) {
    return(iupac)
  }

  reducing_end <- .tokenize_iupac(iupac)[[1]]
  mono <- .extract_substituent(reducing_end)[["mono"]]
  paste0(iupac, "(?", infer_anomer_pos(mono), "-")
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


# Tokenize IUPAC condensed string
.tokenize_iupac <- function(iupac) {
  # Monosaccharide name pattern (including potential substituents)
  # Allow known names that start with digits, e.g. "6dGul" and "4eLeg".
  # Allow letters, digits, ?, and / for substituents like "Man?S",
  # "Glc3Me6S", and "Gal4/6S".
  # Substituents are directly concatenated in IUPAC format, no commas
  mono_pattern <- "([A-Za-z]|[0-9][A-Za-z])[A-Za-z0-9\\?/]*"
  mono_linkage_pattern <- stringr::str_glue(
    "{mono_pattern}(\\({linkage_pattern(anchored = FALSE)}\\))?"
  )

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
  if (!valid_linkages(linkage)) {
    cli::cli_abort(paste0("Invalid linkage format: ", linkage))
  }

  c(mono = mono, sub = sub, linkage = linkage)
}

# Extract substituent from monosaccharide name
.extract_substituent <- function(mono) {
  single_sub_pattern <- substituent_token_pattern(longest_first = TRUE)

  # Handle different types of monosaccharides
  result <- if (stringr::str_starts(mono, "Neu")) {
    # Handle all Neu-based monosaccharides
    .handle_neu_monosaccharide(mono, single_sub_pattern)
  } else {
    # Handle non-Neu monosaccharides
    .handle_general_monosaccharide(mono, single_sub_pattern)
  }

  # Validate that the monosaccharide is known
  if (!is_known_monosaccharide(result[["mono"]])) {
    cli::cli_abort(paste0("Unknown monosaccharide: ", result[["mono"]]))
  }

  result
}

# Handle Neu-based monosaccharides with substituents
# This function determines the correct base monosaccharide (Neu5Ac, Neu5Gc, or Neu)
# based on the presence of 5Ac or 5Gc substituents
.handle_neu_monosaccharide <- function(mono, single_sub_pattern) {
  # Check for conflicting 5Ac5Gc pattern
  if (.has_neu_marker(mono, "5Gc") && .has_neu_marker(mono, "5Ac")) {
    cli::cli_abort(
      "Monosaccharide cannot have both 5Ac and 5Gc substituents: {mono}"
    )
  }

  # Handle all Neu variants containing 5Ac
  if (.has_neu_marker(mono, "5Ac")) {
    base_mono <- if (stringr::str_starts(mono, "Neuf")) {
      "Neuf5Ac"
    } else {
      "Neu5Ac"
    }
    return(.handle_neu5ac_variant(mono, single_sub_pattern, base_mono))
  }

  # Handle all Neu variants containing 5Gc
  if (.has_neu_marker(mono, "5Gc")) {
    base_mono <- if (stringr::str_starts(mono, "Neuf")) {
      "Neuf5Gc"
    } else {
      "Neu5Gc"
    }
    return(.handle_neu5gc_variant(mono, single_sub_pattern, base_mono))
  }

  # Handle other Neu-based monosaccharides (no 5Ac or 5Gc)
  .handle_general_monosaccharide(mono, single_sub_pattern)
}

.neu_marker_pattern <- function(marker) {
  stringr::str_glue("(?<![0-9/]){stringr::str_escape(marker)}")
}

.has_neu_marker <- function(mono, marker) {
  stringr::str_detect(mono, .neu_marker_pattern(marker))
}

# Handle general (non-Neu) monosaccharides with substituents
.handle_general_monosaccharide <- function(mono, single_sub_pattern) {
  # Try to find all substituents in the monosaccharide name
  # Substituents are directly concatenated in IUPAC, e.g., "Glc3Me6S"
  all_subs <- stringr::str_extract_all(mono, single_sub_pattern)[[1]]

  if (length(all_subs) > 0) {
    # Remove all substituents from the mono name to get the base monosaccharide
    clean_mono <- mono
    for (sub in all_subs) {
      clean_mono <- stringr::str_remove(clean_mono, stringr::fixed(sub))
    }

    sub_string <- collapse_substituent_tokens(all_subs)

    c(mono = clean_mono, sub = sub_string)
  } else {
    c(mono = mono, sub = "")
  }
}

# Handle all Neu variants containing 5Ac
.handle_neu5ac_variant <- function(
  mono,
  single_sub_pattern,
  base_mono = "Neu5Ac"
) {
  .handle_neu_with_marker(mono, single_sub_pattern, "5Ac", base_mono)
}

# Handle all Neu variants containing 5Gc
.handle_neu5gc_variant <- function(
  mono,
  single_sub_pattern,
  base_mono = "Neu5Gc"
) {
  .handle_neu_with_marker(mono, single_sub_pattern, "5Gc", base_mono)
}

# Generic helper function for handling Neu variants with specific markers
.handle_neu_with_marker <- function(
  mono,
  single_sub_pattern,
  marker,
  base_mono
) {
  # Remove the marker from the monosaccharide name to get the remaining part
  mono_without_marker <- stringr::str_remove(
    mono,
    .neu_marker_pattern(marker)
  )

  # Extract all substituents from the remaining part using normal logic
  all_subs <- stringr::str_extract_all(
    mono_without_marker,
    single_sub_pattern
  )[[1]]

  # Sort substituents by position
  sub_string <- if (length(all_subs) > 0) {
    collapse_substituent_tokens(all_subs)
  } else {
    ""
  }

  c(mono = base_mono, sub = sub_string)
}
