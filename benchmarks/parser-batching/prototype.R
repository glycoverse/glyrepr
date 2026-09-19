# Standalone R batching prototype. Source after pkgload::load_all().
batch_env <- new.env(parent = asNamespace("glyrepr"))
evalq(
  {
    tokens_to_arrays <- function(tokens) {
      n <- sum(!tokens %in% c("[", "]"))
      mono <- sub <- character(n)
      edges <- integer(2L * (n - 1L))
      linkage <- character(n - 1L)
      first <- .extract_substituent(tokens[[1]])
      mono[[1]] <- first[["mono"]]
      sub[[1]] <- first[["sub"]]
      node_stack <- rstackdeque::insert_top(rstackdeque::rstack(), 1L)
      current <- 1L
      next_id <- 1L
      for (token in tokens[-1L]) {
        if (token == "[") {
          node_stack <- rstackdeque::insert_top(node_stack, current)
        } else if (token == "]") {
          current <- rstackdeque::peek_top(node_stack)
          node_stack <- rstackdeque::without_top(node_stack)
        } else {
          parsed <- .parse_token(token)
          next_id <- next_id + 1L
          mono[[next_id]] <- parsed[["mono"]]
          sub[[next_id]] <- parsed[["sub"]]
          e <- next_id - 1L
          edges[2L * e - 1:0] <- c(current, next_id)
          linkage[[e]] <- parsed[["linkage"]]
          current <- next_id
        }
      }
      list(mono = mono, sub = sub, edges = edges, linkage = linkage)
    }
    graph_from_arrays <- function(a) {
      n <- length(a$mono)
      graph <- igraph::make_empty_graph(n, directed = TRUE)
      graph <- igraph::set_vertex_attr(
        graph,
        "name",
        value = as.character(seq_len(n))
      )
      graph <- igraph::set_vertex_attr(graph, "mono", value = a$mono)
      graph <- igraph::set_vertex_attr(graph, "sub", value = a$sub)
      if (length(a$edges)) {
        graph <- igraph::add_edges(graph, a$edges, linkage = a$linkage)
      } else {
        graph <- igraph::set_edge_attr(graph, "linkage", value = character())
      }
      graph
    }
    graph_incremental <- function(a) {
      graph <- igraph::make_empty_graph()
      for (v in seq_along(a$mono)) {
        graph <- igraph::add_vertices(
          graph,
          1,
          name = as.character(v),
          mono = a$mono[[v]],
          sub = a$sub[[v]]
        )
        if (v > 1L) {
          e <- v - 1L
          graph <- igraph::add_edges(
            graph,
            a$edges[2L * e - 1:0],
            linkage = a$linkage[[e]]
          )
        }
      }
      if (!length(a$edges)) {
        graph <- igraph::set_edge_attr(graph, "linkage", value = character())
      }
      graph
    }
    parse_tree_batch <- function(x) {
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
          alditol_result <- parse_alditol_iupac(x)
          x <- alditol_result$iupac
          alditol <- alditol_result$alditol
          x <- .infer_reducing_end_anomer(x)
          anomer <- .extract_anomer(x)
          x <- stringr::str_sub(x, 1, -stringr::str_length(anomer) - 3)

          tokens <- .tokenize_iupac(x)

          arrays <- tokens_to_arrays(tokens)
          graph <- graph_from_arrays(arrays)

          graph$anomer <- anomer
          graph$alditol <- alditol
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
  },
  batch_env
)
with_batch_parser <- function(code) {
  ns <- asNamespace("glyrepr")
  name <- ".parse_iupac_tree_single"
  original <- get(name, ns)
  locked <- bindingIsLocked(name, ns)
  set <- function(fun) {
    if (locked) {
      unlockBinding(name, ns)
    }
    assign(name, fun, ns)
    if (locked) lockBinding(name, ns)
  }
  on.exit(set(original), add = TRUE)
  set(batch_env$parse_tree_batch)
  force(code)
}
