#' Construct Structures from Residue and Edge Arrays
#'
#' A format-independent entry point for parsers. Node IDs are one-based positions
#' in `mono`; edges point from parent to child. Input order need not be canonical.
#' Each record contains `mono`, `sub`, `edges`, `linkage`, and `anomer`.
#' `edges` is an interleaved integer vector `c(parent, child, ...)`.
#' `sub` contains comma-separated substituent tokens (or empty strings).
#' Optional `alditol` defaults to `FALSE`.
#'
#' Optional `floating_parts` is a list of records containing `root`, `nodes`,
#' `linkage`, and `parents`. Each `nodes` vector must contain exactly the nodes
#' in that disconnected component. Optional `floating_substituents` is a list of
#' records containing `substituent` and `parents`. Empty `parents` means all
#' feasible nodes; all indices refer to the original input arrays. Singleton
#' candidates are resolved during canonicalization.
#'
#' @param x A list of structure records. A `NULL` element represents a missing
#'   structure. Names, missing positions, and duplicate positions are preserved.
#' @param on_failure Either `"error"` (default) or `"na"`. The latter warns and
#'   returns missing structures at invalid positions.
#' @returns A `glyrepr_structure` vector, with graph storage deduplicated by
#'   canonical IUPAC key.
#' @examples
#' structure_from_arrays(list(list(
#'   mono = c("Glc", "Gal"), sub = c("", ""),
#'   edges = c(1L, 2L), linkage = "b1-4", anomer = "?1"
#' )))
#' @export
structure_from_arrays <- function(x, on_failure = c("error", "na")) {
  checkmate::assert_list(x)
  on_failure <- match.arg(on_failure)
  if (!length(x)) {
    return(new_glycan_structure(stats::setNames(character(), names(x))))
  }
  outcomes <- lapply(seq_along(x), function(i) {
    if (is.null(x[[i]])) {
      return(NULL)
    }
    tryCatch(.structure_array_outcome(x[[i]]), error = identity)
  })
  failed <- vapply(outcomes, inherits, logical(1), "error")
  if (any(failed) && on_failure == "error") {
    i <- which(failed)[[1]]
    cli::cli_abort("Invalid structure at position {i}.", parent = outcomes[[i]])
  }
  present <- !vapply(outcomes, is.null, logical(1))
  assemble_recovered_structure_outcomes(
    outcomes[present],
    as.list(which(present)),
    length(x),
    names(x)
  )
}

.structure_array_outcome <- function(a) {
  a <- .validate_structure_arrays(a)
  process_glycan_structure_element(.compact_build_graph(a))
}

.validate_structure_arrays <- function(a) {
  checkmate::assert_list(a, names = "unique")
  required <- c("mono", "sub", "edges", "linkage", "anomer")
  if (!all(required %in% names(a))) {
    cli::cli_abort("Structure records require {.val {required}}.")
  }
  allowed <- c(required, "alditol", "floating_parts", "floating_substituents")
  if (any(!names(a) %in% allowed)) {
    cli::cli_abort(
      "Unknown structure record field: {.val {setdiff(names(a), allowed)}}."
    )
  }
  checkmate::assert_character(a$mono, min.len = 1, any.missing = FALSE)
  n <- length(a$mono)
  checkmate::assert_character(a$sub, len = n, any.missing = FALSE)
  .check_array_indices(a$edges, n)
  if (length(a$edges) %% 2L) {
    cli::cli_abort("Edges must contain parent-child pairs.")
  }
  checkmate::assert_character(
    a$linkage,
    len = length(a$edges) / 2,
    any.missing = FALSE
  )
  checkmate::assert_string(a$anomer)
  if (is.null(a$alditol)) {
    a$alditol <- FALSE
  }
  checkmate::assert_flag(a$alditol)
  for (field in c("floating_parts", "floating_substituents")) {
    if (is.null(a[[field]])) {
      next
    }
    checkmate::assert_list(a[[field]])
    a[[field]] <- lapply(a[[field]], function(p) {
      fields <- if (field == "floating_parts") {
        c("root", "nodes", "linkage", "parents")
      } else {
        c("substituent", "parents")
      }
      checkmate::assert_list(p, names = "unique")
      if (!setequal(names(p), fields)) {
        cli::cli_abort("Invalid {.field {field}} fields.")
      }
      .check_array_indices(p$parents, n, unique = TRUE)
      p$parents <- as.integer(p$parents)
      if (field == "floating_parts") {
        .check_array_indices(p$root, n)
        if (length(p$root) != 1L) {
          cli::cli_abort("Floating root must be one node ID.")
        }
        .check_array_indices(p$nodes, n, unique = TRUE)
        p$root <- as.integer(p$root)
        p$nodes <- as.integer(p$nodes)
        checkmate::assert_string(p$linkage)
      } else {
        checkmate::assert_string(p$substituent)
      }
      p
    })
  }
  a$edges <- as.integer(a$edges)
  a
}

.check_array_indices <- function(x, n, unique = FALSE) {
  checkmate::assert_integerish(x, lower = 1, upper = n, any.missing = FALSE)
  if (unique && anyDuplicated(x)) cli::cli_abort("Node IDs must be unique.")
}
