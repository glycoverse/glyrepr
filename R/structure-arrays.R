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
#' Records above the native backend's size guard use the graph reference path;
#' the guard is not an input-size limit. Only failing or unsupported records
#' use that path. No format-specific strings are generated and reparsed.
#' Failures have class `glyrepr_error_structure_failure` and fields `position`,
#' `input_name`, and `reason`. Recovery warnings have class
#' `glyrepr_warning_structure_failure` and fields `positions` and `reasons`.
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
  outcomes <- .structure_array_outcomes(x)
  failed <- vapply(outcomes, inherits, logical(1), "error")
  if (any(failed) && on_failure == "error") {
    i <- which(failed)[[1]]
    .abort_structure_failure(outcomes[[i]], i, names(x), rlang::current_env())
  }
  present <- !vapply(outcomes, is.null, logical(1))
  if (!any(present)) {
    return(new_glycan_structure(stats::setNames(
      rep(NA_character_, length(x)),
      names(x)
    )))
  }
  assemble_recovered_structure_outcomes(
    outcomes[present],
    as.list(which(present)),
    length(x),
    names(x)
  )
}

.structure_array_outcomes <- function(x) {
  checked <- lapply(x, function(a) {
    if (is.null(a)) {
      return(NULL)
    }
    tryCatch(.validate_structure_arrays(a), error = identity)
  })
  valid <- !vapply(checked, inherits, logical(1), "error")
  native <- .compact_arrays_native(
    checked[valid],
    base::order,
    .compact_bliss_labels,
    identical(Sys.getlocale("LC_COLLATE"), "C")
  )
  cache <- new.env(hash = TRUE, parent = emptyenv())
  valid_indices <- which(valid)
  for (j in seq_along(native)) {
    i <- valid_indices[[j]]
    a <- native[[j]]
    if (a$status == "missing") {
      next
    }
    checked[i] <- list(tryCatch(
      {
        if (a$status != "ok") {
          # Reference graph path handles unsupported sizes and supplies established
          # semantic errors without reparsing any format-specific string.
          process_glycan_structure_element(.compact_build_graph(checked[[i]]))
        } else {
          if (!exists(a$iupac, cache, inherits = FALSE)) {
            cache[[a$iupac]] <- .compact_build_graph(a)
          }
          list(graph = cache[[a$iupac]], iupac = a$iupac)
        }
      },
      error = identity
    ))
  }
  checked
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
  if (!all(is_known_mono(a$mono))) {
    cli::cli_abort("Unknown monosaccharide.")
  }
  if (!all(valid_substituent(a$sub))) {
    cli::cli_abort("Invalid substituent.")
  }
  if (!all(valid_linkages(a$linkage))) {
    cli::cli_abort("Invalid linkage.")
  }
  if (!valid_anomer(a$anomer)) {
    cli::cli_abort("Invalid anomer.")
  }
  for (p in a$floating_parts) {
    if (!valid_linkages(p$linkage)) cli::cli_abort("Invalid floating linkage.")
  }
  for (p in a$floating_substituents) {
    if (
      !grepl(
        substituent_token_pattern(anchored = TRUE),
        p$substituent,
        perl = TRUE
      ) ||
        !valid_substituent(p$substituent)
    ) {
      cli::cli_abort("Invalid floating substituent.")
    }
  }
  a$edges <- as.integer(a$edges)
  a
}

.check_array_indices <- function(x, n, unique = FALSE) {
  checkmate::assert_integerish(x, lower = 1, upper = n, any.missing = FALSE)
  if (unique && anyDuplicated(x)) cli::cli_abort("Node IDs must be unique.")
}

.abort_structure_failure <- function(cnd, position, input_names, call) {
  input_name <- if (is.null(input_names)) {
    NA_character_
  } else {
    input_names[[position]]
  }
  cli::cli_abort(
    "Invalid structure at position {position}.",
    class = "glyrepr_error_structure_failure",
    position = position,
    input_name = input_name,
    reason = normalize_structure_failure_reason(cnd),
    parent = cnd,
    call = call
  )
}
