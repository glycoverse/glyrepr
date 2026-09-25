# Construct Structures from Residue and Edge Arrays

A format-independent entry point for parsers. Node IDs are one-based
positions in `mono`; edges point from parent to child. Input order need
not be canonical. Each record contains `mono`, `sub`, `edges`,
`linkage`, and `anomer`. `edges` is an interleaved integer vector
`c(parent, child, ...)`. `sub` contains comma-separated substituent
tokens (or empty strings). Optional `alditol` defaults to `FALSE`.

## Usage

``` r
structure_from_arrays(x, on_failure = c("error", "na"), progress = FALSE)
```

## Arguments

- x:

  A list of structure records. A `NULL` element represents a missing
  structure. Names, missing positions, and duplicate positions are
  preserved.

- on_failure:

  Either `"error"` (default) or `"na"`. The latter warns and returns
  missing structures at invalid positions.

- progress:

  `FALSE` (default) for silent operation, `TRUE` for a cli progress bar,
  or a function with arguments `stage`, `current`, and `total`. A
  callback owns its display; no cli bar is created. Counts are local to
  each stage, not an overall percentage. Character parsing counts unique
  non-missing inputs; graph construction may count deduplicated
  canonical structures. Updates are throttled, with stage boundaries
  always reported. Empty or already-converted inputs may report no
  stages. Callback errors abort the operation, including with
  `on_failure = "na"`. The cli bar is closed when the call exits,
  including on errors or interrupts, and follows cli display settings
  (short or non-interactive calls may show no visible bar).

## Value

A `glyrepr_structure` vector, with graph storage deduplicated by
canonical IUPAC key.

## Details

Optional `floating_parts` is a list of records containing `root`,
`nodes`, `linkage`, and `parents`. Each `nodes` vector must contain
exactly the nodes in that disconnected component. Optional
`floating_substituents` is a list of records containing `substituent`
and `parents`. Empty `parents` means all feasible nodes; all indices
refer to the original input arrays. Singleton candidates are resolved
during canonicalization.

Records above the native backend's size guard use the graph reference
path; the guard is not an input-size limit. Only failing or unsupported
records use that path. No format-specific strings are generated and
reparsed. Failures have class `glyrepr_error_structure_failure` and
fields `position`, `input_name`, and `reason`. Recovery warnings have
class `glyrepr_warning_structure_failure` and fields `positions` and
`reasons`.

## Examples

``` r
structure_from_arrays(list(list(
  mono = c("Glc", "Gal"), sub = c("", ""),
  edges = c(1L, 2L), linkage = "b1-4", anomer = "?1"
)))
#> <glycan_structure[1]>
#> [1] Gal(b1-4)Glc(?1-
#> # Unique structures: 1
```
