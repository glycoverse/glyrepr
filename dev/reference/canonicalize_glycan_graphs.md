# Validate and Canonicalize a Batch of Glycan Graphs

Return canonical graphs and their IUPAC keys together. For ordinary
trees, ordering and key generation share one traversal. Graph, vertex,
and edge attributes are retained with their corresponding objects.
Unlike structure vector construction, this function does not deduplicate
graphs: equal structures can retain different source attributes.

## Usage

``` r
canonicalize_glycan_graphs(
  graphs,
  validate = TRUE,
  on_failure = c("error", "na")
)
```

## Arguments

- graphs:

  A list of glycan `igraph` objects. `NULL` elements represent missing
  structures. List names and positions are preserved.

- validate:

  Whether to validate each graph before canonicalization. Set to `FALSE`
  only for graphs already validated with
  [`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md).
  This skips semantic validation, not canonicalization. Array records
  should use
  [`structure_from_arrays()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_from_arrays.md).

- on_failure:

  Either `"error"` (default) or `"na"`. With `"na"`, invalid graphs
  produce a warning and missing output, with details in `reason`.

## Value

A list containing aligned, named vectors `iupac`, `status`, and
`reason`, and an aligned named list `graphs`. Status is `"ok"`,
`"missing"`, or `"invalid"`. Missing and invalid entries have a `NULL`
graph and `NA` key. Reasons are `NA` except for invalid entries.

## Details

Strict failures have class `glyrepr_error_structure_failure` with
`position`, `input_name`, and `reason` fields. Recovery warnings have
class `glyrepr_warning_structure_failure` with `positions` and `reasons`
fields.

## Examples

``` r
graphs <- as.list(n_glycan_core())
result <- canonicalize_glycan_graphs(graphs)
result$iupac
#> [1] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
new_glycan_structure(result$iupac, stats::setNames(result$graphs, result$iupac))
#> <glycan_structure[1]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1
```
