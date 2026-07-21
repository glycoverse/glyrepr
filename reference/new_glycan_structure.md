# Construct a Glycan Structure Vector from Trusted Data

Assemble a glycan structure vector from IUPAC-condensed values and a
graph lookup table without graph validation, canonicalization, IUPAC
generation, vector-level compatibility checks, or graph deduplication.

## Usage

``` r
new_glycan_structure(iupac = character(), graphs = list())
```

## Arguments

- iupac:

  A character vector of canonical IUPAC-condensed strings. Missing
  values are allowed, and names are preserved exactly.

- graphs:

  A named list of valid, canonical, mutually compatible `igraph` glycan
  graphs keyed by IUPAC-condensed strings.

## Value

A `glyrepr_structure` vector.

## Details

`graphs` must be a named list containing one graph for every distinct,
non-missing value in `iupac`. Additional named graphs are allowed so
that vctrs prototypes can retain their graph lookup tables. Graph names
must be unique and non-missing. This function checks these inexpensive
representation invariants but trusts that each graph matches its key.

## Low-level API warning

These functions are low-level, developer-facing APIs. Calling them
directly is usually not a good idea unless you understand and can
guarantee all glycan graph and `glyrepr_structure` invariants. Prefer
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.md)
for ordinary construction. Incorrect use of these functions can create
invalid structure vectors that fail in later operations.

## Name-preserving manual construction

The five low-level functions can reproduce strict graph-based
construction while preserving the names of the input graph list:

    input_names <- names(graphs)
    graphs <- unname(graphs)

    graphs <- purrr::map(graphs, validate_glycan_graph)
    graphs <- purrr::map(graphs, canonicalize_glycan_graph)
    validate_glycan_graph_vector(graphs)

    iupacs <- purrr::map_chr(graphs, graph_to_iupac)
    names(iupacs) <- input_names

    unique <- !duplicated(unname(iupacs))
    unique_graphs <- graphs[unique]
    names(unique_graphs) <- unname(iupacs[unique])

    new_glycan_structure(iupacs, unique_graphs)

Unlike `as_glycan_structure(graphs, on_failure = "na")`, this strict
pipeline stops at the first invalid graph.

## See also

Other low-level glycan structure functions:
[`canonicalize_glycan_graph()`](https://glycoverse.github.io/glyrepr/reference/canonicalize_glycan_graph.md),
[`graph_to_iupac()`](https://glycoverse.github.io/glyrepr/reference/graph_to_iupac.md),
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/reference/validate_glycan_graph.md),
[`validate_glycan_graph_vector()`](https://glycoverse.github.io/glyrepr/reference/validate_glycan_graph_vector.md)
