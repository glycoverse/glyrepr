# Map Functions Over Two Glycan Structure Vectors

These functions apply a function to each unique structure combination in
two glycan structure vectors, taking advantage of hash-based
deduplication to avoid redundant computation. Similar to purrr map2
functions, but optimized for glycan structure vectors.

## Usage

``` r
smap2(.x, .y, .f, ..., .parallel = FALSE)

smap2_vec(.x, .y, .f, ..., .ptype = NULL, .parallel = FALSE)

smap2_lgl(.x, .y, .f, ..., .parallel = FALSE)

smap2_int(.x, .y, .f, ..., .parallel = FALSE)

smap2_dbl(.x, .y, .f, ..., .parallel = FALSE)

smap2_chr(.x, .y, .f, ..., .parallel = FALSE)

smap2_structure(.x, .y, .f, ..., .parallel = FALSE)
```

## Arguments

- .x:

  A glycan structure vector (glyrepr_structure).

- .y:

  A vector of the same length as `.x`, or length 1 (will be recycled).

- .f:

  A function that takes an igraph object (from `.x`) and a value (from
  `.y`) and returns a result. Can be a function, purrr-style lambda
  (`~ .x + .y`), or a character string naming a function.

- ...:

  Additional arguments passed to `.f`.

- .parallel:

  Logical; whether to use parallel processing. If `FALSE` (default),
  parallel processing is disabled. Set to `TRUE` to enable parallel
  processing. See examples in
  [`smap`](https://glycoverse.github.io/glyrepr/reference/smap.md) for
  how to set up and use parallel processing.

- .ptype:

  A prototype for the return type (for `smap2_vec`).

## Value

- `smap2()`: A list

- `smap2_vec()`: An atomic vector of type specified by `.ptype`

- `smap2_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type

- `smap2_structure()`: A new glyrepr_structure object

## Details

These functions only compute `.f` once for each unique combination of
structure and corresponding `.y` value, then map the results back to the
original vector positions. This is much more efficient than applying
`.f` to each element pair individually when there are duplicate
structure-value combinations.

**Return Types:**

- `smap2()`: Returns a list with the same length as `.x`

- `smap2_vec()`: Returns an atomic vector with the same length as `.x`

- `smap2_lgl()`: Returns a logical vector

- `smap2_int()`: Returns an integer vector

- `smap2_dbl()`: Returns a double vector

- `smap2_chr()`: Returns a character vector

- `smap2_structure()`: Returns a new glycan structure vector (`.f` must
  return igraph objects)

## Examples

``` r
# Create structure vectors with duplicates
core1 <- o_glycan_core_1()
core2 <- n_glycan_core()
structures <- glycan_structure(core1, core2, core1)  # core1 appears twice
weights <- c(1.0, 2.0, 1.0)  # corresponding weights

# Map a function that uses both structure and weight
smap2_dbl(structures, weights, function(g, w) igraph::vcount(g) * w)
#> [1]  2 10  2

# Use purrr-style lambda functions  
smap2_dbl(structures, weights, ~ igraph::vcount(.x) * .y)
#> [1]  2 10  2

# Test with recycling (single weight for all structures)
smap2_dbl(structures, 2.5, ~ igraph::vcount(.x) * .y)
#> [1]  5.0 12.5  5.0

# Map a function that modifies structure based on second argument
# This example adds a graph attribute instead of modifying topology
add_weight_attr <- function(g, weight) {
  igraph::set_graph_attr(g, "weight", weight)
}
weights_to_add <- c(1.5, 2.5, 1.5)
smap2_structure(structures, weights_to_add, add_weight_attr)
#> <glycan_structure[3]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [3] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 2
```
