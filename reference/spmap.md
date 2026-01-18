# Map Functions Over Glycan Structure Vectors and Multiple Arguments

These functions apply a function to each unique structure in a glycan
structure vector along with corresponding elements from multiple other
vectors, taking advantage of hash-based deduplication to avoid redundant
computation. Similar to purrr pmap functions, but optimized for glycan
structure vectors.

## Usage

``` r
spmap(.l, .f, ..., .parallel = FALSE)

spmap_vec(.l, .f, ..., .ptype = NULL, .parallel = FALSE)

spmap_lgl(.l, .f, ..., .parallel = FALSE)

spmap_int(.l, .f, ..., .parallel = FALSE)

spmap_dbl(.l, .f, ..., .parallel = FALSE)

spmap_chr(.l, .f, ..., .parallel = FALSE)

spmap_structure(.l, .f, ..., .parallel = FALSE)
```

## Arguments

- .l:

  A list where the first element is a glycan structure vector
  (glyrepr_structure) and the remaining elements are vectors of the same
  length or length 1 (will be recycled).

- .f:

  A function that takes an igraph object (from first element of `.l`)
  and values from other elements, returning a result. Can be a function,
  purrr-style lambda (`~ .x + .y + .z`), or a character string naming a
  function.

- ...:

  Additional arguments passed to `.f`.

- .parallel:

  Logical; whether to use parallel processing. If `FALSE` (default),
  parallel processing is disabled. Set to `TRUE` to enable parallel
  processing. See examples in
  [`smap`](https://glycoverse.github.io/glyrepr/reference/smap.md) for
  how to set up and use parallel processing.

- .ptype:

  A prototype for the return type (for `spmap_vec`).

## Value

- `spmap()`: A list

- `spmap_vec()`: An atomic vector of type specified by `.ptype`

- `spmap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type

- `spmap_structure()`: A new glyrepr_structure object

## Details

These functions only compute `.f` once for each unique combination of
structure and corresponding values from other vectors, then map the
results back to the original vector positions. This is much more
efficient than applying `.f` to each element combination individually
when there are duplicate combinations.

**Time Complexity Performance:**

Performance scales with unique combinations of all arguments rather than
total vector length. When argument vectors are highly redundant,
performance approaches O(unique_structures). Scaling factor shows time
increase when vector size increases 20x.

**Return Types:**

- `spmap()`: Returns a list with the same length as the input vectors

- `spmap_vec()`: Returns an atomic vector with the same length as the
  input vectors

- `spmap_lgl()`: Returns a logical vector

- `spmap_int()`: Returns an integer vector

- `spmap_dbl()`: Returns a double vector

- `spmap_chr()`: Returns a character vector

- `spmap_structure()`: Returns a new glycan structure vector (`.f` must
  return igraph objects)

## Examples

``` r
# Create structure vectors with duplicates
core1 <- o_glycan_core_1()
core2 <- n_glycan_core()
structures <- c(core1, core2, core1)  # core1 appears twice
weights <- c(1.0, 2.0, 1.0)  # corresponding weights
factors <- c(2, 3, 2)  # corresponding factors

# Map a function that uses structure, weight, and factor
spmap_dbl(list(structures, weights, factors), 
          function(g, w, f) igraph::vcount(g) * w * f)
#> [1]  4 30  4

# Use purrr-style lambda functions
spmap_dbl(list(structures, weights, factors), ~ igraph::vcount(..1) * ..2 * ..3)
#> [1]  4 30  4

# Test with recycling
spmap_dbl(list(structures, 2.0, 3), ~ igraph::vcount(..1) * ..2 * ..3)
#> [1] 12 30 12
```
