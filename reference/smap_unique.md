# Apply Function to Unique Structures Only

Apply a function only to the unique structures in a glycan structure
vector, returning results in the same order as the unique structures
appear. This is useful when you need to perform expensive computations
but only care about unique results.

## Usage

``` r
smap_unique(.x, .f, ..., .parallel = FALSE)
```

## Arguments

- .x:

  A glycan structure vector (glyrepr_structure).

- .f:

  A function that takes an igraph object and returns a result. Can be a
  function, purrr-style lambda (`~ .x$attr`), or a character string
  naming a function.

- ...:

  Additional arguments passed to `.f`.

- .parallel:

  Logical; whether to use parallel processing. If `FALSE` (default),
  parallel processing is disabled. Set to `TRUE` to enable parallel
  processing. See examples in
  [`smap`](https://glycoverse.github.io/glyrepr/reference/smap.md) for
  how to set up and use parallel processing.

## Value

A list with results for each unique structure, named by their hash
codes.

## Examples

``` r
# Create a structure vector with duplicates
core1 <- o_glycan_core_1()
structures <- c(core1, core1, core1)  # same structure 3 times

# Only compute once for the unique structure
unique_results <- smap_unique(structures, igraph::vcount)
length(unique_results)  # 1, not 3
#> [1] 1

# Use purrr-style lambda
unique_results2 <- smap_unique(structures, ~ igraph::vcount(.x))
length(unique_results2)  # 1, not 3
#> [1] 1
```
