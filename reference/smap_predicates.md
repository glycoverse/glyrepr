# Test Predicates on Glycan Structure Vectors

These functions test predicates on unique structures in a glycan
structure vector, taking advantage of hash-based deduplication to avoid
redundant computation. Similar to purrr predicate functions, but
optimized for glycan structure vectors.

## Usage

``` r
ssome(.x, .p, ...)

severy(.x, .p, ...)

snone(.x, .p, ...)
```

## Arguments

- .x:

  A glycan structure vector (glyrepr_structure).

- .p:

  A predicate function that takes an igraph object and returns a logical
  value. Can be a function, purrr-style lambda (`~ .x$attr`), or a
  character string naming a function.

- ...:

  Additional arguments passed to `.p`.

## Value

A single logical value.

## Details

These functions only evaluate `.p` once for each unique structure,
making them much more efficient than applying `.p` to each element
individually when there are duplicate structures.

**Return Values:**

- `ssome()`: Returns `TRUE` if at least one unique structure satisfies
  the predicate

- `severy()`: Returns `TRUE` if all unique structures satisfy the
  predicate

- `snone()`: Returns `TRUE` if no unique structures satisfy the
  predicate

## Examples

``` r
# Create a structure vector with duplicates
core1 <- o_glycan_core_1()
core2 <- n_glycan_core()
structures <- glycan_structure(core1, core2, core1)  # core1 appears twice

# Test if some structures have more than 5 vertices
ssome(structures, function(g) igraph::vcount(g) > 5)
#> [1] FALSE

# Test if all structures have at least 3 vertices
severy(structures, function(g) igraph::vcount(g) >= 3)
#> [1] FALSE

# Test if no structures have more than 20 vertices
snone(structures, function(g) igraph::vcount(g) > 20)
#> [1] TRUE

# Use purrr-style lambda functions
ssome(structures, ~ igraph::vcount(.x) > 5)
#> [1] FALSE
severy(structures, ~ igraph::vcount(.x) >= 3)
#> [1] FALSE
snone(structures, ~ igraph::vcount(.x) > 20)
#> [1] TRUE
```
