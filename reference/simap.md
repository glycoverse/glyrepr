# Map Functions Over Glycan Structure Vectors with Indices

These functions apply a function to each unique structure in a glycan
structure vector along with their corresponding indices, taking
advantage of hash-based deduplication to avoid redundant computation.
Similar to purrr imap functions, but optimized for glycan structure
vectors.

## Usage

``` r
simap(.x, .f, ...)

simap_vec(.x, .f, ..., .ptype = NULL)

simap_lgl(.x, .f, ...)

simap_int(.x, .f, ...)

simap_dbl(.x, .f, ...)

simap_chr(.x, .f, ...)

simap_structure(.x, .f, ...)
```

## Arguments

- .x:

  A glycan structure vector (glyrepr_structure).

- .f:

  A function that takes an igraph object (from `.x`) and an index/name,
  returning a result. Can be a function, purrr-style lambda
  (`~ paste(.x, .y)`), or a character string naming a function.

- ...:

  Additional arguments passed to `.f`.

- .ptype:

  A prototype for the return type (for `simap_vec`).

## Value

- `simap()`: A list

- `simap_vec()`: An atomic vector of type specified by `.ptype`

- `simap_lgl()`: Returns a logical vector

- `simap_int()`: Returns an integer vector

- `simap_dbl()`: Returns a double vector

- `simap_chr()`: Returns a character vector

- `simap_structure()`: A new glyrepr_structure object

## Details

These functions only compute `.f` once for each unique combination of
structure and corresponding index/name, then map the results back to the
original vector positions. This is much more efficient than applying
`.f` to each element individually when there are duplicate structures.

**IMPORTANT PERFORMANCE NOTE:** Due to the inclusion of position
indices, `simap` functions have **O(total_structures)** time complexity
because each position creates a unique combination, even with identical
structures.

**Alternative:** Consider
[`smap()`](https://glycoverse.github.io/glyrepr/reference/smap.md)
functions if position information is not required.

The index passed to `.f` is the position in the original vector
(1-based). If the vector has names, the names are passed instead of
indices.

**Return Types:**

- `simap()`: Returns a list with the same length as `.x`

- `simap_vec()`: Returns an atomic vector with the same length as `.x`

- `simap_lgl()`: Returns a logical vector

- `simap_int()`: Returns an integer vector

- `simap_dbl()`: Returns a double vector

- `simap_chr()`: Returns a character vector

- `simap_structure()`: Returns a new glycan structure vector (`.f` must
  return igraph objects)

## Examples

``` r
# Create structure vectors with duplicates
core1 <- o_glycan_core_1()
core2 <- n_glycan_core()
structures <- c(core1, core2, core1)  # core1 appears twice

# Map a function that uses both structure and index
simap_chr(structures, function(g, i) paste0("Structure_", i, "_vcount_", igraph::vcount(g)))
#> [1] "Structure_1_vcount_2" "Structure_2_vcount_5" "Structure_3_vcount_2"

# Use purrr-style lambda functions
simap_chr(structures, ~ paste0("Pos", .y, "_vertices", igraph::vcount(.x)))
#> [1] "Pos1_vertices2" "Pos2_vertices5" "Pos3_vertices2"
```
