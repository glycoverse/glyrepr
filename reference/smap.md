# Map Functions Over Glycan Structure Vectors

These functions apply a function to each unique structure in a glycan
structure vector, taking advantage of hash-based deduplication to avoid
redundant computation. Similar to purrr mapping functions, but optimized
for glycan structure vectors.

## Usage

``` r
smap(.x, .f, ..., .parallel = FALSE)

smap_vec(.x, .f, ..., .ptype = NULL, .parallel = FALSE)

smap_lgl(.x, .f, ..., .parallel = FALSE)

smap_int(.x, .f, ..., .parallel = FALSE)

smap_dbl(.x, .f, ..., .parallel = FALSE)

smap_chr(.x, .f, ..., .parallel = FALSE)

smap_structure(.x, .f, ..., .parallel = FALSE)
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
  processing.

- .ptype:

  A prototype for the return type (for `smap_vec`).

## Value

- `smap()`: A list

- `smap_vec()`: An atomic vector of type specified by `.ptype`

- `smap_lgl/int/dbl/chr()`: Atomic vectors of the corresponding type

- `smap_structure()`: A new glyrepr_structure object

## Details

These functions only compute `.f` once for each unique structure, then
map the results back to the original vector positions. This is much more
efficient than applying `.f` to each element individually when there are
duplicate structures.

**Return Types:**

- `smap()`: Returns a list with the same length as `.x`

- `smap_vec()`: Returns an atomic vector with the same length as `.x`

- `smap_lgl()`: Returns a logical vector

- `smap_int()`: Returns an integer vector

- `smap_dbl()`: Returns a double vector

- `smap_chr()`: Returns a character vector

- `smap_structure()`: Returns a new glycan structure vector (`.f` must
  return igraph objects)

## Examples

``` r
# Create a structure vector with duplicates
core1 <- o_glycan_core_1()
core2 <- n_glycan_core()
structures <- glycan_structure(core1, core2, core1)  # core1 appears twice

# Map a function that counts vertices - only computed twice, not three times
smap_int(structures, igraph::vcount)
#> [1] 2 5 2

# Map a function that returns logical
smap_lgl(structures, function(g) igraph::vcount(g) > 5)
#> [1] FALSE FALSE FALSE

# Use purrr-style lambda functions  
smap_int(structures, ~ igraph::vcount(.x))
#> [1] 2 5 2
smap_lgl(structures, ~ igraph::vcount(.x) > 5)
#> [1] FALSE FALSE FALSE

# Map a function that modifies structure (must return igraph)
add_vertex_names <- function(g) {
  if (!("name" %in% igraph::vertex_attr_names(g))) {
    igraph::set_vertex_attr(g, "name", value = paste0("v", seq_len(igraph::vcount(g))))
  } else {
    g
  }
}
smap_structure(structures, add_vertex_names)
#> <glycan_structure[3]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [3] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 2
```
