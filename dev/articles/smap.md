# Efficient Glycan Manipulation with smap

## Overview

This vignette introduces the `smap` family of functions. These functions
are useful when you want to apply custom `igraph`-based operations to
glycan structure vectors.

This guide assumes you are comfortable with R programming and have some
familiarity with graph concepts. If you are just getting started, read
the “Getting Started with glyrepr” vignette first.

``` r
library(glyrepr)
```

## Unique Structure Optimization

Before using `smap`, it helps to understand why these functions exist.

### The Problem

Working with glycan structures means working with graphs, and graph
operations are computationally expensive. When you are analyzing
thousands of glycans from a large-scale study, this becomes a real
bottleneck.

### The Solution

`glyrepr` implements an optimization called **unique structure
storage**. Instead of storing thousands of identical graphs, it stores
only the unique ones and keeps track of which original positions they
belong to.

Let’s see this in action:

``` r
# Our test data: some common glycan structures
iupacs <- c(
  "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # N-glycan core
  "Gal(b1-3)GalNAc(a1-",                                    # O-glycan core 1
  "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",                     # O-glycan core 2
  "Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-",          # Branched mannose
  "GlcNAc6Ac(b1-4)Glc3Me(a1-"                              # With decorations
)

struc <- as_glycan_structure(iupacs)

# Now create a realistic dataset with lots of repetition.
large_struc <- rep(struc, 1000)  # 5,000 total structures
large_struc
#> <glycan_structure[5000]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [2] Gal(b1-3)GalNAc(a1-
#> [3] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [4] Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-
#> [5] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> [6] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [7] Gal(b1-3)GalNAc(a1-
#> [8] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [9] Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-
#> [10] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> ... (4990 more not shown)
#> # Unique structures: 5
```

Notice that the object reports only 5 unique structures. The vector has
5,000 elements, but only 5 unique graphs are stored internally.

We can verify that directly:

``` r
# Only 5 unique graphs are stored internally
length(attr(large_struc, "structures"))
#> [1] 0

# But we have 5,000 total elements
length(large_struc)
#> [1] 5000
```

### Memory Savings

``` r
library(lobstr)
obj_sizes(struc, large_struc)
#> * 15.18 kB
#> * 40.69 kB
```

The memory difference can be substantial. For repeated structures, the
optimized representation can be much smaller than storing every graph
independently.

## The `smap` Family

There is one important consequence of this internal representation:
regular [`lapply()`](https://rdrr.io/r/base/lapply.html) or
[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
functions do not operate directly on a glycan structure vector as if it
were a list of graphs.

``` r
# This will not work and will raise an error.
tryCatch(
  purrr::map_int(large_struc, ~ igraph::vcount(.x)),
  error = function(e) cat("Error:", rlang::cnd_message(e))
)
#> Error: ℹ In index: 1.
#> Caused by error in `ensure_igraph()`:
#> ! Must provide a graph object (provided wrong object type).
```

**Why does this fail?** Because `purrr` functions don’t understand the
internal structure optimization of `glycan_structure` objects.

### Structure-Aware Mapping

The `smap` functions are structure-aware alternatives to `purrr` mapping
functions. They understand the unique structure optimization and work
directly with the underlying graph objects.

``` r
vertex_counts <- smap_int(large_struc, ~ igraph::vcount(.x))
vertex_counts[1:10]
#>  [1] 5 2 3 5 2 5 2 3 5 2
```

The “s” stands for “structure”: these functions operate on the
underlying `igraph` objects that represent glycan structures.

## The `smap` Toolkit

The `smap` family provides glycan-aware equivalents for virtually all
`purrr` functions:

| purrr                                                         | smap                                                                                | purrr                                                           | smap                                                                         |
|---------------------------------------------------------------|-------------------------------------------------------------------------------------|-----------------------------------------------------------------|------------------------------------------------------------------------------|
| [`map()`](https://purrr.tidyverse.org/reference/map.html)     | [`smap()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)              | [`map2()`](https://purrr.tidyverse.org/reference/map2.html)     | [`smap2()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)     |
| [`map_lgl()`](https://purrr.tidyverse.org/reference/map.html) | [`smap_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)          | [`map2_lgl()`](https://purrr.tidyverse.org/reference/map2.html) | [`smap2_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md) |
| [`map_int()`](https://purrr.tidyverse.org/reference/map.html) | [`smap_int()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)          | [`map2_int()`](https://purrr.tidyverse.org/reference/map2.html) | [`smap2_int()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md) |
| [`map_dbl()`](https://purrr.tidyverse.org/reference/map.html) | [`smap_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)          | [`map2_dbl()`](https://purrr.tidyverse.org/reference/map2.html) | [`smap2_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md) |
| [`map_chr()`](https://purrr.tidyverse.org/reference/map.html) | [`smap_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)          | [`map2_chr()`](https://purrr.tidyverse.org/reference/map2.html) | [`smap2_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md) |
| [`some()`](https://purrr.tidyverse.org/reference/every.html)  | [`ssome()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md)  | [`pmap()`](https://purrr.tidyverse.org/reference/pmap.html)     | [`spmap()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)     |
| [`every()`](https://purrr.tidyverse.org/reference/every.html) | [`severy()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md) | `pmap_*()`                                                      | `spmap_*()`                                                                  |
| [`none()`](https://purrr.tidyverse.org/reference/every.html)  | [`snone()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md)  | [`imap()`](https://purrr.tidyverse.org/reference/imap.html)     | [`simap()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)     |
|                                                               |                                                                                     | `imap_*()`                                                      | `simap_*()`                                                                  |

As a simple rule, replace `map` with `smap`, `pmap` with `spmap`, and
`imap` with `simap`. The function signatures are designed to feel
familiar if you already use `purrr`.

### Basic Examples

**Count vertices in each structure:**

``` r
vertex_counts <- smap_int(large_struc, igraph::vcount)
summary(vertex_counts)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     2.0     2.0     3.0     3.4     5.0     5.0
```

**Find structures with more than 4 vertices:**

``` r
has_many_vertices <- smap_lgl(large_struc, ~ igraph::vcount(.x) > 4)
sum(has_many_vertices)
#> [1] 2000
```

**Get the degree sequence of each structure:**

``` r
degree_sequences <- smap(large_struc, ~ igraph::degree(.x))
degree_sequences[1:3]
#> [[1]]
#> 1 2 3 4 5 
#> 1 1 3 2 1 
#> 
#> [[2]]
#> 1 2 
#> 1 1 
#> 
#> [[3]]
#> 1 2 3 
#> 1 1 2
```

**Check if any structure has isolated vertices:**

``` r
ssome(large_struc, ~ any(igraph::degree(.x) == 0))
#> [1] FALSE
```

**Verify all structures are connected:**

``` r
severy(large_struc, ~ igraph::is_connected(.x))
#> [1] TRUE
```

### Beyond Basic `smap()`

Quick examples of the extended family:

``` r
# smap2: Apply function with additional parameters
thresholds <- c(3, 4, 5)
large_enough <- smap2_lgl(struc[1:3], thresholds, function(g, threshold) {
  igraph::vcount(g) >= threshold
})
large_enough
#> [1]  TRUE FALSE FALSE
```

``` r
# simap: Include position information
indexed_report <- simap_chr(large_struc[1:3], function(g, i) {
  paste0("#", i, ": ", igraph::vcount(g), " vertices")
})
indexed_report
#> [1] "#1: 5 vertices" "#2: 2 vertices" "#3: 3 vertices"
```

**Performance note:** `simap` functions do not benefit from the unique
structure optimization. Since each element has a different index, the
combination of `(structure, index)` is always unique, breaking the
deduplication that makes other `smap` functions fast. Use `simap` only
when you truly need position information.

## Performance

The main performance benefit of `smap` functions comes from automatic
deduplication:

``` r
# Create a large dataset with high redundancy
huge_struc <- rep(struc, 5000)  # 25,000 structures, only 5 unique

cat("Dataset size:", length(huge_struc), "structures\n")
#> Dataset size: 25000 structures
cat("Unique structures:", length(attr(huge_struc, "structures")), "\n")
#> Unique structures: 0
cat("Redundancy factor:", length(huge_struc) / length(attr(huge_struc, "structures")), "x\n")
#> Redundancy factor: Inf x

library(tictoc)

# Optimized approach: smap only processes 5 unique structures
tic("smap_int (optimized)")
vertex_counts_optimized <- smap_int(huge_struc, igraph::vcount)
toc()
#> smap_int (optimized): 0.003 sec elapsed

# Naive approach: extract all graphs and process each one
tic("Naive approach (all graphs)")
all_graphs <- get_structure_graphs(huge_struc)  # Extracts all 25,000 graphs
vertex_counts_naive <- purrr::map_int(all_graphs, igraph::vcount)
toc()
#> Naive approach (all graphs): 0.201 sec elapsed

# Verify results are equivalent (though data types may differ)
all.equal(vertex_counts_optimized, vertex_counts_naive)
#> [1] TRUE
```

The higher the redundancy, the larger the performance gain. In real
glycoproteomics datasets with repeated structures, this optimization can
provide about 10x speedups.

## Additional Patterns

### Working with Complex Functions

The function you pass to `smap` must accept an `igraph` object as its
first argument. You can use purrr-style lambda notation:

``` r
# Calculate clustering coefficient for each structure
clustering_coeffs <- smap_dbl(large_struc, ~ igraph::transitivity(.x, type = "global"))
summary(clustering_coeffs)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     NAs 
#>       0       0       0       0       0       0    2000
```

### Combining Multiple Metrics

``` r
# Create a compact structure summary.
structure_analysis <- smap(large_struc, function(g) {
  list(
    vertices = igraph::vcount(g),
    edges = igraph::ecount(g),
    diameter = ifelse(igraph::is_connected(g), igraph::diameter(g), NA),
    clustering = igraph::transitivity(g, type = "global")
  )
})

# Convert to a more usable format
analysis_df <- do.call(rbind, lapply(structure_analysis, data.frame))
head(analysis_df)
#>   vertices edges diameter clustering
#> 1        5     4        3          0
#> 2        2     1        1        NaN
#> 3        3     2        1          0
#> 4        5     4        2          0
#> 5        2     1        1        NaN
#> 6        5     4        3          0
```

### Memory-Efficient Filtering

``` r
# Find only structures with exactly 5 vertices
has_five_vertices <- smap_lgl(large_struc, ~ igraph::vcount(.x) == 5)
five_vertex_structures <- large_struc[has_five_vertices]

cat("Found", sum(has_five_vertices), "structures with exactly 5 vertices\n")
#> Found 2000 structures with exactly 5 vertices
```

## When to Use `smap` Functions

**Use `smap` functions when:**

- You need to apply `igraph`-based functions to glycan structures.
- You want better performance with datasets containing repeated
  structures.
- You are building custom glycan analysis pipelines.

**Use regular R functions when:**

- You are working with compositions.
- You are operating on string representations.

**Special note on `simap`:**

While `simap` functions are convenient for position-aware operations,
they do not provide performance benefits over regular `imap` functions.
The inclusion of index information breaks the unique structure
optimization, making each `(structure, index)` pair unique even when
structures are identical.

## Example: Custom Motif Detection

Here’s how you might build a custom glycan analysis pipeline using
`smap` functions:

``` r
# Custom motif detector
detect_branching <- function(g) {
  degrees <- igraph::degree(g)
  any(degrees >= 3)
}

# Apply to a large dataset using unique structure optimization.
has_branching <- smap_lgl(large_struc, detect_branching)
cat("Structures with branching:", sum(has_branching), "out of", length(large_struc), "\n")
#> Structures with branching: 2000 out of 5000

# Use smap2 to check structures against complexity thresholds
complexity_thresholds <- rep(c(3, 4, 5, 2, 4), 1000)  # Thresholds for each structure
meets_threshold <- smap2_lgl(large_struc, complexity_thresholds, function(g, threshold) {
  igraph::vcount(g) >= threshold
})
cat("Structures meeting complexity threshold:", sum(meets_threshold), "out of", length(large_struc), "\n")
#> Structures meeting complexity threshold: 2000 out of 5000
```

## Summary

The `smap` family provides structure-aware mapping functions for glycan
structure vectors. It lets you write custom graph-based analyses while
preserving the unique structure optimization used by `glyrepr`.

Key takeaways:

- Unique structure optimization stores repeated structures efficiently.
- `smap` functions are `purrr`-like tools that understand glycan
  structure vectors.
- Performance gains are strongest when datasets contain repeated
  structures.
- Use `smap` for structures, and use regular R or `purrr` functions for
  other data types.

## Session Information

``` r
sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] tictoc_1.2.1        lobstr_1.2.1        glyrepr_0.10.1.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0    dplyr_1.2.1       compiler_4.6.0    tidyselect_1.2.1 
#>  [5] stringr_1.6.0     jquerylib_0.1.4   systemfonts_1.3.2 textshaping_1.0.5
#>  [9] yaml_2.3.12       fastmap_1.2.0     R6_2.6.1          generics_0.1.4   
#> [13] igraph_2.3.0      knitr_1.51        backports_1.5.1   checkmate_2.3.4  
#> [17] tibble_3.3.1      rstackdeque_1.1.1 desc_1.4.3        bslib_0.10.0     
#> [21] pillar_1.11.1     rlang_1.2.0       cachem_1.1.0      stringi_1.8.7    
#> [25] xfun_0.57         fs_2.1.0          sass_0.4.10       cli_3.6.6        
#> [29] pkgdown_2.2.0     magrittr_2.0.5    digest_0.6.39     lifecycle_1.0.5  
#> [33] prettyunits_1.2.0 vctrs_0.7.3       evaluate_1.0.5    glue_1.8.1       
#> [37] ragg_1.5.2        rmarkdown_2.31    purrr_1.2.2       tools_4.6.0      
#> [41] pkgconfig_2.0.3   htmltools_0.5.9
```
