# Access Individual Glycan Structures

Extract individual glycan structure graphs from a glycan structure
vector.

## Usage

``` r
get_structure_graphs(x, return_list = NULL)
```

## Arguments

- x:

  A glycan structure vector.

- return_list:

  If `TRUE`, always returns a list. If `FALSE` and `x` has a length of
  1, return the igraph object directly. If not provided (default),
  `FALSE` when `x` has a length of 1 and `TRUE` otherwise.

## Value

A list of igraph objects or an igraph object directly (see `return_list`
parameter).

## Examples

``` r
structures <- c(o_glycan_core_1(), n_glycan_core())
get_structure_graphs(structures)
#> [[1]]
#> IGRAPH 73f8b0f DN-- 2 1 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edge from 73f8b0f (vertex names):
#> [1] 2->1
#> 
#> [[2]]
#> IGRAPH 847df3f DN-- 5 4 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edges from 847df3f (vertex names):
#> [1] 3->1 3->2 4->3 5->4
#> 
get_structure_graphs(structures)
#> [[1]]
#> IGRAPH 73f8b0f DN-- 2 1 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edge from 73f8b0f (vertex names):
#> [1] 2->1
#> 
#> [[2]]
#> IGRAPH 847df3f DN-- 5 4 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edges from 847df3f (vertex names):
#> [1] 3->1 3->2 4->3 5->4
#> 
```
