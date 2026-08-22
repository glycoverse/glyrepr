# Access Individual Glycan Structures

Extract individual glycan structure graphs from a glycan structure
vector. A structure with floating parts is returned as one annotated,
weakly disconnected `igraph`: its main tree and floating components
share the graph, and the `floating_parts` graph attribute records each
component's node indices, virtual attachment, and candidate parents. See
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
for the metadata schema. A structure with floating substituents carries
a `floating_substituents` graph attribute containing their tokens and
candidate parent indices.

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
  `FALSE` when `x` has a length of 1 and `TRUE` otherwise, including for
  an empty vector.

## Value

A list of igraph objects or an igraph object directly (see `return_list`
parameter).

## Examples

``` r
structures <- c(o_glycan_core_1(), n_glycan_core())
get_structure_graphs(structures)
#> [[1]]
#> IGRAPH 9a04905 DN-- 2 1 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edge from 9a04905 (vertex names):
#> [1] 2->1
#> 
#> [[2]]
#> IGRAPH 941c0b5 DN-- 5 4 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from 941c0b5 (vertex names):
#> [1] 3->1 3->2 4->3 5->4
#> 
get_structure_graphs(structures)
#> [[1]]
#> IGRAPH 9a04905 DN-- 2 1 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edge from 9a04905 (vertex names):
#> [1] 2->1
#> 
#> [[2]]
#> IGRAPH 941c0b5 DN-- 5 4 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from 941c0b5 (vertex names):
#> [1] 3->1 3->2 4->3 5->4
#> 
```
