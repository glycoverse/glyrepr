# List Potential Virtual Edges for Floating Parts

`structure_candidate_edges()` represents every potential floating-part
attachment as an explicit virtual edge. `from_node` is a candidate
parent in the main tree and `to_node` is the root of the floating part.

The rows correspond one-to-one with the floating-part rows from
[`structure_floating_candidates()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_floating_candidates.md).
Floating substituents do not create virtual graph edges. For
unrestricted `{<floating>}` parts, every original main-tree node is
returned and `scope` is `"all"`. For explicitly restricted parts, only
the declared parent nodes are returned and `scope` is `"explicit"`.

Node indices refer to `structure_nodes()$node_id` for the same glycan.
Missing structures and structures without floating parts contribute no
rows. Duplicate structures are expanded to their original vector
positions. For graph input, node indices are current numeric vertex
positions and `glycan_id` is `1L`. If vector input is named, the result
also contains a `glycan_name` column.

## Usage

``` r
structure_candidate_edges(x)
```

## Arguments

- x:

  A glycan structure vector or one glycan `igraph`.

## Value

A tibble with columns `glycan_id`, `part_id`, `from_node`, `to_node`,
`linkage`, and `scope`, plus `glycan_name` when `x` is named.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
)
structure_candidate_edges(glycan)
#> # A tibble: 2 × 6
#>   glycan_id part_id from_node to_node linkage scope   
#>       <int>   <int>     <int>   <int> <chr>   <chr>   
#> 1         1       1         2       1 a2-6    explicit
#> 2         1       1         3       1 a2-6    explicit
```
