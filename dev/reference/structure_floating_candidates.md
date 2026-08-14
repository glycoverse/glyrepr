# List Candidate Parents for Floating Metadata

`structure_floating_candidates()` expands floating-part and
floating-substituent metadata to one row per candidate parent. This
provides a uniform representation for explicitly restricted and
unrestricted floating metadata.

For an unrestricted floating part, every feasible node outside its own
component is returned; for an unrestricted floating substituent, every
feasible residue node in the complete structure is returned. These rows
use `scope = "all"`. For metadata written with an explicit `|<parents>`
suffix, only the declared parent nodes are returned and `scope` is
`"explicit"`.

Floating-part rows have a non-missing `part_id`, `root_node`, and
`linkage`. Floating-substituent rows instead have a non-missing
`substituent_id` and `substituent`. Exactly one of `part_id` and
`substituent_id` is non-missing in each row.

Node indices refer to `structure_nodes()$node_id` for the same glycan.
Missing structures and structures without floating metadata contribute
no rows. Duplicate structures are expanded to their original vector
positions. For graph input, node indices are current numeric vertex
positions and `glycan_id` is `1L`. If vector input is named, the result
also contains a `glycan_name` column.

## Usage

``` r
structure_floating_candidates(x)
```

## Arguments

- x:

  A glycan structure vector or one glycan `igraph`.

## Value

A tibble with columns `glycan_id`, `part_id`, `root_node`,
`parent_node`, `linkage`, `scope`, `substituent_id`, and `substituent`,
plus `glycan_name` when `x` is named.

## Examples

``` r
glycans <- as_glycan_structure(c(
  unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
  restricted = "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
  substituent = "{6S}Gal(a1-3)Gal(a1-"
))
structure_floating_candidates(glycans)
#> # A tibble: 6 × 9
#>   glycan_id glycan_name  part_id root_node parent_node linkage scope   
#>       <int> <chr>          <int>     <int>       <int> <chr>   <chr>   
#> 1         1 unrestricted       1         1           2 a2-3    all     
#> 2         1 unrestricted       1         1           3 a2-3    all     
#> 3         2 restricted         1         1           2 a2-6    explicit
#> 4         2 restricted         1         1           3 a2-6    explicit
#> 5         3 substituent       NA        NA           1 NA      all     
#> 6         3 substituent       NA        NA           2 NA      all     
#> # ℹ 2 more variables: substituent_id <int>, substituent <chr>
```
