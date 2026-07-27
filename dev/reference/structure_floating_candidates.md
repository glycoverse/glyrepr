# List Candidate Attachments for Floating Parts

`structure_floating_candidates()` expands floating-part attachment
metadata to one row per candidate parent. This provides a uniform
representation for both explicitly restricted floating parts and
unrestricted `{<floating>}` parts.

For an unrestricted floating part, every node in the original main tree
is returned as a candidate and `scope` is `"all"`. For a floating part
written with `{<floating>|<parents>}`, only the declared parent nodes
are returned and `scope` is `"explicit"`.

Node indices refer to `structure_nodes()$node_id` for the same glycan.
Missing structures and structures without floating parts contribute no
rows. Duplicate structures are expanded to their original vector
positions. If `x` is named, the result also contains a `glycan_name`
column.

## Usage

``` r
structure_floating_candidates(x)
```

## Arguments

- x:

  A glycan structure vector.

## Value

A tibble with columns `glycan_id`, `part_id`, `root_node`,
`parent_node`, `linkage`, and `scope`, plus `glycan_name` when `x` is
named.

## Examples

``` r
glycans <- as_glycan_structure(c(
  unrestricted = "{Neu5Ac(a2-3)}Gal(b1-3)GalNAc(a1-",
  restricted = "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
))
structure_floating_candidates(glycans)
#> # A tibble: 4 × 7
#>   glycan_id glycan_name  part_id root_node parent_node linkage scope   
#>       <int> <chr>          <int>     <int>       <int> <chr>   <chr>   
#> 1         1 unrestricted       1         1           2 a2-3    all     
#> 2         1 unrestricted       1         1           3 a2-3    all     
#> 3         2 restricted         1         1           2 a2-6    explicit
#> 4         2 restricted         1         1           3 a2-6    explicit
```
