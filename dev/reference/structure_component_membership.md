# Identify Main and Floating Structure Components

`structure_component_membership()` identifies whether each glycan node
belongs to the main tree or to a floating part. It provides a public,
normalized alternative to reading the private floating-part graph
attribute.

Main-tree nodes have `component_type = "main"` and a missing `part_id`.
Nodes in a floating subtree have `component_type = "floating"` and the
corresponding `part_id` from
[`structure_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md).
Floating substituents do not introduce vertices, so they do not create
additional component-membership rows.

Node indices refer to `structure_nodes()$node_id` for the same glycan.
Missing structures contribute no rows. Duplicate structures are expanded
to their original vector positions. For graph input, node indices are
current numeric vertex positions and `glycan_id` is `1L`. If vector
input is named, the result also contains a `glycan_name` column.

## Usage

``` r
structure_component_membership(x)
```

## Arguments

- x:

  A glycan structure vector or one glycan `igraph`.

## Value

A tibble with columns `glycan_id`, `node_id`, `component_type`, and
`part_id`, plus `glycan_name` when `x` is named.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
)
structure_component_membership(glycan)
#> # A tibble: 3 × 4
#>   glycan_id node_id component_type part_id
#>       <int>   <int> <chr>            <int>
#> 1         1       1 floating             1
#> 2         1       2 main                NA
#> 3         1       3 main                NA
```
