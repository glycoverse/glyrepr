# Convert Glycan Structures to Graph Tables

`structure_nodes()` and `structure_edges()` convert a glycan structure
vector to node and edge tibbles. `structure_from_tibbles()` rebuilds a
`glyrepr_structure` vector from those tibbles and a vector of
reducing-end anomers.

The `glycan_id` column is the integer position of each glycan in the
input vector. Duplicate structures are expanded to their original vector
positions. Missing structures have no node or edge rows and are
reconstructed from missing values in `anomers`. If `x` is named, the
node and edge tibbles also contain a `glycan_name` column.
`structure_from_tibbles()` uses `glycan_name` as output names when that
column is present.

## Usage

``` r
structure_nodes(x)

structure_edges(x)

structure_from_tibbles(nodes, edges, anomers)
```

## Arguments

- x:

  A glycan structure vector.

- nodes:

  A data frame with columns `glycan_id`, `node_id`, `mono`, and `sub`,
  and optionally `glycan_name`.

- edges:

  A data frame with columns `glycan_id`, `edge_id`, `from_node`,
  `to_node`, and `linkage`, and optionally `glycan_name`.

- anomers:

  A character vector of reducing-end anomers, one per glycan.

## Value

- `structure_nodes()` returns a tibble with columns `glycan_id`,
  `node_id`, `mono`, and `sub`.

- `structure_edges()` returns a tibble with columns `glycan_id`,
  `edge_id`, `from_node`, `to_node`, and `linkage`.

- `structure_from_tibbles()` returns a `glyrepr_structure` vector.

## Examples

``` r
glycans <- c(o_glycan_core_1(), o_glycan_core_1())
nodes <- structure_nodes(glycans)
edges <- structure_edges(glycans)
structure_from_tibbles(nodes, edges, get_anomer(glycans))
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```
