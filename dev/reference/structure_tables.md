# Convert Glycan Structures to Graph Tables

`structure_nodes()`, `structure_edges()`, `structure_floating_parts()`,
and `structure_floating_substituents()` convert a glycan structure
vector or one glycan `igraph` to normalized graph tables.
`structure_from_tibbles()` rebuilds a `glyrepr_structure` vector from
those tables, a vector of reducing-end anomers, and optional alditol
status.

The `glycan_id` column is the integer position of each glycan in the
input vector. Duplicate structures are expanded to their original vector
positions. Missing structures have no node or edge rows and are
reconstructed from missing values in `anomers`. If `x` is named, all
four tibbles also contain a `glycan_name` column.
`structure_from_tibbles()` uses `glycan_name` as output names when that
column is present.

For structure-vector input, `structure_nodes()$node_id` follows residue
order in the complete canonical IUPAC-condensed string. For graph input,
it follows the graph's current numeric vertex positions without
canonicalizing or renumbering them. A graph is represented with
`glycan_id = 1L` and no `glycan_name` column.

In `structure_floating_parts()`, `root_node` and every integer in the
`nodes` and `parents` list-columns refer to `structure_nodes()$node_id`
for the same glycan. `nodes` contains every node in the floating
component. An empty `parents` vector means all feasible main-tree nodes
are candidates. The `linkage` column describes the virtual attachment to
the main tree; this attachment is intentionally absent from
`structure_edges()`. During reconstruction, a row with exactly one
effective candidate parent is normalized to an ordinary edge and is
therefore absent from the resulting `structure_floating_parts()` table.

Parent indices written after `|` in an IUPAC-condensed floating part are
local to the main tree. `structure_floating_parts()` translates those
values to global `structure_nodes()$node_id` values, so they can differ
when floating nodes precede the main tree. `structure_from_tibbles()`
expects these global node IDs and translates them back during
serialization.

In `structure_floating_substituents()`, each row describes one
unresolved substituent. `substituent` is its canonical position-and-name
token, and the `parents` list-column contains candidate global node IDs.
An empty vector means all feasible main-tree nodes are candidates. A
singleton candidate is normalized into `structure_nodes()$sub`, so it
does not remain in the floating-substituent table.

## Usage

``` r
structure_nodes(x)

structure_edges(x)

structure_floating_parts(x)

structure_floating_substituents(x)

structure_from_tibbles(
  nodes,
  edges,
  anomers,
  floating_parts = NULL,
  floating_substituents = NULL,
  alditols = FALSE
)
```

## Arguments

- x:

  A glycan structure vector or one glycan `igraph`.

- nodes:

  A data frame with columns `glycan_id`, `node_id`, `mono`, and `sub`,
  and optionally `glycan_name`.

- edges:

  A data frame with columns `glycan_id`, `edge_id`, `from_node`,
  `to_node`, and `linkage`, and optionally `glycan_name`.

- anomers:

  A character vector of reducing-end anomers, one per glycan.

- floating_parts:

  A data frame returned by `structure_floating_parts()`, or `NULL` when
  no floating parts are present.

- floating_substituents:

  A data frame returned by `structure_floating_substituents()`, or
  `NULL` when no floating substituents are present.

- alditols:

  A logical vector indicating alditol status, either one value or one
  per glycan. Missing values are allowed only for missing glycans.
  Defaults to `FALSE` for backward compatibility.

## Value

- `structure_nodes()` returns a tibble with columns `glycan_id`,
  `node_id`, `mono`, and `sub`.

- `structure_edges()` returns a tibble with columns `glycan_id`,
  `edge_id`, `from_node`, `to_node`, and `linkage`.

- `structure_floating_parts()` returns a tibble with columns
  `glycan_id`, `part_id`, `root_node`, the list-column `nodes`,
  `linkage`, and the list-column `parents`.

- `structure_floating_substituents()` returns a tibble with columns
  `glycan_id`, `substituent_id`, `substituent`, and the list-column
  `parents`.

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

floating <- as_glycan_structure(
  "{8S|1,2}{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
)
floating_parts <- structure_floating_parts(floating)
floating_substituents <- structure_floating_substituents(floating)
structure_from_tibbles(
  structure_nodes(floating),
  structure_edges(floating),
  get_anomer(floating),
  floating_parts,
  floating_substituents
)
#> <glycan_structure[1]>
#> [1] {8S|1,2}{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```
