# Enumerate Floating Graph Localizations

`enumerate_floating_graph_localizations()` generates every
conflict-free, fully localized graph permitted by the candidate-parent
domains in `graph`. It is the graph-level counterpart to
[`enumerate_floating_localizations()`](https://glycoverse.github.io/glyrepr/dev/reference/enumerate_floating_localizations.md).

The returned graphs are validated but are not canonicalized. Adding the
selected attachment edges does not reorder vertices, so every vertex
keeps the same integer ID, name, and attributes as in `graph`.
Assignment `parent_node` values therefore refer directly to vertex IDs
in both the input and localized graphs.

Floating parts are localized as ordinary edges, while floating
substituents are localized into the selected parent vertex's `sub`
attribute. Every conflict-free acyclic assignment that connects all
components to the main tree is retained, even when multiple assignments
would produce the same canonical IUPAC-condensed structure. A graph
without floating metadata produces one identity row with an empty
assignment table.

`max_variants` is a conservative safeguard. It limits the raw Cartesian
product before conflict filtering, so the function may ask for a higher
bound even when fewer variants would ultimately remain.

## Usage

``` r
enumerate_floating_graph_localizations(graph, max_variants = 256)
```

## Arguments

- graph:

  A valid glycan `igraph`, optionally containing unresolved floating
  parts or substituents.

- max_variants:

  A positive integer giving the maximum raw candidate combinations
  allowed.

## Value

A tibble with columns:

- `variant_id`: the sequential localization identifier.

- `graph`: a list-column of fully localized `igraph` objects whose
  vertex IDs are identical to those in `graph`.

- `assignments`: a list-column of tibbles with `glycan_id`, `part_id`,
  `parent_node`, and `substituent_id`. Exactly one of `part_id` and
  `substituent_id` is non-missing in each row. `glycan_id` is always
  `1L`.

## Low-level API warning

These functions are low-level, developer-facing APIs. Calling them
directly is usually not a good idea unless you understand and can
guarantee all glycan graph and `glyrepr_structure` invariants. Prefer
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)
for ordinary construction. Incorrect use of these functions can create
invalid structure vectors that fail in later operations.

## Floating graph schemas

A floating structure is one weakly disconnected graph with exactly one
main outward tree and one outward tree per floating part. Its
`floating_parts` graph attribute is a list of entries with integer
`root`, integer `nodes`, character `linkage`, and integer `parents`
fields. `nodes` contains every vertex in that floating component.
`parents = integer()` means all feasible main-tree nodes. Legacy input
graphs may omit `nodes`; canonical output graphs always contain it.
During canonicalization, a part with exactly one effective candidate
parent is converted to an ordinary graph edge and its floating metadata
is removed. Otherwise, the virtual attachment linkage is not a graph
edge. See
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
for the complete contract.

A graph may also have a `floating_substituents` attribute. It is a list
of entries with character `substituent` and integer `parents` fields. An
empty parent vector means all feasible main-tree nodes. A singleton
candidate is moved into the selected vertex's `sub` attribute during
canonicalization.

## Name-preserving manual construction

The five low-level functions can reproduce strict graph-based
construction while preserving the names of the input graph list:

    input_names <- names(graphs)
    graphs <- unname(graphs)

    graphs <- purrr::map(graphs, validate_glycan_graph)
    graphs <- purrr::map(graphs, canonicalize_glycan_graph)
    validate_glycan_graph_vector(graphs)

    iupacs <- purrr::map_chr(graphs, graph_to_iupac)
    names(iupacs) <- input_names

    unique <- !duplicated(unname(iupacs))
    unique_graphs <- graphs[unique]
    names(unique_graphs) <- unname(iupacs[unique])

    new_glycan_structure(iupacs, unique_graphs)

Unlike `as_glycan_structure(graphs, on_failure = "na")`, this strict
pipeline stops at the first invalid graph.

## See also

[`enumerate_floating_localizations()`](https://glycoverse.github.io/glyrepr/dev/reference/enumerate_floating_localizations.md)

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
)
graph <- get_structure_graphs(glycan, return_list = FALSE)
localizations <- enumerate_floating_graph_localizations(graph)
localizations$graph
#> [[1]]
#> IGRAPH 887d0d7 DN-- 3 2 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from 887d0d7 (vertex names):
#> [1] 3->2 2->1
#> 
#> [[2]]
#> IGRAPH a86191e DN-- 3 2 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from a86191e (vertex names):
#> [1] 3->2 3->1
#> 
```
