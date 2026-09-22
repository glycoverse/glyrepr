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

## Choosing a construction pipeline

Prefer
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)
for ordinary construction from IUPAC strings or graphs. Parsers that
already produce residue and edge arrays should use
[`structure_from_arrays()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_from_arrays.md)
to validate and canonicalize records through the shared compact backend,
then build deduplicated graph storage.

For existing graph parsers,
[`canonicalize_glycan_graphs()`](https://glycoverse.github.io/glyrepr/dev/reference/canonicalize_glycan_graphs.md)
validates a batch and returns aligned canonical graphs, IUPAC keys, and
per-element status. It preserves input names and arbitrary graph,
vertex, and edge attributes. Equal structures retain separate graphs
until the caller explicitly deduplicates them. This is useful when
source attributes differ.

The individual low-level functions remain available when an intermediate
graph is needed.
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md)
validates one graph;
[`canonicalize_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/canonicalize_glycan_graph.md)
assumes a valid graph;
[`graph_to_iupac()`](https://glycoverse.github.io/glyrepr/dev/reference/graph_to_iupac.md)
assumes a valid, canonical graph.
[`validate_glycan_graph_vector()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph_vector.md)
checks the graph-list container, not each graph's semantics.
[`new_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/new_glycan_structure.md)
trusts that graphs are valid and canonical and match their keys; it does
not validate, canonicalize, or deduplicate them. Incorrect use can
create vectors that fail in later operations.

## Floating graph schemas

A floating structure is one weakly disconnected graph with exactly one
main outward tree and one outward tree per floating part. Its
`floating_parts` graph attribute is a list of entries with integer
`root`, integer `nodes`, character `linkage`, and integer `parents`
fields. `nodes` contains every vertex in that floating component.
`parents = integer()` means all feasible nodes outside that component,
including nodes in other floating components. Legacy input graphs may
omit `nodes`; canonical output graphs always contain it. During
canonicalization, a part with exactly one effective candidate parent is
converted to an ordinary graph edge and its floating metadata is
removed. Otherwise, the virtual attachment linkage is not a graph edge.
See
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
for the complete contract.

A graph may also have a `floating_substituents` attribute. It is a list
of entries with character `substituent` and integer `parents` fields. An
empty parent vector means all feasible residue nodes in the complete
structure. A singleton candidate is moved into the selected vertex's
`sub` attribute during canonicalization. All parent indices refer to the
complete graph's vertex positions, not positions within a component.

## Name-preserving graph construction

Use the batch interface to obtain canonical graphs and keys together,
then explicitly deduplicate graph storage when constructing a structure
vector:

    batch <- canonicalize_glycan_graphs(graphs)
    keep <- batch$status == "ok" & !duplicated(batch$iupac)
    unique_graphs <- batch$graphs[keep]
    names(unique_graphs) <- unname(batch$iupac[keep])

    new_glycan_structure(batch$iupac, unique_graphs)

Names and missing positions are preserved; use `NULL` for missing input
graphs. With the default `on_failure = "error"`, an invalid graph raises
an error identifying its original position. Set `on_failure = "na"` to
warn and retain successful entries, with invalid positions represented
by missing keys and `NULL` graphs. The same assembly code works in both
modes. Deduplication keeps the first graph for each key, including its
source attributes. Keep `batch$graphs` instead when per-input provenance
is needed. Use `validate = FALSE` only when every non-missing graph has
already passed
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md).

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
#> IGRAPH 9d71828 DN-- 3 2 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from 9d71828 (vertex names):
#> [1] 3->2 2->1
#> 
#> [[2]]
#> IGRAPH cee27b0 DN-- 3 2 -- 
#> + attr: anomer (g/c), alditol (g/l), name (v/c), mono (v/c), sub (v/c),
#> | linkage (e/c)
#> + edges from cee27b0 (vertex names):
#> [1] 3->2 3->1
#> 
```
