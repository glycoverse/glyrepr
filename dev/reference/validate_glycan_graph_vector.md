# Validate a List of Glycan Graphs

Check the container used to store individually valid glycan graphs in
one glycan structure vector. Generic and concrete residues may coexist
within a graph and across graphs.

## Usage

``` r
validate_glycan_graph_vector(graphs, label = NULL)
```

## Arguments

- graphs:

  A list of individually valid `igraph` glycan graphs.

- label:

  An optional label retained for backward compatibility.

## Value

`NULL`, invisibly.

## Details

This function assumes that every element has already passed
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md).
It does not repeat scalar graph validation.

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
source attributes differ. Graph and array construction share native
processing. Unsupported inputs use the R path, which also supplies the
established failure diagnostics.
[`structure_from_tibbles()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
converts table rows directly to array records.

The individual low-level functions remain available when an intermediate
graph is needed.
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md)
validates one graph;
[`canonicalize_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/canonicalize_glycan_graph.md)
assumes a valid graph;
[`graph_to_iupac()`](https://glycoverse.github.io/glyrepr/dev/reference/graph_to_iupac.md)
assumes a valid, canonical graph. `validate_glycan_graph_vector()`
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

Other low-level glycan structure functions:
[`canonicalize_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/canonicalize_glycan_graph.md),
[`graph_to_iupac()`](https://glycoverse.github.io/glyrepr/dev/reference/graph_to_iupac.md),
[`new_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/new_glycan_structure.md),
[`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md)
