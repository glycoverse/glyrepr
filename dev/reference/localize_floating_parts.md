# Localize Floating Glycan Parts

`localize_floating_parts()` attaches selected floating parts to
caller-supplied parent nodes. The assignments are interpreted against
the floating-part rows and canonical node identifiers returned by
[`structure_floating_candidates()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_floating_candidates.md).

Each selected parent must belong to the floating part's declared
candidate domain. For an unrestricted `{<floating>}` part, that domain
contains every feasible node outside its own component, including nodes
in other floating components. Assignments must also be simultaneously
compatible with occupied and potential acceptor linkage positions,
acyclic, and ultimately connected to the main tree.

Selected virtual attachments become ordinary graph edges. Unassigned
floating parts remain floating. The resulting structures are
canonicalized, and candidate-parent indices for remaining parts are
remapped to the new canonical complete-sequence order. Attaching one
floating component to another merges their component metadata and can
iteratively resolve newly singleton domains.

Missing values, vector positions, and names in structure-vector input
are preserved. For graph input, `glycan_id` must be `1L`, and selected
edges are appended without canonicalizing or renumbering vertices.

## Usage

``` r
localize_floating_parts(x, assignments)
```

## Arguments

- x:

  A glycan structure vector or a glycan `igraph`.

- assignments:

  A data frame with integer columns `glycan_id`, `part_id`, and
  `parent_node`. Each `(glycan_id, part_id)` pair may occur at most
  once.

## Value

An object of the same representation as `x`. Structure-vector output has
the same length and names as `x`; graph output retains its vertex IDs
and order.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
)
assignments <- tibble::tibble(
  glycan_id = 1L,
  part_id = 1L,
  parent_node = 2L
)
localize_floating_parts(glycan, assignments)
#> <glycan_structure[1]>
#> [1] Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```
