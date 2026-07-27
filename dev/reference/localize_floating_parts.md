# Localize Floating Glycan Parts

`localize_floating_parts()` attaches selected floating parts to
caller-supplied parent nodes. The assignments are interpreted against
the canonical node and part identifiers returned by
[`structure_floating_candidates()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_floating_candidates.md).

Each selected parent must belong to the floating part's declared
candidate domain. For an unrestricted `{<floating>}` part, that domain
contains every node in the original main tree. Assignments must also be
simultaneously compatible with occupied and potential acceptor linkage
positions.

Selected virtual attachments become ordinary graph edges. Unassigned
floating parts remain floating. The resulting structures are
canonicalized, and candidate-parent indices for remaining parts are
remapped to the new canonical main-tree order.

Missing values, vector positions, and names in `x` are preserved.

## Usage

``` r
localize_floating_parts(x, assignments)
```

## Arguments

- x:

  A glycan structure vector.

- assignments:

  A data frame with integer columns `glycan_id`, `part_id`, and
  `parent_node`. Each `(glycan_id, part_id)` pair may occur at most
  once.

## Value

A glycan structure vector with the same length and names as `x`.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
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
