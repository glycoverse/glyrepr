# Detect Floating Glycan Parts

Test whether each glycan structure contains one or more unresolved
floating parts.

## Usage

``` r
has_floating_parts(x)
```

## Arguments

- x:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  vector or a glycan `igraph`.

## Value

A logical vector with the same length and names as vector input, or a
logical scalar for graph input. Missing structures produce `NA`.

## Details

A floating part is a known glycan residue or substructure whose parent
residue on the main glycan tree is not fully localized. For example, a
bi-antennary N-glycan may contain one sialic acid while the available
evidence cannot determine which of its two terminal galactoses carries
that residue. The sialic acid can then be represented as a floating part
with both galactoses as candidate parents.

In the `glyrepr` IUPAC extension, floating parts appear in braces before
the main glycan:

- `{<floating>}<main>` means every feasible main-tree node is a
  candidate parent.

- `{<floating>|<parents>}<main>` restricts the candidates to the
  comma-separated main-tree node indices in `<parents>`.

Main-tree indices follow residue order within the main IUPAC-condensed
tree and are numbered from 1 independently of preceding floating parts.
A floating part may contain one residue or an entire subtree, and its
virtual attachment linkage may be fully known, partially known, or
unknown.

Internally, an unresolved structure is an annotated forest containing
one main tree and one disconnected tree per floating part. The virtual
linkage and candidate parents are stored as graph metadata rather than
as an edge. Each floating-part metadata entry also stores `nodes`, the
complete vector of vertex indices in that floating component, so
downstream graph operations do not need to rediscover its membership by
traversal.
[`structure_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
exposes attachment metadata in tabular form, and
[`structure_component_membership()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_component_membership.md)
exposes component membership.

An explicit singleton parent list fully localizes the attachment, as
does an omitted parent list when the main tree has only one node. Such a
part is normalized to an ordinary graph edge, so `has_floating_parts()`
returns `FALSE` for the normalized structure.

## See also

[`structure_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md),
[`has_floating_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/has_floating_substituents.md),
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)

## Examples

``` r
main <- paste0(
  "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
  "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
  "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
)
ambiguous <- as_glycan_structure(
  paste0("{Neu5Ac(a2-3)|1,4}", main)
)
glycans <- c(ambiguous = ambiguous, ordinary = as_glycan_structure(main))
has_floating_parts(glycans)
#> ambiguous  ordinary 
#>      TRUE     FALSE 
structure_floating_parts(ambiguous)
#> # A tibble: 1 × 6
#>   glycan_id part_id root_node nodes     linkage parents  
#>       <int>   <int>     <int> <list>    <chr>   <list>   
#> 1         1       1         1 <int [1]> a2-3    <int [2]>

localized <- as_glycan_structure(
  "{Neu5Ac(a2-3)|1}Gal(b1-4)GlcNAc(b1-"
)
has_floating_parts(localized)
#> [1] FALSE
```
