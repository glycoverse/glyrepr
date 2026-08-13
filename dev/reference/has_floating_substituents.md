# Detect Floating Substituents

Test whether each glycan structure contains one or more substituents
whose parent residue is unresolved.

## Usage

``` r
has_floating_substituents(x)
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

In the `glyrepr` IUPAC extension, floating substituents appear in braces
before the main glycan. `{6S}<main>` allows every feasible main-tree
residue to carry the substituent, `{6S|1,2}<main>` restricts it to
main-tree residues 1 and 2, and `{?S}<main>` also leaves the carbon
position unknown.

Internally, unresolved substituents are stored in the
`floating_substituents` graph attribute. It is a list with one entry per
substituent. Each entry contains a canonical `substituent` token and an
integer `parents` vector of candidate main-tree vertex indices. An empty
vector means every feasible main-tree vertex is a candidate.

A singleton candidate set is normalized into that vertex's ordinary
`sub` attribute. This also happens for an unrestricted floating
substituent when the main tree has only one residue.
[`structure_floating_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
exposes the metadata as a normalized table.

## See also

[`structure_floating_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md),
[`has_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/has_floating_parts.md),
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)

## Examples

``` r
glycans <- as_glycan_structure(c(
  floating = "{6S}Gal(a1-3)Gal(a1-",
  restricted = "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-",
  ordinary = "Gal6S(a1-"
))
has_floating_substituents(glycans)
#>   floating restricted   ordinary 
#>       TRUE       TRUE      FALSE 
```
