# Convert to Glycan Structure Vector

Convert an object to a glycan structure vector.

## Usage

``` r
as_glycan_structure(x, on_failure = c("error", "na"))
```

## Arguments

- x:

  An object to convert to a glycan structure vector. Can be an igraph
  object, a list of igraph objects, a character vector of
  IUPAC-condensed strings, or an existing glyrepr_structure object.

- on_failure:

  The failure policy for element-local parsing, validation, and
  canonicalization errors. `"error"` preserves the default strict
  behavior. `"na"` replaces failed elements with `NA` and emits one
  warning that reports their positions and failure reasons. Existing
  missing elements remain missing without a warning. Vector-level
  incompatibilities still produce an error.

## Value

A glyrepr_structure object.

## Details

Character input assumes the natural absolute configuration for
unprefixed monosaccharides. Less common configurations use a leading
`D-` or `L-`, such as `D-Fuc`, `L-Gul`, and `D-Fucf`. Alditols use `-ol`
on the main reducing-end residue, for example `Gal(b1-4)GlcNAc-ol(a1-`.
The reducing-end anomer annotation remains part of the canonical
representation.

Character input supports floating-part blocks before the main
IUPAC-condensed structure. `{Neu5Ac(a2-3)}<main>` allows every feasible
node outside its own component as a candidate parent, while an explicit
`|<parents>` suffix restricts that domain. Parent indices follow residue
order in the complete supplied sequence: residues in floating blocks are
counted left to right before the main glycan, and substituent blocks add
no indices. A floating part may target another floating component or the
main tree, but cannot target itself. Indices are remapped to canonical
complete sequence order in the result. The suffix is a `glyrepr`
extension to curly-brace IUPAC notation. A singleton candidate set is
accepted as input but fully localizes the attachment, so
`{Neu5Ac(a2-3)|2}Gal(b1-4)GlcNAc(b1-` canonicalizes to the ordinary
structure `Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-`.

Floating substituents use the same leading-brace and candidate-parent
syntax. For example, `{6S}<main>` leaves the sulfated residue
unrestricted across all residue nodes, `{6S|1,2}<main>` restricts it to
complete-sequence nodes 1 and 2, and `{?S}<main>` also leaves the carbon
position unknown. A singleton candidate is normalized into the selected
residue's ordinary `sub` attribute.

## Examples

``` r
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union

# Convert a single igraph
graph <- make_graph(~ 1-+2)
V(graph)$mono <- c("GlcNAc", "GlcNAc")
V(graph)$sub <- ""
E(graph)$linkage <- "b1-4"
graph$anomer <- "a1"
as_glycan_structure(graph)
#> <glycan_structure[1]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> # Unique structures: 1

# Convert a list of igraphs
o_glycan_vec <- o_glycan_core_1()
o_glycan_graph <- get_structure_graphs(o_glycan_vec)
as_glycan_structure(list(graph, o_glycan_graph))
#> <glycan_structure[2]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> [2] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 2

# Convert a character vector of IUPAC-condensed strings
as_glycan_structure(c("GlcNAc(b1-4)GlcNAc(b1-", "Man(a1-2)GlcNAc(b1-"))
#> <glycan_structure[2]>
#> [1] GlcNAc(b1-4)GlcNAc(b1-
#> [2] Man(a1-2)GlcNAc(b1-
#> # Unique structures: 2
as_glycan_structure(c("D-Fuc(a1-", "L-Gul(b1-", "D-Fucf(a1-"))
#> <glycan_structure[3]>
#> [1] D-Fuc(a1-
#> [2] L-Gul(b1-
#> [3] D-Fucf(a1-
#> # Unique structures: 3
as_glycan_structure("Gal(b1-4)GlcNAc-ol(a1-")
#> <glycan_structure[1]>
#> [1] Gal(b1-4)GlcNAc-ol(a1-
#> # Unique structures: 1

# Parse a floating residue with two candidate parents
floating_iupac <- paste0(
  "{Neu5Ac(a2-3)|2,5}",
  "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
  "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
  "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
)
as_glycan_structure(floating_iupac)
#> <glycan_structure[1]>
#> [1] {Neu5Ac(a2-3)|2,5}Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1

# Preserve valid elements while replacing an invalid element with NA
as_glycan_structure(
  c(valid = "Glc(?1-", invalid = "not-a-structure"),
  on_failure = "na"
)
#> Warning: 1 structure failed validation and was replaced with `NA`.
#> ✖ Position 2 (`invalid`): Could not parse IUPAC-condensed string:
#>   "not-a-structure" ℹ Invalid characters or format in IUPAC-condensed string
#> <glycan_structure[2]>
#> [1] valid    Glc(?1-
#> [2] invalid  NA
#> # Unique structures: 1
```
