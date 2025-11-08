# Convert Glycan Structure to IUPAC-like Sequence

Convert a glycan structure to a sequence representation in the form of
mono(linkage)mono, with branches represented by square brackets \[\].
The backbone is chosen as the longest path, and for branches, linkages
are ordered lexicographically with smaller linkages on the backbone.

## Usage

``` r
structure_to_iupac(glycan)
```

## Arguments

- glycan:

  A glyrepr_structure vector.

## Value

A character vector representing the IUPAC sequences.

## Sequence Format

The sequence follows the format mono(linkage)mono, where:

- mono: monosaccharide name with optional substituents (e.g., Glc,
  GlcNAc, Glc3Me)

- linkage: glycosidic linkage (e.g., b1-4, a1-3)

- Branches are enclosed in square brackets \[\]

- Substituents are appended directly to monosaccharide names (e.g.,
  Glc3Me for Glc with 3Me substituent)

## Backbone Selection

The backbone is selected as the longest path in the tree. For branches,
the same rule applies recursively.

## Linkage Comparison

Linkages are compared lexicographically:

1.  First by anomeric configuration: ? \> b \> a

2.  Then by first position: ? \> numbers (numerically)

3.  Finally by second position: ? \> numbers (numerically)

Smaller linkages are placed on the backbone, larger ones in branches.

## Examples

``` r
# Simple linear structure
structure_to_iupac(o_glycan_core_1())
#> [1] "Gal(b1-3)GalNAc(a1-"

# Branched structure  
structure_to_iupac(n_glycan_core())
#> [1] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

# Structure with substituents
graph <- igraph::make_graph(~ 1-+2)
igraph::V(graph)$mono <- c("Glc", "GlcNAc")
igraph::V(graph)$sub <- c("3Me", "6Ac")
igraph::E(graph)$linkage <- "b1-4"
graph$anomer <- "a1"
glycan <- glycan_structure(graph)
structure_to_iupac(glycan)  # Returns "GlcNAc6Ac(b1-4)Glc3Me(a1-"
#> [1] "GlcNAc6Ac(b1-4)Glc3Me(a1-"

# Vectorized structures
structs <- glycan_structure(o_glycan_core_1(), n_glycan_core())
structure_to_iupac(structs)
#> [1] "Gal(b1-3)GalNAc(a1-"                                
#> [2] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
```
