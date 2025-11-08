# Example Glycan Structures

Create example glycan structures for testing and demonstration. Includes
**N-glycan core** and **O-glycan core 1** and **core 2**.

## Usage

``` r
n_glycan_core(linkage = TRUE, mono_type = "concrete")

o_glycan_core_1(linkage = TRUE, mono_type = "concrete")

o_glycan_core_2(linkage = TRUE, mono_type = "concrete")
```

## Arguments

- linkage:

  A logical indicating whether to include linkages (e.g. "b1-4").
  Default is `TRUE`.

- mono_type:

  A character string specifying the type of monosaccharides. Can be
  "generic" (Hex, HexNAc, dHex, NeuAc, etc.) or "concrete" (Man, Gal,
  HexNAc, Fuc, etc.). Default is "concrete".

## Value

A glycan structure (igraph) object.

## N-Glycan Core

**N-Glycans** are branched oligosaccharides that are bound, most
commonly, via GlcNAc to an Asn residue of the protein backbone. A common
motif of all N-glycans is the **chitobiose core**, composed of three
mannose and two GlcNAc moieties, which is commonly attached to the
protein backbone via GlcNAc. The mannose residue is branched and
connected via a1,3- and a1,6-glycosidic linkages to the two other
mannose building blocks.

        Man
      a1-6 \   b1-4      b1-4      b1-
            Man -- GlcNAc -- GlcNAc -
      a1-3 /
        Man

## O-Glycan Core

**O-Glycans** are highly abundant in extracellular proteins. Generally,
O-glycans are extended following four major core structures: **core 1**,
**core 2**, core 3, and core 4. The first two are by far the most common
core structures in O-glycosylation and are found throughout the body.

**core 1**:

              a1-
        GalNAc -
       / b1-3
    Gal

**core 2**:

    GlcNAc
          \ b1-6 a1-
           GalNAc -
          / b1-3
       Gal

## Examples

``` r
print(n_glycan_core(), verbose = TRUE)
#> <glycan_structure[1]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1
print(o_glycan_core_1(), verbose = TRUE)
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```
