# Get Available Monosaacharides

This function returns a character vector of monosaccharide names of the
given type. See
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/reference/get_mono_type.md)
for monosaacharide types.

## Usage

``` r
available_monosaccharides(mono_type = "all")
```

## Arguments

- mono_type:

  A character string specifying the type of monosaccharides. Can be
  "all", "generic", or "concrete". Default is "all".

## Value

A character vector of monosaccharide names.

## Examples

``` r
available_monosaccharides()
#>  [1] "Hex"      "HexNAc"   "HexN"     "HexA"     "dHex"     "dHexNAc" 
#>  [7] "ddHex"    "Pen"      "NeuAc"    "NeuGc"    "gNeu"     "gKdn"    
#> [13] "gPse"     "gLeg"     "gAci"     "g4eLeg"   "gBac"     "Hep"     
#> [19] "gKdo"     "HepA"     "MurAc"    "MurGc"    "gMur"     "Glc"     
#> [25] "Man"      "Gal"      "Gul"      "Alt"      "All"      "Tal"     
#> [31] "Ido"      "GlcNAc"   "GalNAc"   "ManNAc"   "GulNAc"   "AltNAc"  
#> [37] "AllNAc"   "TalNAc"   "IdoNAc"   "GlcN"     "ManN"     "GalN"    
#> [43] "GulN"     "AltN"     "AllN"     "TalN"     "IdoN"     "GlcA"    
#> [49] "ManA"     "GalA"     "GulA"     "AltA"     "AllA"     "TalA"    
#> [55] "IdoA"     "Fuc"      "Qui"      "Rha"      "6dGul"    "6dAlt"   
#> [61] "6dTal"    "QuiNAc"   "RhaNAc"   "6dAltNAc" "6dTalNAc" "FucNAc"  
#> [67] "Oli"      "Tyv"      "Abe"      "Par"      "Dig"      "Col"     
#> [73] "Ara"      "Lyx"      "Xyl"      "Rib"      "Neu5Ac"   "Neu5Gc"  
#> [79] "Neu"      "Kdn"      "Pse"      "Leg"      "Aci"      "4eLeg"   
#> [85] "Bac"      "LDmanHep" "Kdo"      "Dha"      "DDmanHep" "MurNAc"  
#> [91] "MurNGc"   "Mur"      "Api"      "Fru"      "Tag"      "Sor"     
#> [97] "Psi"     
```
