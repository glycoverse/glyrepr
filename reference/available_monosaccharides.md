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
#>  [7] "ddHex"    "Pent"     "NeuAc"    "NeuGc"    "Glc"      "Man"     
#> [13] "Gal"      "Gul"      "Alt"      "All"      "Tal"      "Ido"     
#> [19] "GlcNAc"   "GalNAc"   "ManNAc"   "GulNAc"   "AltNAc"   "AllNAc"  
#> [25] "TalNAc"   "IdoNAc"   "GlcN"     "ManN"     "GalN"     "GulN"    
#> [31] "AltN"     "AllN"     "TalN"     "IdoN"     "GlcA"     "ManA"    
#> [37] "GalA"     "GulA"     "AltA"     "AllA"     "TalA"     "IdoA"    
#> [43] "Fuc"      "Qui"      "Rha"      "6dGul"    "6dAlt"    "6dTal"   
#> [49] "QuiNAc"   "RhaNAc"   "6dAltNAc" "6dTalNAc" "FucNAc"   "Oli"     
#> [55] "Tyv"      "Abe"      "Par"      "Dig"      "Col"      "Ara"     
#> [61] "Lyx"      "Xyl"      "Rib"      "Neu5Ac"   "Neu5Gc"   "Sia"     
#> [67] "Neu"      "Kdn"      "Pse"      "Leg"      "Aci"      "4eLeg"   
#> [73] "Bac"      "LDmanHep" "Kdo"      "Dha"      "DDmanHep" "MurNAc"  
#> [79] "MurNGc"   "Mur"      "Api"      "Fru"      "Tag"      "Sor"     
#> [85] "Psi"     
```
