# Get Available Monosaacharides

This function returns a character vector of monosaccharide names of the
given type. See
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)
for monosaacharide types. Concrete furanose forms use an `f` after the
monosaccharide stem, such as `Galf` and `GlcfNAc`. Generic names do not
encode ring form.

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
#>   [1] "Hex"       "HexNAc"    "HexN"      "HexA"      "dHex"      "dHexNAc"  
#>   [7] "ddHex"     "Pen"       "NeuAc"     "NeuGc"     "gNeu"      "gKdn"     
#>  [13] "gPse"      "gLeg"      "gAci"      "g4eLeg"    "gBac"      "Hep"      
#>  [19] "gKdo"      "HepA"      "MurAc"     "MurGc"     "gMur"      "Glc"      
#>  [25] "Man"       "Gal"       "Gul"       "Alt"       "All"       "Tal"      
#>  [31] "Ido"       "GlcNAc"    "GalNAc"    "ManNAc"    "GulNAc"    "AltNAc"   
#>  [37] "AllNAc"    "TalNAc"    "IdoNAc"    "GlcN"      "ManN"      "GalN"     
#>  [43] "GulN"      "AltN"      "AllN"      "TalN"      "IdoN"      "GlcA"     
#>  [49] "ManA"      "GalA"      "GulA"      "AltA"      "AllA"      "TalA"     
#>  [55] "IdoA"      "Fuc"       "Qui"       "Rha"       "6dGul"     "6dAlt"    
#>  [61] "6dTal"     "QuiNAc"    "RhaNAc"    "6dAltNAc"  "6dTalNAc"  "FucNAc"   
#>  [67] "Oli"       "Tyv"       "Abe"       "Par"       "Dig"       "Col"      
#>  [73] "Ara"       "Lyx"       "Xyl"       "Rib"       "Neu5Ac"    "Neu5Gc"   
#>  [79] "Neu"       "Kdn"       "Pse"       "Leg"       "Aci"       "4eLeg"    
#>  [85] "Bac"       "LDmanHep"  "Kdo"       "Dha"       "DDmanHep"  "MurNAc"   
#>  [91] "MurNGc"    "Mur"       "Api"       "Fru"       "Tag"       "Sor"      
#>  [97] "Psi"       "Glcf"      "Manf"      "Galf"      "Gulf"      "Altf"     
#> [103] "Allf"      "Talf"      "Idof"      "GlcfNAc"   "GalfNAc"   "ManfNAc"  
#> [109] "GulfNAc"   "AltfNAc"   "AllfNAc"   "TalfNAc"   "IdofNAc"   "GlcfN"    
#> [115] "ManfN"     "GalfN"     "GulfN"     "AltfN"     "AllfN"     "TalfN"    
#> [121] "IdofN"     "GlcfA"     "ManfA"     "GalfA"     "GulfA"     "AltfA"    
#> [127] "AllfA"     "TalfA"     "IdofA"     "Fucf"      "Quif"      "Rhaf"     
#> [133] "6dGulf"    "6dAltf"    "6dTalf"    "QuifNAc"   "RhafNAc"   "6dAltfNAc"
#> [139] "6dTalfNAc" "FucfNAc"   "Olif"      "Tyvf"      "Abef"      "Parf"     
#> [145] "Digf"      "Colf"      "Araf"      "Lyxf"      "Xylf"      "Ribf"     
#> [151] "Neuf5Ac"   "Neuf5Gc"   "Neuf"      "Kdnf"      "Psef"      "Legf"     
#> [157] "Acif"      "4eLegf"    "Bacf"      "LDmanHepf" "Kdof"      "Dhaf"     
#> [163] "DDmanHepf" "MurfNAc"   "MurfNGc"   "Murf"      "Apif"      "Fruf"     
#> [169] "Tagf"      "Sorf"      "Psif"     
available_monosaccharides("concrete")
#>   [1] "Glc"       "Man"       "Gal"       "Gul"       "Alt"       "All"      
#>   [7] "Tal"       "Ido"       "GlcNAc"    "GalNAc"    "ManNAc"    "GulNAc"   
#>  [13] "AltNAc"    "AllNAc"    "TalNAc"    "IdoNAc"    "GlcN"      "ManN"     
#>  [19] "GalN"      "GulN"      "AltN"      "AllN"      "TalN"      "IdoN"     
#>  [25] "GlcA"      "ManA"      "GalA"      "GulA"      "AltA"      "AllA"     
#>  [31] "TalA"      "IdoA"      "Fuc"       "Qui"       "Rha"       "6dGul"    
#>  [37] "6dAlt"     "6dTal"     "QuiNAc"    "RhaNAc"    "6dAltNAc"  "6dTalNAc" 
#>  [43] "FucNAc"    "Oli"       "Tyv"       "Abe"       "Par"       "Dig"      
#>  [49] "Col"       "Ara"       "Lyx"       "Xyl"       "Rib"       "Neu5Ac"   
#>  [55] "Neu5Gc"    "Neu"       "Kdn"       "Pse"       "Leg"       "Aci"      
#>  [61] "4eLeg"     "Bac"       "LDmanHep"  "Kdo"       "Dha"       "DDmanHep" 
#>  [67] "MurNAc"    "MurNGc"    "Mur"       "Api"       "Fru"       "Tag"      
#>  [73] "Sor"       "Psi"       "Glcf"      "Manf"      "Galf"      "Gulf"     
#>  [79] "Altf"      "Allf"      "Talf"      "Idof"      "GlcfNAc"   "GalfNAc"  
#>  [85] "ManfNAc"   "GulfNAc"   "AltfNAc"   "AllfNAc"   "TalfNAc"   "IdofNAc"  
#>  [91] "GlcfN"     "ManfN"     "GalfN"     "GulfN"     "AltfN"     "AllfN"    
#>  [97] "TalfN"     "IdofN"     "GlcfA"     "ManfA"     "GalfA"     "GulfA"    
#> [103] "AltfA"     "AllfA"     "TalfA"     "IdofA"     "Fucf"      "Quif"     
#> [109] "Rhaf"      "6dGulf"    "6dAltf"    "6dTalf"    "QuifNAc"   "RhafNAc"  
#> [115] "6dAltfNAc" "6dTalfNAc" "FucfNAc"   "Olif"      "Tyvf"      "Abef"     
#> [121] "Parf"      "Digf"      "Colf"      "Araf"      "Lyxf"      "Xylf"     
#> [127] "Ribf"      "Neuf5Ac"   "Neuf5Gc"   "Neuf"      "Kdnf"      "Psef"     
#> [133] "Legf"      "Acif"      "4eLegf"    "Bacf"      "LDmanHepf" "Kdof"     
#> [139] "Dhaf"      "DDmanHepf" "MurfNAc"   "MurfNGc"   "Murf"      "Apif"     
#> [145] "Fruf"      "Tagf"      "Sorf"      "Psif"     
```
