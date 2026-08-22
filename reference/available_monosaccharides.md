# Get Available Monosaacharides

This function returns a character vector of monosaccharide names of the
given type. See
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/reference/get_mono_type.md)
for monosaacharide types. Concrete furanose forms use an `f` after the
monosaccharide stem, such as `Galf` and `GlcfNAc`. Generic names do not
encode ring form. Less common absolute configurations use a leading `D-`
or `L-`, such as `D-Fuc`, `L-Gul`, and `D-Fucf`. Unprefixed names retain
their natural configurations.

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
#>   [1] "Hex"         "HexNAc"      "HexN"        "HexA"        "dHex"       
#>   [6] "dHexNAc"     "ddHex"       "Pen"         "NeuAc"       "NeuGc"      
#>  [11] "gNeu"        "gKdn"        "gPse"        "gLeg"        "gAci"       
#>  [16] "g4eLeg"      "gBac"        "Hep"         "gKdo"        "HepA"       
#>  [21] "MurAc"       "MurGc"       "gMur"        "Glc"         "Man"        
#>  [26] "Gal"         "Gul"         "Alt"         "All"         "Tal"        
#>  [31] "Ido"         "GlcNAc"      "GalNAc"      "ManNAc"      "GulNAc"     
#>  [36] "AltNAc"      "AllNAc"      "TalNAc"      "IdoNAc"      "GlcN"       
#>  [41] "ManN"        "GalN"        "GulN"        "AltN"        "AllN"       
#>  [46] "TalN"        "IdoN"        "GlcA"        "ManA"        "GalA"       
#>  [51] "GulA"        "AltA"        "AllA"        "TalA"        "IdoA"       
#>  [56] "Fuc"         "Qui"         "Rha"         "6dGul"       "6dAlt"      
#>  [61] "6dTal"       "QuiNAc"      "RhaNAc"      "6dAltNAc"    "6dTalNAc"   
#>  [66] "FucNAc"      "Oli"         "Tyv"         "Abe"         "Par"        
#>  [71] "Dig"         "Col"         "Ara"         "Lyx"         "Xyl"        
#>  [76] "Rib"         "Neu5Ac"      "Neu5Gc"      "Neu"         "Kdn"        
#>  [81] "Pse"         "Leg"         "Aci"         "4eLeg"       "Bac"        
#>  [86] "LDmanHep"    "Kdo"         "Dha"         "DDmanHep"    "MurNAc"     
#>  [91] "MurNGc"      "Mur"         "Api"         "Fru"         "Tag"        
#>  [96] "Sor"         "Psi"         "Glcf"        "Manf"        "Galf"       
#> [101] "Gulf"        "Altf"        "Allf"        "Talf"        "Idof"       
#> [106] "GlcfNAc"     "GalfNAc"     "ManfNAc"     "GulfNAc"     "AltfNAc"    
#> [111] "AllfNAc"     "TalfNAc"     "IdofNAc"     "GlcfN"       "ManfN"      
#> [116] "GalfN"       "GulfN"       "AltfN"       "AllfN"       "TalfN"      
#> [121] "IdofN"       "GlcfA"       "ManfA"       "GalfA"       "GulfA"      
#> [126] "AltfA"       "AllfA"       "TalfA"       "IdofA"       "Fucf"       
#> [131] "Quif"        "Rhaf"        "6dGulf"      "6dAltf"      "6dTalf"     
#> [136] "QuifNAc"     "RhafNAc"     "6dAltfNAc"   "6dTalfNAc"   "FucfNAc"    
#> [141] "Olif"        "Tyvf"        "Abef"        "Parf"        "Digf"       
#> [146] "Colf"        "Araf"        "Lyxf"        "Xylf"        "Ribf"       
#> [151] "Neuf5Ac"     "Neuf5Gc"     "Neuf"        "Kdnf"        "Psef"       
#> [156] "Legf"        "Acif"        "4eLegf"      "Bacf"        "LDmanHepf"  
#> [161] "Kdof"        "Dhaf"        "DDmanHepf"   "MurfNAc"     "MurfNGc"    
#> [166] "Murf"        "Apif"        "Fruf"        "Tagf"        "Sorf"       
#> [171] "Psif"        "L-Glc"       "L-Man"       "L-Gal"       "L-Gul"      
#> [176] "D-Alt"       "L-All"       "L-Tal"       "D-Ido"       "L-GlcNAc"   
#> [181] "L-GalNAc"    "L-ManNAc"    "L-GulNAc"    "D-AltNAc"    "L-AllNAc"   
#> [186] "L-TalNAc"    "D-IdoNAc"    "L-GlcN"      "L-ManN"      "L-GalN"     
#> [191] "L-GulN"      "D-AltN"      "L-AllN"      "L-TalN"      "D-IdoN"     
#> [196] "L-GlcA"      "L-ManA"      "L-GalA"      "L-GulA"      "D-AltA"     
#> [201] "L-AllA"      "L-TalA"      "D-IdoA"      "D-Fuc"       "L-Qui"      
#> [206] "D-Rha"       "L-6dGul"     "D-6dAlt"     "L-6dTal"     "L-QuiNAc"   
#> [211] "D-RhaNAc"    "D-6dAltNAc"  "L-6dTalNAc"  "D-FucNAc"    "L-Oli"      
#> [216] "L-Tyv"       "L-Abe"       "L-Par"       "L-Dig"       "D-Col"      
#> [221] "D-Ara"       "L-Lyx"       "L-Xyl"       "L-Rib"       "L-Neu5Ac"   
#> [226] "L-Neu5Gc"    "L-Kdn"       "L-Bac"       "L-Kdo"       "L-Dha"      
#> [231] "L-MurNAc"    "L-MurNGc"    "L-Mur"       "D-Api"       "L-Fru"      
#> [236] "L-Tag"       "D-Sor"       "L-Psi"       "L-Glcf"      "L-Manf"     
#> [241] "L-Galf"      "L-Gulf"      "D-Altf"      "L-Allf"      "L-Talf"     
#> [246] "D-Idof"      "L-GlcfNAc"   "L-GalfNAc"   "L-ManfNAc"   "L-GulfNAc"  
#> [251] "D-AltfNAc"   "L-AllfNAc"   "L-TalfNAc"   "D-IdofNAc"   "L-GlcfN"    
#> [256] "L-ManfN"     "L-GalfN"     "L-GulfN"     "D-AltfN"     "L-AllfN"    
#> [261] "L-TalfN"     "D-IdofN"     "L-GlcfA"     "L-ManfA"     "L-GalfA"    
#> [266] "L-GulfA"     "D-AltfA"     "L-AllfA"     "L-TalfA"     "D-IdofA"    
#> [271] "D-Fucf"      "L-Quif"      "D-Rhaf"      "L-6dGulf"    "D-6dAltf"   
#> [276] "L-6dTalf"    "L-QuifNAc"   "D-RhafNAc"   "D-6dAltfNAc" "L-6dTalfNAc"
#> [281] "D-FucfNAc"   "L-Olif"      "L-Tyvf"      "L-Abef"      "L-Parf"     
#> [286] "L-Digf"      "D-Colf"      "D-Araf"      "L-Lyxf"      "L-Xylf"     
#> [291] "L-Ribf"      "L-Neuf5Ac"   "L-Neuf5Gc"   "L-Kdnf"      "L-Bacf"     
#> [296] "L-Kdof"      "L-Dhaf"      "L-MurfNAc"   "L-MurfNGc"   "L-Murf"     
#> [301] "D-Apif"      "L-Fruf"      "L-Tagf"      "D-Sorf"      "L-Psif"     
available_monosaccharides("concrete")
#>   [1] "Glc"         "Man"         "Gal"         "Gul"         "Alt"        
#>   [6] "All"         "Tal"         "Ido"         "GlcNAc"      "GalNAc"     
#>  [11] "ManNAc"      "GulNAc"      "AltNAc"      "AllNAc"      "TalNAc"     
#>  [16] "IdoNAc"      "GlcN"        "ManN"        "GalN"        "GulN"       
#>  [21] "AltN"        "AllN"        "TalN"        "IdoN"        "GlcA"       
#>  [26] "ManA"        "GalA"        "GulA"        "AltA"        "AllA"       
#>  [31] "TalA"        "IdoA"        "Fuc"         "Qui"         "Rha"        
#>  [36] "6dGul"       "6dAlt"       "6dTal"       "QuiNAc"      "RhaNAc"     
#>  [41] "6dAltNAc"    "6dTalNAc"    "FucNAc"      "Oli"         "Tyv"        
#>  [46] "Abe"         "Par"         "Dig"         "Col"         "Ara"        
#>  [51] "Lyx"         "Xyl"         "Rib"         "Neu5Ac"      "Neu5Gc"     
#>  [56] "Neu"         "Kdn"         "Pse"         "Leg"         "Aci"        
#>  [61] "4eLeg"       "Bac"         "LDmanHep"    "Kdo"         "Dha"        
#>  [66] "DDmanHep"    "MurNAc"      "MurNGc"      "Mur"         "Api"        
#>  [71] "Fru"         "Tag"         "Sor"         "Psi"         "Glcf"       
#>  [76] "Manf"        "Galf"        "Gulf"        "Altf"        "Allf"       
#>  [81] "Talf"        "Idof"        "GlcfNAc"     "GalfNAc"     "ManfNAc"    
#>  [86] "GulfNAc"     "AltfNAc"     "AllfNAc"     "TalfNAc"     "IdofNAc"    
#>  [91] "GlcfN"       "ManfN"       "GalfN"       "GulfN"       "AltfN"      
#>  [96] "AllfN"       "TalfN"       "IdofN"       "GlcfA"       "ManfA"      
#> [101] "GalfA"       "GulfA"       "AltfA"       "AllfA"       "TalfA"      
#> [106] "IdofA"       "Fucf"        "Quif"        "Rhaf"        "6dGulf"     
#> [111] "6dAltf"      "6dTalf"      "QuifNAc"     "RhafNAc"     "6dAltfNAc"  
#> [116] "6dTalfNAc"   "FucfNAc"     "Olif"        "Tyvf"        "Abef"       
#> [121] "Parf"        "Digf"        "Colf"        "Araf"        "Lyxf"       
#> [126] "Xylf"        "Ribf"        "Neuf5Ac"     "Neuf5Gc"     "Neuf"       
#> [131] "Kdnf"        "Psef"        "Legf"        "Acif"        "4eLegf"     
#> [136] "Bacf"        "LDmanHepf"   "Kdof"        "Dhaf"        "DDmanHepf"  
#> [141] "MurfNAc"     "MurfNGc"     "Murf"        "Apif"        "Fruf"       
#> [146] "Tagf"        "Sorf"        "Psif"        "L-Glc"       "L-Man"      
#> [151] "L-Gal"       "L-Gul"       "D-Alt"       "L-All"       "L-Tal"      
#> [156] "D-Ido"       "L-GlcNAc"    "L-GalNAc"    "L-ManNAc"    "L-GulNAc"   
#> [161] "D-AltNAc"    "L-AllNAc"    "L-TalNAc"    "D-IdoNAc"    "L-GlcN"     
#> [166] "L-ManN"      "L-GalN"      "L-GulN"      "D-AltN"      "L-AllN"     
#> [171] "L-TalN"      "D-IdoN"      "L-GlcA"      "L-ManA"      "L-GalA"     
#> [176] "L-GulA"      "D-AltA"      "L-AllA"      "L-TalA"      "D-IdoA"     
#> [181] "D-Fuc"       "L-Qui"       "D-Rha"       "L-6dGul"     "D-6dAlt"    
#> [186] "L-6dTal"     "L-QuiNAc"    "D-RhaNAc"    "D-6dAltNAc"  "L-6dTalNAc" 
#> [191] "D-FucNAc"    "L-Oli"       "L-Tyv"       "L-Abe"       "L-Par"      
#> [196] "L-Dig"       "D-Col"       "D-Ara"       "L-Lyx"       "L-Xyl"      
#> [201] "L-Rib"       "L-Neu5Ac"    "L-Neu5Gc"    "L-Kdn"       "L-Bac"      
#> [206] "L-Kdo"       "L-Dha"       "L-MurNAc"    "L-MurNGc"    "L-Mur"      
#> [211] "D-Api"       "L-Fru"       "L-Tag"       "D-Sor"       "L-Psi"      
#> [216] "L-Glcf"      "L-Manf"      "L-Galf"      "L-Gulf"      "D-Altf"     
#> [221] "L-Allf"      "L-Talf"      "D-Idof"      "L-GlcfNAc"   "L-GalfNAc"  
#> [226] "L-ManfNAc"   "L-GulfNAc"   "D-AltfNAc"   "L-AllfNAc"   "L-TalfNAc"  
#> [231] "D-IdofNAc"   "L-GlcfN"     "L-ManfN"     "L-GalfN"     "L-GulfN"    
#> [236] "D-AltfN"     "L-AllfN"     "L-TalfN"     "D-IdofN"     "L-GlcfA"    
#> [241] "L-ManfA"     "L-GalfA"     "L-GulfA"     "D-AltfA"     "L-AllfA"    
#> [246] "L-TalfA"     "D-IdofA"     "D-Fucf"      "L-Quif"      "D-Rhaf"     
#> [251] "L-6dGulf"    "D-6dAltf"    "L-6dTalf"    "L-QuifNAc"   "D-RhafNAc"  
#> [256] "D-6dAltfNAc" "L-6dTalfNAc" "D-FucfNAc"   "L-Olif"      "L-Tyvf"     
#> [261] "L-Abef"      "L-Parf"      "L-Digf"      "D-Colf"      "D-Araf"     
#> [266] "L-Lyxf"      "L-Xylf"      "L-Ribf"      "L-Neuf5Ac"   "L-Neuf5Gc"  
#> [271] "L-Kdnf"      "L-Bacf"      "L-Kdof"      "L-Dhaf"      "L-MurfNAc"  
#> [276] "L-MurfNGc"   "L-Murf"      "D-Apif"      "L-Fruf"      "L-Tagf"     
#> [281] "D-Sorf"      "L-Psif"     
```
