# Get Available Monosaacharides

This function returns a character vector of monosaccharide names of the
given type. See
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)
for monosaacharide types. Concrete furanose forms use an `f` after the
monosaccharide stem, such as `Galf` and `GlcfNAc`. Generic names do not
encode ring form. Less common absolute configurations use a leading `D`
or `L` without a separator, such as `DFuc`, `LGul`, and `DFucf`.
Unprefixed names retain their natural configurations.

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
#>   [1] "Hex"        "HexNAc"     "HexN"       "HexA"       "dHex"      
#>   [6] "dHexNAc"    "ddHex"      "Pen"        "NeuAc"      "NeuGc"     
#>  [11] "gNeu"       "gKdn"       "gPse"       "gLeg"       "gAci"      
#>  [16] "g4eLeg"     "gBac"       "Hep"        "gKdo"       "HepA"      
#>  [21] "MurAc"      "MurGc"      "gMur"       "Glc"        "Man"       
#>  [26] "Gal"        "Gul"        "Alt"        "All"        "Tal"       
#>  [31] "Ido"        "GlcNAc"     "GalNAc"     "ManNAc"     "GulNAc"    
#>  [36] "AltNAc"     "AllNAc"     "TalNAc"     "IdoNAc"     "GlcN"      
#>  [41] "ManN"       "GalN"       "GulN"       "AltN"       "AllN"      
#>  [46] "TalN"       "IdoN"       "GlcA"       "ManA"       "GalA"      
#>  [51] "GulA"       "AltA"       "AllA"       "TalA"       "IdoA"      
#>  [56] "Fuc"        "Qui"        "Rha"        "6dGul"      "6dAlt"     
#>  [61] "6dTal"      "QuiNAc"     "RhaNAc"     "6dAltNAc"   "6dTalNAc"  
#>  [66] "FucNAc"     "Oli"        "Tyv"        "Abe"        "Par"       
#>  [71] "Dig"        "Col"        "Ara"        "Lyx"        "Xyl"       
#>  [76] "Rib"        "Neu5Ac"     "Neu5Gc"     "Neu"        "Kdn"       
#>  [81] "Pse"        "Leg"        "Aci"        "4eLeg"      "Bac"       
#>  [86] "LDmanHep"   "Kdo"        "Dha"        "DDmanHep"   "MurNAc"    
#>  [91] "MurNGc"     "Mur"        "Api"        "Fru"        "Tag"       
#>  [96] "Sor"        "Psi"        "Glcf"       "Manf"       "Galf"      
#> [101] "Gulf"       "Altf"       "Allf"       "Talf"       "Idof"      
#> [106] "GlcfNAc"    "GalfNAc"    "ManfNAc"    "GulfNAc"    "AltfNAc"   
#> [111] "AllfNAc"    "TalfNAc"    "IdofNAc"    "GlcfN"      "ManfN"     
#> [116] "GalfN"      "GulfN"      "AltfN"      "AllfN"      "TalfN"     
#> [121] "IdofN"      "GlcfA"      "ManfA"      "GalfA"      "GulfA"     
#> [126] "AltfA"      "AllfA"      "TalfA"      "IdofA"      "Fucf"      
#> [131] "Quif"       "Rhaf"       "6dGulf"     "6dAltf"     "6dTalf"    
#> [136] "QuifNAc"    "RhafNAc"    "6dAltfNAc"  "6dTalfNAc"  "FucfNAc"   
#> [141] "Olif"       "Tyvf"       "Abef"       "Parf"       "Digf"      
#> [146] "Colf"       "Araf"       "Lyxf"       "Xylf"       "Ribf"      
#> [151] "Neuf5Ac"    "Neuf5Gc"    "Neuf"       "Kdnf"       "Psef"      
#> [156] "Legf"       "Acif"       "4eLegf"     "Bacf"       "LDmanHepf" 
#> [161] "Kdof"       "Dhaf"       "DDmanHepf"  "MurfNAc"    "MurfNGc"   
#> [166] "Murf"       "Apif"       "Fruf"       "Tagf"       "Sorf"      
#> [171] "Psif"       "LGlc"       "LMan"       "LGal"       "LGul"      
#> [176] "DAlt"       "LAll"       "LTal"       "DIdo"       "LGlcNAc"   
#> [181] "LGalNAc"    "LManNAc"    "LGulNAc"    "DAltNAc"    "LAllNAc"   
#> [186] "LTalNAc"    "DIdoNAc"    "LGlcN"      "LManN"      "LGalN"     
#> [191] "LGulN"      "DAltN"      "LAllN"      "LTalN"      "DIdoN"     
#> [196] "LGlcA"      "LManA"      "LGalA"      "LGulA"      "DAltA"     
#> [201] "LAllA"      "LTalA"      "DIdoA"      "DFuc"       "LQui"      
#> [206] "DRha"       "L6dGul"     "D6dAlt"     "L6dTal"     "LQuiNAc"   
#> [211] "DRhaNAc"    "D6dAltNAc"  "L6dTalNAc"  "DFucNAc"    "LOli"      
#> [216] "LTyv"       "LAbe"       "LPar"       "LDig"       "DCol"      
#> [221] "DAra"       "LLyx"       "LXyl"       "LRib"       "LNeu5Ac"   
#> [226] "LNeu5Gc"    "LKdn"       "LBac"       "LKdo"       "LDha"      
#> [231] "LMurNAc"    "LMurNGc"    "LMur"       "DApi"       "LFru"      
#> [236] "LTag"       "DSor"       "LPsi"       "LGlcf"      "LManf"     
#> [241] "LGalf"      "LGulf"      "DAltf"      "LAllf"      "LTalf"     
#> [246] "DIdof"      "LGlcfNAc"   "LGalfNAc"   "LManfNAc"   "LGulfNAc"  
#> [251] "DAltfNAc"   "LAllfNAc"   "LTalfNAc"   "DIdofNAc"   "LGlcfN"    
#> [256] "LManfN"     "LGalfN"     "LGulfN"     "DAltfN"     "LAllfN"    
#> [261] "LTalfN"     "DIdofN"     "LGlcfA"     "LManfA"     "LGalfA"    
#> [266] "LGulfA"     "DAltfA"     "LAllfA"     "LTalfA"     "DIdofA"    
#> [271] "DFucf"      "LQuif"      "DRhaf"      "L6dGulf"    "D6dAltf"   
#> [276] "L6dTalf"    "LQuifNAc"   "DRhafNAc"   "D6dAltfNAc" "L6dTalfNAc"
#> [281] "DFucfNAc"   "LOlif"      "LTyvf"      "LAbef"      "LParf"     
#> [286] "LDigf"      "DColf"      "DAraf"      "LLyxf"      "LXylf"     
#> [291] "LRibf"      "LNeuf5Ac"   "LNeuf5Gc"   "LKdnf"      "LBacf"     
#> [296] "LKdof"      "LDhaf"      "LMurfNAc"   "LMurfNGc"   "LMurf"     
#> [301] "DApif"      "LFruf"      "LTagf"      "DSorf"      "LPsif"     
available_monosaccharides("concrete")
#>   [1] "Glc"        "Man"        "Gal"        "Gul"        "Alt"       
#>   [6] "All"        "Tal"        "Ido"        "GlcNAc"     "GalNAc"    
#>  [11] "ManNAc"     "GulNAc"     "AltNAc"     "AllNAc"     "TalNAc"    
#>  [16] "IdoNAc"     "GlcN"       "ManN"       "GalN"       "GulN"      
#>  [21] "AltN"       "AllN"       "TalN"       "IdoN"       "GlcA"      
#>  [26] "ManA"       "GalA"       "GulA"       "AltA"       "AllA"      
#>  [31] "TalA"       "IdoA"       "Fuc"        "Qui"        "Rha"       
#>  [36] "6dGul"      "6dAlt"      "6dTal"      "QuiNAc"     "RhaNAc"    
#>  [41] "6dAltNAc"   "6dTalNAc"   "FucNAc"     "Oli"        "Tyv"       
#>  [46] "Abe"        "Par"        "Dig"        "Col"        "Ara"       
#>  [51] "Lyx"        "Xyl"        "Rib"        "Neu5Ac"     "Neu5Gc"    
#>  [56] "Neu"        "Kdn"        "Pse"        "Leg"        "Aci"       
#>  [61] "4eLeg"      "Bac"        "LDmanHep"   "Kdo"        "Dha"       
#>  [66] "DDmanHep"   "MurNAc"     "MurNGc"     "Mur"        "Api"       
#>  [71] "Fru"        "Tag"        "Sor"        "Psi"        "Glcf"      
#>  [76] "Manf"       "Galf"       "Gulf"       "Altf"       "Allf"      
#>  [81] "Talf"       "Idof"       "GlcfNAc"    "GalfNAc"    "ManfNAc"   
#>  [86] "GulfNAc"    "AltfNAc"    "AllfNAc"    "TalfNAc"    "IdofNAc"   
#>  [91] "GlcfN"      "ManfN"      "GalfN"      "GulfN"      "AltfN"     
#>  [96] "AllfN"      "TalfN"      "IdofN"      "GlcfA"      "ManfA"     
#> [101] "GalfA"      "GulfA"      "AltfA"      "AllfA"      "TalfA"     
#> [106] "IdofA"      "Fucf"       "Quif"       "Rhaf"       "6dGulf"    
#> [111] "6dAltf"     "6dTalf"     "QuifNAc"    "RhafNAc"    "6dAltfNAc" 
#> [116] "6dTalfNAc"  "FucfNAc"    "Olif"       "Tyvf"       "Abef"      
#> [121] "Parf"       "Digf"       "Colf"       "Araf"       "Lyxf"      
#> [126] "Xylf"       "Ribf"       "Neuf5Ac"    "Neuf5Gc"    "Neuf"      
#> [131] "Kdnf"       "Psef"       "Legf"       "Acif"       "4eLegf"    
#> [136] "Bacf"       "LDmanHepf"  "Kdof"       "Dhaf"       "DDmanHepf" 
#> [141] "MurfNAc"    "MurfNGc"    "Murf"       "Apif"       "Fruf"      
#> [146] "Tagf"       "Sorf"       "Psif"       "LGlc"       "LMan"      
#> [151] "LGal"       "LGul"       "DAlt"       "LAll"       "LTal"      
#> [156] "DIdo"       "LGlcNAc"    "LGalNAc"    "LManNAc"    "LGulNAc"   
#> [161] "DAltNAc"    "LAllNAc"    "LTalNAc"    "DIdoNAc"    "LGlcN"     
#> [166] "LManN"      "LGalN"      "LGulN"      "DAltN"      "LAllN"     
#> [171] "LTalN"      "DIdoN"      "LGlcA"      "LManA"      "LGalA"     
#> [176] "LGulA"      "DAltA"      "LAllA"      "LTalA"      "DIdoA"     
#> [181] "DFuc"       "LQui"       "DRha"       "L6dGul"     "D6dAlt"    
#> [186] "L6dTal"     "LQuiNAc"    "DRhaNAc"    "D6dAltNAc"  "L6dTalNAc" 
#> [191] "DFucNAc"    "LOli"       "LTyv"       "LAbe"       "LPar"      
#> [196] "LDig"       "DCol"       "DAra"       "LLyx"       "LXyl"      
#> [201] "LRib"       "LNeu5Ac"    "LNeu5Gc"    "LKdn"       "LBac"      
#> [206] "LKdo"       "LDha"       "LMurNAc"    "LMurNGc"    "LMur"      
#> [211] "DApi"       "LFru"       "LTag"       "DSor"       "LPsi"      
#> [216] "LGlcf"      "LManf"      "LGalf"      "LGulf"      "DAltf"     
#> [221] "LAllf"      "LTalf"      "DIdof"      "LGlcfNAc"   "LGalfNAc"  
#> [226] "LManfNAc"   "LGulfNAc"   "DAltfNAc"   "LAllfNAc"   "LTalfNAc"  
#> [231] "DIdofNAc"   "LGlcfN"     "LManfN"     "LGalfN"     "LGulfN"    
#> [236] "DAltfN"     "LAllfN"     "LTalfN"     "DIdofN"     "LGlcfA"    
#> [241] "LManfA"     "LGalfA"     "LGulfA"     "DAltfA"     "LAllfA"    
#> [246] "LTalfA"     "DIdofA"     "DFucf"      "LQuif"      "DRhaf"     
#> [251] "L6dGulf"    "D6dAltf"    "L6dTalf"    "LQuifNAc"   "DRhafNAc"  
#> [256] "D6dAltfNAc" "L6dTalfNAc" "DFucfNAc"   "LOlif"      "LTyvf"     
#> [261] "LAbef"      "LParf"      "LDigf"      "DColf"      "DAraf"     
#> [266] "LLyxf"      "LXylf"      "LRibf"      "LNeuf5Ac"   "LNeuf5Gc"  
#> [271] "LKdnf"      "LBacf"      "LKdof"      "LDhaf"      "LMurfNAc"  
#> [276] "LMurfNGc"   "LMurf"      "DApif"      "LFruf"      "LTagf"     
#> [281] "DSorf"      "LPsif"     
```
