# glyrepr 0.6.0

## Breaking changes

* Remove the "alditol" attribute from `glycan_structure()` objects. This information is rarely used in glycomics and glycoproteomics data analysis. It is removed according to the razor principle.
* `as_glycan_structure()` now doesn't allow the input IUPAC-condensed strings to omit the anomer information. Previously, something like "Glc(a1-3)GlcNAc" is valid. `as_glycan_structure()` assumed that the core "GlcNAc" has a "?1-" anomer and added it automatically. The problem is that this behavior was not easily awared by users and might cause confusion. Again, less is more, so we remove it.

# glyrepr 0.5.2

## Minor improvements and bug fixes

* Fix some error message format errors.
* Update the abbreviated type name of `glyrepr_structure` from `structure` to `struct`.

# glyrepr 0.5.0

## Major changes

* Glycan structures now support multiple substituents on a single monosaccharide.
  Substituents are stored as comma-separated strings internally and concatenated
  in IUPAC format for display.

* Glycan compositions now support substituents. The `glycan_composition` class
  can now represent and count substituents alongside monosaccharides.
