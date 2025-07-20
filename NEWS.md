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
