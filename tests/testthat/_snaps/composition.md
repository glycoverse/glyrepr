# truncation works in tibble for compositions

    Code
      print(tibble, width = 30)
    Output
      # A tibble: 3 x 2
        comp                a
        <comp>          <dbl>
      1 Hex(1)HexNAc(1)     1
      2 Hex(2)HexNAc(1)     1
      3 Hex(3)dHex(1)       1

# as_glycan_composition works for E and L

    Code
      comps <- as_glycan_composition(chars)
    Condition
      Warning:
      Simple composition codes "E" and/or "L" are parsed as "NeuAc".
      i Linkage-specific Neu5Ac information is discarded in glycan compositions.

