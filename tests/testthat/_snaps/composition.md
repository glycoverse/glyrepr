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

# print.glyrepr_composition supports n

    <glycan_composition[11]>
    [1] Hex(5)
    [2] Hex(5)
    [3] Hex(5)
    [4] Hex(5)
    [5] Hex(5)
    [6] Hex(5)
    [7] Hex(5)
    [8] Hex(5)
    [9] Hex(5)
    [10] Hex(5)
    ... (1 more not shown)

---

    <glycan_composition[11]>
    [1] Hex(5)
    [2] Hex(5)
    [3] Hex(5)
    [4] Hex(5)
    [5] Hex(5)
    [6] Hex(5)
    [7] Hex(5)
    [8] Hex(5)
    [9] Hex(5)
    [10] Hex(5)
    [11] Hex(5)

