# truncation works in tibble for compositions

    Code
      print(tibble, width = 30)
    Output
      # A tibble: 3 x 2
        comp                a
        <comp>          <dbl>
      1 Hex(2)HexNAc(1)     1
      2 Glc(1)Gal(2)        1
      3 Hex(3)dHex(1)       1

