# format.glyrepr_structure includes names with tab separation

    Code
      format(glycans)
    Output
      [1] "A\tGal(b1-3)GalNAc(a1-                                "
      [2] "B\tMan(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

# format.glyrepr_structure without names works correctly

    Code
      format(glycans)
    Output
      [1] "Gal(b1-3)GalNAc(a1-                                "
      [2] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

# truncation works in tibble

    Code
      print(tibble, width = 30)
    Output
      # A tibble: 3 x 2
        struc                      a
        <struct>               <dbl>
      1 Man(a1-3)[Man(a1-6)]M~     1
      2 Man(a1-3)[Man(a1-6)]M~     1
      3 Man(a1-3)[Man(a1-6)]M~     1

