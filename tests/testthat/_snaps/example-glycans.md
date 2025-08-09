# N-glycan core structure

    Code
      print(glycan, verbose = TRUE)
    Output
      <glycan_structure[1]>
      [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
      # Unique structures: 1

# N-glycan core structure without linkages

    Code
      print(glycan, verbose = TRUE)
    Output
      <glycan_structure[1]>
      [1] Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
      # Unique structures: 1

# N-glycan core structure with generic monosaccharides

    Code
      print(glycan, verbose = TRUE)
    Output
      <glycan_structure[1]>
      [1] Hex(a1-3)[Hex(a1-6)]Hex(b1-4)HexNAc(b1-4)HexNAc(b1-
      # Unique structures: 1

# O-glycan core 1

    Code
      print(glycan, verbose = TRUE)
    Output
      <glycan_structure[1]>
      [1] Gal(b1-3)GalNAc(a1-
      # Unique structures: 1

# O-glycan core 2

    Code
      print(glycan, verbose = TRUE)
    Output
      <glycan_structure[1]>
      [1] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
      # Unique structures: 1

