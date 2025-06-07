# N-glycan core graph

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      Man(3)GlcNAc(2)
      ------------------
      GlcNAc (?1-)
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# N-glycan core graph without linkages

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      Man(3)GlcNAc(2)
      ------------------
      GlcNAc (?1-)
      └─GlcNAc
        └─Man
          ├─Man
          └─Man

# N-glycan core graph with simple monosaccharides

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      H3N2
      ------------------
      N (?1-)
      └─N (b1-4)
        └─H (b1-4)
          ├─H (a1-3)
          └─H (a1-6)

# N-glycan core graph with generic monosaccharides

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      Hex(3)HexNAc(2)
      ------------------
      HexNAc (?1-)
      └─HexNAc (b1-4)
        └─Hex (b1-4)
          ├─Hex (a1-3)
          └─Hex (a1-6)

# O-glycan core 1

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      Gal(1)GalNAc(1)
      ------------------
      GalNAc (a1-)
      └─Gal (b1-3)

# O-glycan core 2

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph
      Gal(1)GlcNAc(1)GalNAc(1)
      ------------------
      GalNAc (a1-)
      ├─Gal (b1-3)
      └─GlcNAc (b1-6)

