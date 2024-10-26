# N-glycan core NE graph

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# N-glycan core DN graph

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (DN)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
      └─b1-4
        └─GlcNAc
          └─b1-4
            └─Man
              ├─a1-3
              │ └─Man
              └─a1-6
                └─Man

# N-glycan core NE graph without linkages

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
      └─GlcNAc (??-?)
        └─Man (??-?)
          ├─Man (??-?)
          └─Man (??-?)

# N-glycan core NE graph with simple monosaccharides

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      H: 3, N: 2
      ------------------
      N
      └─N (b1-4)
        └─H (b1-4)
          ├─H (a1-3)
          └─H (a1-6)

# N-glycan core NE graph with generic monosaccharides

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      Hex: 3, HexNAc: 2
      ------------------
      HexNAc
      └─HexNAc (b1-4)
        └─Hex (b1-4)
          ├─Hex (a1-3)
          └─Hex (a1-6)

# O-glycan core 1

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc
      └─Gal (b1-3)

# O-glycan core 2

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      Gal: 1, GalNAc: 1, GlcNAc: 1
      ------------------
      GalNAc
      ├─Gal (b1-3)
      └─GlcNAc (b1-6)

