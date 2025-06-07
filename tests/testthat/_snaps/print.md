# print works for glycan graphs

    Code
      print(x)
    Output
      Glycan Graph
      Man(3)GlcNAc(2)
      ------------------
      GlcNAc (?1-)
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# print works for glycan graphs with verbose = FALSE

    Code
      print(x, verbose = FALSE)
    Output
      Glycan Graph
      Man(3)GlcNAc(2)

# print works for glycan graphs without linkages

    Code
      print(x)
    Output
      Glycan Graph
      Man(3)GlcNAc(2)
      ------------------
      GlcNAc (?1-)
      └─GlcNAc
        └─Man
          ├─Man
          └─Man

# print works for one-mono glycan graphs

    Code
      print(glycan)
    Output
      Glycan Graph
      N1
      ------------------
      N (a1-)

# print works for graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph
      Gal(1)GalNAc(1)
      ------------------
      GalNAc-6S (a1-)
      └─Gal (b1-3)

# print works for one-node graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph
      N1
      ------------------
      N-6S (a1-)

# print works for graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph
      Gal(1)GalNAc(1)
      ------------------
      GalNAc-ol (a1-)
      └─Gal (b1-3)

# print works for one-node graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph
      N1
      ------------------
      N-ol (a1-)

# print works for graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph
      Gal(1)GalNAc(1)
      ------------------
      GalNAc-ol-6S (a1-)
      └─Gal (b1-3)

# print works for one-node graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph
      N1
      ------------------
      N-ol-6S (a1-)

