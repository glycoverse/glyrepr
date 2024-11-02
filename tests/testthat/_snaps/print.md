# print works for NE glycan graphs

    Code
      print(x)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc (?1-)
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# print works for DN glycan graphs

    Code
      print(x)
    Output
      Glycan Graph (DN)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc (?1-)
      └─b1-4
        └─GlcNAc
          └─b1-4
            └─Man
              ├─a1-3
              │ └─Man
              └─a1-6
                └─Man

# print works for NE glycan graphs with verbose = FALSE

    Code
      print(x, verbose = FALSE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3

# print works for DN glycan graphs with verbose = FALSE

    Code
      print(x, verbose = FALSE)
    Output
      Glycan Graph (DN)
      GlcNAc: 2, Man: 3

# print works for NE glycan graphs without linkages

    Code
      print(x)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc (?1-)
      └─GlcNAc
        └─Man
          ├─Man
          └─Man

# print works for DN glycan graphs without linkages

    Code
      print(x)
    Output
      Glycan Graph (DN)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc (?1-)
      └─??-?
        └─GlcNAc
          └─??-?
            └─Man
              ├─??-?
              │ └─Man
              └─??-?
                └─Man

# print works for one-mono NE glycan graphs

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N (a1-)

# print works for one-mono DN glycan graphs

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N (a1-)

# print works for NE graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-6S (a1-)
      └─Gal-6S (b1-3)

# print works for DN graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-6S (a1-)
      └─b1-3
        └─Gal-6S

# print works for one-node NE graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N-6S (a1-)

# print works for one-node DN graphs with substituent

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N-6S (a1-)

# print works for NE graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-ol (a1-)
      └─Gal (b1-3)

# print works for DN graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-ol (a1-)
      └─b1-3
        └─Gal

# print works for one-node NE graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N-ol (a1-)

# print works for one-node DN graphs with alditol

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N-ol (a1-)

# print works for NE graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-ol-6S (a1-)
      └─Gal (b1-3)

# print works for DN graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      Gal: 1, GalNAc: 1
      ------------------
      GalNAc-ol-6S (a1-)
      └─b1-3
        └─Gal

# print works for one-node NE graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N-ol-6S (a1-)

# print works for one-node DN graph with alditol and substituent on root

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N-ol-6S (a1-)

