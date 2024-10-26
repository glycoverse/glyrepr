# print works for NE glycan graphs

    Code
      print(x)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
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
      GlcNAc
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
      GlcNAc
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
      GlcNAc
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
      N

# print works for one-mono DN glycan graphs

    Code
      print(glycan)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N

