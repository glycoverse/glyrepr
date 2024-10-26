# print works for NE glycan graphs

    Code
      print(x)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3

# print works for DN glycan graphs

    Code
      print(x)
    Output
      Glycan Graph (DN)
      GlcNAc: 2, Man: 3

# verbose print works for NE glycan graphs

    Code
      print(x, verbose = TRUE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# verbose print works for DN glycan graphs

    Code
      print(x, verbose = TRUE)
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

# verbose print works for NE glycan graphs without linkages

    Code
      print(x, verbose = TRUE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc
      └─GlcNAc (??-?)
        └─Man (??-?)
          ├─Man (??-?)
          └─Man (??-?)

# verbose print works for DN glycan graphs without linkages

    Code
      print(x, verbose = TRUE)
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

# verbose print works for one-mono NE glycan graphs

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N

# verbose print works for one-mono DN glycan graphs

    Code
      print(glycan, verbose = TRUE)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N

