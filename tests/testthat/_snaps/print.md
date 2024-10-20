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

