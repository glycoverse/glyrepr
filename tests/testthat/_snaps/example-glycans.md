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
      └─GlcNAc
        └─Man
          ├─Man
          └─Man

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

