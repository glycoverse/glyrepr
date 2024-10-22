# converting NE to DN graph works on simple example

    Code
      print(dn_graph, verbose = TRUE)
    Output
      Glycan Graph (DN)
      H: 1, N: 2
      ------------------
      N
      ├─b1-4
      │ └─N
      └─a1-3
        └─H

# converting NE to DN graph works on complex example

    Code
      print(dn_graph, verbose = TRUE)
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

# converting one-node graph to DN graph works

    Code
      print(dn_graph, verbose = TRUE)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N

