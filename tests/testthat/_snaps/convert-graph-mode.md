# converting NE to DN graph works on simple example

    Code
      print(dn_graph, verbose = TRUE)
    Output
      Glycan Graph (DN)
      H: 1, N: 2
      ------------------
      N-S
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
      GlcNAc-S
      └─b1-4
        └─GlcNAc
          └─b1-4
            └─Man
              ├─a1-3
              │ └─Man
              └─a1-6
                └─Man

# converting to DN graph fails for DN graphs

    Code
      convert_graph_mode(glycan, "dn")
    Condition
      Error in `convert_graph_mode()`:
      ! The glycan is already in "DN" mode.

# converting one-node graph to DN graph works

    Code
      print(dn_graph, verbose = TRUE)
    Output
      Glycan Graph (DN)
      N: 1
      ------------------
      N

# converting DN to NE graph works on simple example

    Code
      print(ne_graph, verbose = TRUE)
    Output
      Glycan Graph (NE)
      H: 2, N: 1
      ------------------
      N-S
      ├─H (b1-4)
      └─H (a1-3)

# converting DN to NE graph works on complex example

    Code
      print(ne_graph, verbose = TRUE)
    Output
      Glycan Graph (NE)
      GlcNAc: 2, Man: 3
      ------------------
      GlcNAc-S
      └─GlcNAc (b1-4)
        └─Man (b1-4)
          ├─Man (a1-3)
          └─Man (a1-6)

# converting one-node graph to NE graph works

    Code
      print(ne_graph, verbose = TRUE)
    Output
      Glycan Graph (NE)
      N: 1
      ------------------
      N

