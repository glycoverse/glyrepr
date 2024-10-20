# N-glycan core NE graph

    Code
      str(igraph::edge_attr(glycan))
    Output
      List of 1
       $ linkage: chr [1:4] "b1-4" "b1-4" "a1-3" "a1-6"

---

    Code
      str(igraph::vertex_attr(glycan))
    Output
      List of 2
       $ name: chr [1:5] "1" "2" "3" "4" ...
       $ mono: chr [1:5] "GlcNAc" "GlcNAc" "Man" "Man" ...

# N-glycan core DN graph

    Code
      str(igraph::edge_attr(glycan))
    Output
       Named list()

---

    Code
      str(igraph::vertex_attr(glycan))
    Output
      List of 4
       $ name   : chr [1:9] "1" "2" "3" "4" ...
       $ type   : chr [1:9] "mono" "linkage" "mono" "linkage" ...
       $ mono   : chr [1:9] "GlcNAc" NA "GlcNAc" NA ...
       $ linkage: chr [1:9] NA "b1-4" NA "b1-4" ...

