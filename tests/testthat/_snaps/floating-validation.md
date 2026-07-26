# floating metadata rejects unsupported fields

    Code
      glycan_structure(graph)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part metadata contains unsupported fields.
      x Unsupported field: probability.

# floating metadata rejects duplicated fields

    Code
      glycan_structure(graph)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part metadata contains duplicated fields.
      x Duplicated field: root.

# floating metadata requires every core field

    Code
      glycan_structure(graph)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part metadata is incomplete.
      x Missing field: parents.

# floating node metadata must describe exactly one component

    Code
      glycan_structure(duplicate_nodes)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part node indices must be unique.

---

    Code
      glycan_structure(missing_root)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part nodes must contain its root.

---

    Code
      glycan_structure(wrong_component)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `validate_floating_graph_shape()`:
      ! Floating part node metadata does not match its graph component.
      x Floating part 1 must contain exactly node 2.

# explicit parents reject definitely occupied acceptor slots

    Code
      glycan_structure(graph)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `floating_attachment_domain()`:
      ! Floating part has impossible explicit parent metadata.
      x Floating part 1 cannot use explicit parent node 1 because every acceptor position declared by linkage "a2-3" is already occupied.

# floating parts require a simultaneous attachment assignment

    Code
      glycan_structure(conflict)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `validate_floating_attachment_slots()`:
      ! Floating parts cannot be attached simultaneously.
      x No conflict-free assignment exists for the declared parent and acceptor positions.

# ambiguous main edges participate in slot matching

    Code
      glycan_structure(conflict)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `validate_floating_attachment_slots()`:
      ! Floating parts cannot be attached simultaneously.
      x No conflict-free assignment exists for the declared parent and acceptor positions.

# unrestricted parents use all available main-tree slots

    Code
      glycan_structure(conflict)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `validate_floating_attachment_slots()`:
      ! Floating parts cannot be attached simultaneously.
      x No conflict-free assignment exists for the declared parent and acceptor positions.

