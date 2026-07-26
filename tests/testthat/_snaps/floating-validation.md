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

# unrestricted parents use all available main-tree slots

    Code
      glycan_structure(conflict)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `validate_floating_attachment_slots()`:
      ! Floating parts cannot be attached simultaneously.
      x No conflict-free assignment exists for the declared parent and acceptor positions.

