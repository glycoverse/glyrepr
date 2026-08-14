# localize_floating_parts validates assignment tables

    Code
      localize_floating_parts(glycan, tibble::tibble())
    Condition
      Error in `validate_structure_table()`:
      ! `assignments` must contain the required columns.
      x Missing columns: glycan_id, part_id, and parent_node.

---

    Code
      localize_floating_parts(glycan, tibble::tibble(glycan_id = c(1L, 1L), part_id = c(
        1L, 1L), parent_node = c(1L, 2L)))
    Condition
      Error in `validate_floating_assignments()`:
      ! Each floating part can be assigned at most once.
      x Duplicated (glycan_id, part_id) pair: "(1, 1)".

---

    Code
      localize_floating_parts(glycan, tibble::tibble(glycan_id = 2L, part_id = 1L,
        parent_node = 1L))
    Condition
      Error in `validate_floating_assignments()`:
      ! `assignments` contains out-of-range glycan identifiers.
      x glycan_id must be between 1 and 1.

# localize_floating_parts rejects invalid targets

    Code
      localize_floating_parts(glycans, tibble::tibble(glycan_id = 1L, part_id = 1L,
        parent_node = 1L))
    Condition
      Error in `localize_floating_graph()`:
      ! Assignment selects a parent outside the candidate domain.
      x Glycan 1 floating part 1 cannot attach to parent node 1.
      i Candidate parent domain: 2 and 3.

---

    Code
      localize_floating_parts(glycans, tibble::tibble(glycan_id = 2L, part_id = 1L,
        parent_node = 1L))
    Condition
      Error in `localize_floating_parts()`:
      ! Cannot localize a missing glycan structure.
      x Glycan 2 is missing.

---

    Code
      localize_floating_parts(glycans, tibble::tibble(glycan_id = 3L, part_id = 1L,
        parent_node = 1L))
    Condition
      Error in `localize_floating_graph()`:
      ! Assignment refers to a nonexistent floating part.
      x Glycan 3 has no floating part with part_id 1.

# localize_floating_parts validates simultaneous slot conflicts

    Code
      localize_floating_parts(glycan, tibble::tibble(glycan_id = c(1L, 1L), part_id = c(
        1L, 2L), parent_node = c(3L, 3L)))
    Condition
      Error in `validate_floating_metadata_assignments()`:
      ! Floating parts cannot be attached simultaneously.
      x No conflict-free acyclic assignment connects every floating component to the main tree.

# localize_floating_parts checks occupied main-tree slots

    Code
      localize_floating_parts(glycan, tibble::tibble(glycan_id = 1L, part_id = 1L,
        parent_node = 3L))
    Condition
      Error in `floating_attachment_domain()`:
      ! Floating part has impossible explicit parent metadata.
      x Floating part 1 cannot use explicit parent node 3 because every acceptor position declared by linkage "a2-3" is already occupied.

# graph localization validates inputs and its conservative bound

    Code
      enumerate_floating_graph_localizations(graph, max_variants = 3)
    Condition
      Error in `enumerate_floating_graph_localizations()`:
      ! Floating localization count exceeds `max_variants`.
      x Input 1 has 4 raw candidate combinations; the limit is 3.
      i Increase `max_variants` to enumerate this input.

# enumerate_floating_localizations enforces a conservative bound

    Code
      enumerate_floating_localizations(glycan, max_variants = 3)
    Condition
      Error in `enumerate_floating_localizations_one()`:
      ! Floating localization count exceeds `max_variants`.
      x Input 1 has 4 raw candidate combinations; the limit is 3.
      i Increase `max_variants` to enumerate this input.

