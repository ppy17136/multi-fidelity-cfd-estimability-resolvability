# Version 2.1.9

This release corrects numerical certificate boundary handling in 2.1.8.

- Two-label minimum cuts use exact integer capacities derived from stored
  binary edge costs, with strict forcing and independent label/cost checks.
- Scalar baseline and independent-endpoint comparisons use exact rational
  representations of stored values. Impossible crossings are rejected before
  MILP execution; non-crossing or inconsistent solver witnesses raise an error.
- An initially indeterminate scalar case has zero additional relaxation cost,
  with a distinct status and no invented crossing witness.
- Thirteen guard tests extend the core/weight suite from 48 to 61 tests.
  Eight theory, six comparison-policy and four routing tests remain separate.
- The forest recurrence and stored timing records are unchanged. Its
  floating-point cost accumulation is not an exact-arithmetic optimality proof.
- The independent diffusion verifier explicitly locates its sibling module,
  including in deeply nested Windows extraction directories.

Original CFD vectors, diffusion data, source records and scientific figures are
unchanged. Existing timing records describe the original benchmark runs, not
a new runtime comparison of this correction. The immutable 2.1.8 record
remains unchanged.
