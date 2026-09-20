# Sparse linear solver: build-vs-buy, and what's implemented so far

## Decision: don't hand-roll sparse factorization -- link a library

Writing a production-grade sparse direct solver for complex-symmetric indefinite
systems (numerically stable pivoting, à la Bunch-Kaufman, plus fill-reducing
reordering) is a large undertaking on its own -- comparable in scope to what the
MUMPS, SuiteSparse, and PARDISO teams have each spent years building. For a
startup with limited runway, that effort is better spent on the parts of this
solver that are actually differentiated (the A-Phi formulation, gauge handling,
EDA/scattering-specific modeling) than on re-deriving general sparse linear
algebra that mature libraries already do well.

**Recommendation: MUMPS**, pending your own evaluation once real problems exist to
benchmark against.

- MUMPS supports complex-symmetric indefinite systems directly (the structure
  this project's A-Phi matrices actually have -- see `docs/CONDITIONING.md`),
  with scaling and fill-reducing reordering (AMD, and METIS if linked in)
  built in.
- MUMPS also supports a general (non-symmetric) matrix type directly, as one
  of three modes alongside symmetric-positive-definite and
  symmetric-indefinite (confirmed via MFEM's `MUMPSSolver` wrapper
  documentation, docs.mfem.org, Sept 2026). This matters because
  `docs/CONDITIONING.md`'s "Formulation 3" (skip the frequency-dependent
  symmetrizing scaling and solve the naturally non-symmetric A-Phi system
  directly) doesn't require a different solver library -- just a different
  mode of the same one.
- Licensing: MUMPS is distributed under **CeCILL-C**, a permissive/weak-copyleft
  license compatible with linking into closed-source commercial software (its
  bundled AMD ordering is separately BSD-3-clause). See the license page at
  mumps-solver.org.
- By contrast, **SuiteSparse's CHOLMOD** (specifically the supernodal module,
  which is what you'd actually want for performance) is **GPL** -- that would
  force a choice between open-sourcing this solver or real legal exposure
  distributing it closed-source. There's a documented case of this catching
  robotics companies shipping Ceres+SuiteSparse without realizing it
  (ceres-solver/ceres-solver#1026 on GitHub). Ruled out for a closed-source
  commercial product on that basis alone, independent of technical merit.
- Intel MKL PARDISO is a reasonable third option technically, but its
  commercial licensing terms have changed hands (Panua Technologies now
  maintains it) -- verify current terms directly before relying on it, don't
  assume last year's terms still hold.

This is a decision to revisit once Phase 04 produces real assembled systems to
benchmark MUMPS against on your actual meshes -- treat "MUMPS" above as the
leading candidate, not a final commitment made in the abstract.

## Planned: solver-mode selection in the input file (design, Sept 2026)

Once Phase 04+ defines a real input-file format, it should expose two
independent choices -- `gauge: albanese_rubinacci | munteanu_unsymmetric`
(`docs/TREE_COTREE_GAUGE.md`) and `frequency_scaling: none | row_scaling |
scaled_phi` (`docs/CONDITIONING.md`'s Formulations 3/1/2) -- rather than
hard-coding either. Both are genuine, measured trade-offs (conditioning vs.
sparsity/fill-in for the gauge; DC-degeneracy vs. memory/factorization cost for
the scaling), not settled questions, and this project's whole practice has
been to let the user pick from real numbers rather than assume a default.

This is *not* a third independent "symmetric vs. non-symmetric solve" flag,
though. Whether the reduced system actually ends up symmetric is fully
determined -- by proof, not measurement -- by the (gauge, frequency_scaling)
combination -- see `docs/CONDITIONING.md`'s decision matrix -- and one
combination (Munteanu-unsymmetric gauge + scaling on) looks like it should
give a symmetric system but doesn't. So the solver-mode (MUMPS
symmetric-indefinite vs. general) should be *derived* from the combination
actually selected via a cheap static lookup on those two settings, not taken
as a third free-standing input-file field a user could set inconsistently
with the other two, and not re-verified numerically on every solve -- the
numerical `||A - A^T|| / ||A||` check belongs in the test suite, as a
one-time-per-code-change guard against an implementation bug breaking the
proven invariant, not on the solve path (`docs/CONDITIONING.md`,
"Interaction with the tree-cotree gauge choice").

## What's implemented now vs. what's waiting on real assembly

Two preprocessing concerns for factorization, and they're handled separately in
this codebase because they solve different problems (see `docs/CONDITIONING.md`
for a fuller discussion of the distinction):

- **Physics-motivated conditioning** (`conditioning.hpp`) -- symmetric row
  scaling, the Phi = j*omega*Phi' substitution. Implemented and tested.
- **Numerical diagonal equilibration** (`equilibration.hpp`) -- symmetric
  Ruiz-style scaling (`diag(d) * A * diag(d)`, kept symmetric on purpose so a
  complex-symmetric A stays complex-symmetric) to bring row/column magnitudes
  into a numerically reasonable range before factorization. Implemented and
  tested, including a check that it composes correctly on top of the physics
  conditioning layer (`tests/test_equilibration.cpp`,
  `test_composes_with_physics_conditioning`) and that it measurably improves the
  estimated condition number on a deliberately badly-scaled test matrix.

Both of the above operate on small hand-built dense matrices right now, the same
way `conditioning.hpp` did when it was first written -- there's no real mesh or
assembly to test against yet.

**Deliberately not implemented yet, and why:**

- **Sparse matrix storage** (CSR/COO) -- straightforward to write, but there's no
  point committing to a specific layout before knowing what MUMPS (or whatever
  solver is chosen) actually wants as input.
- **Fill-reducing reordering** (AMD/METIS) -- if MUMPS is used, it does this
  internally, so writing our own would likely be wasted effort. Even if we
  eventually hand-roll it, reordering only means something relative to a *real*
  mesh connectivity graph (locality, element adjacency) -- testing it against a
  synthetic matrix would validate that the code runs, not that it reduces fill-in
  the way it's supposed to. This waits for Phase 02/03's real mesh graph.
- **The factorization itself** -- waiting on the MUMPS integration (or
  equivalent), which in turn waits on Phase 04 producing a real system to feed
  it.
