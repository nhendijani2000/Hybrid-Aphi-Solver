# Matrix conditioning for the A-Phi system

`include/aphi_solver/conditioning.hpp` implements two algebraic conditioning
transforms for the coupled frequency-domain A-Phi block system

```
[ K_AA    K_APhi  ] [ a   ]   [ rhs_A   ]
[ K_PhiA  K_PhiPhi] [ Phi ] = [ rhs_Phi ]
```

Both are standard results, re-derived here (not lifted from any single source), and
both leave the physical solution unchanged -- the unit tests in
`tests/test_conditioning.cpp` verify that directly by solving a small system both
ways and checking the recovered `(a, Phi)` match.

## Formulation 1 -- symmetric row scaling

Divide the entire scalar-potential row by `j*omega`:

```
row_Phi_new = row_Phi_old / (j*omega)
```

Why this helps: in the natural assembly, `K_APhi` and `K_PhiA` are built from the
same underlying integral (`<grad Phi, A'>`), so up to their scalar coefficients
they are transposes of one another. But the coefficient in front of `K_APhi` is
`(sigma + j*omega*eps)` while the coefficient in front of `K_PhiA` is
`(j*omega*sigma - omega^2*eps)` -- different, so the *assembled* off-diagonal
blocks are not transposes and the system isn't symmetric. Dividing the Phi row by
`j*omega` turns `(j*omega*sigma - omega^2*eps)` into `(sigma + j*omega*eps)`,
matching the other block exactly. `apply_symmetric_row_scaling` implements this;
`test_symmetric_row_scaling_arithmetic` checks the two blocks become transposes on
a small worked example.

Trade-off: dividing by `j*omega` is undefined at DC and increasingly ill-behaved
as `omega -> 0` (that's a division by a shrinking number). `apply_symmetric_row_scaling`
throws if called with `omega == 0`.

## Formulation 2 -- scaled scalar potential

Substitute `Phi = j*omega*Phi'` throughout the system. Wherever `Phi` appeared,
`j*omega*Phi'` now does, so the coupling block `K_APhi` and the scalar block
`K_PhiPhi` each pick up a factor of `j*omega`, while `K_AA`, `K_PhiA`, and both
right-hand sides are untouched:

```
K_APhi_new   = j*omega * K_APhi
K_PhiPhi_new = j*omega * K_PhiPhi
```

Solve for `Phi'`, then recover the physical potential with
`recover_scaled_scalar_potential(phi_prime, omega) = j*omega * phi_prime`.

Trade-off: this rescales columns instead of dividing a row, so it doesn't have the
DC singularity of Formulation 1 -- but it degenerates at low frequency in its own
way, and the two directions are worth keeping straight.

**Corrected 25 Sept 2026.** This paragraph previously said the risk was "very
small `Phi'` values at very low frequency, if `Phi` itself doesn't vanish as
`omega -> 0`". That has the direction backwards. Since `Phi = j*omega*Phi'`,

```
Phi' = Phi / (j*omega)
```

so if `Phi` stays O(1) as `omega -> 0` -- which it does, a terminal held at 1 V
is 1 V at any frequency -- then `Phi'` **diverges**. What becomes small is the
matrix: the `Phi'` diagonal block is `j*omega*K_PhiPhi`, which vanishes, and a
vanishing block against an unchanged right-hand side is precisely what produces
a large `Phi'`. Written in terms of the element integrals (see
`docs/ASSEMBLY_PLAN.md` Sec. 2), the `Phi'` block is `alpha * L` with
`alpha = j*omega*sigma - omega^2*eps`, which is O(omega).

So the low-frequency failure modes are opposite, not shared:

| | `Phi` block scales as | as `omega -> 0` |
|---|---|---|
| Formulation 1 | `beta/(j*omega) * L`, i.e. O(sigma/omega) | blows up |
| Formulation 2 | `alpha * L`, i.e. O(omega) | vanishes |
| Formulation 3 | `beta * L`, i.e. O(sigma) | stays put |

Both 1 and 2 are ill-conditioned near DC relative to the A block (which is O(1)
through its curl-curl part); only the direction differs. That is what the
"Formulations 1-2's DC-degeneracy trade-off" below refers to, and it is
unaffected by this correction.

## Formulation 3 -- solve the natural non-symmetric system directly

Skip both of the above and solve the coupled A-Phi system exactly as assembled,
with `K_APhi != K_PhiA^T`. This has no `j*omega` anywhere in the conditioning
step, so it has none of Formulations 1-2's DC-degeneracy trade-off -- there's no
extra transform that needs special-casing as `omega -> 0`.

Trade-off: you need a general (non-symmetric) sparse solver rather than a
symmetric-indefinite one. `docs/LINEAR_SOLVER.md` confirms MUMPS supports this
directly (a general-unsymmetric matrix type, alongside its symmetric modes), so
this doesn't cost you the solver choice already made -- but a general LU
factorization typically costs roughly 2x the memory and factorization time of a
symmetric-indefinite factorization on the same matrix (full partial pivoting
over the whole matrix instead of exploiting one triangle), and it forecloses
COCG/COCR-type short-recurrence iterative solvers if this project ever goes
iterative instead of direct (those require complex-symmetric structure; a
non-symmetric system would need GMRES or BiCGStab instead).

## Choosing between them

Decide by measurement, the same way as always in this project -- there is no
built-in default crossover frequency in this code, and there shouldn't be one
assumed from outside your own measurements. Use the Phase 05 frequency sweep
from the project roadmap: for your actual mesh and materials, plot
`estimate_condition_number` (or, at real problem sizes, iterative-solver
iteration count and factorization memory) against frequency for all three
formulations, and read off where the curves cross. `recommend_strategy` takes
measured values as explicit arguments -- it does not embed one. Formulation 3
adds a genuine axis to that comparison, not just a fallback: it may be the
right choice even away from `omega -> 0` if factorization memory at your
target mesh sizes turns out to be the binding constraint, per
`docs/ENGINEERING_STANDARDS.md`'s "speed first, then memory, but it's a
trade-off."

## Interaction with the tree-cotree gauge choice (Sept 2026)

The frequency-scaling choice above and the tree-cotree gauge choice
(`docs/TREE_COTREE_GAUGE.md`, Albanese-Rubinacci vs. Munteanu unsymmetric) both
affect whether the final reduced system is symmetric, but they act at
different points in the pipeline, and they do not combine independently:

- **Albanese-Rubinacci (Method A)** eliminates tree-edge DOFs by pure
  row/column restriction (deleting rows and columns of the assembled system).
  Restriction never introduces asymmetry that wasn't already there, so Method A
  is symmetry-neutral: whatever symmetry Formulation 1/2/3 gave the coupled
  system going in, Method A's reduction preserves coming out.
- **Munteanu unsymmetric (Method D)** eliminates tree-edge DOFs by an oblique
  (Petrov-Galerkin) projection -- a different subspace selects rows than
  substitutes columns. `docs/TREE_COTREE_GAUGE.md` Sec. 5 proves this makes the
  reduced matrix non-symmetric even when the input is symmetric. This is not a
  side effect; it is what "unsymmetric" in the method's own name refers to.

That gives four combinations, only three of which are meaningful:

| Gauge | Frequency scaling | Result |
|---|---|---|
| Albanese-Rubinacci | Formulation 1 or 2 (on) | Fully symmetric reduced system -- the only combination that actually achieves this. |
| Albanese-Rubinacci | Formulation 3 (off) | Non-symmetric (from the coupling-block mismatch alone); simplest code, no DC-fragile transform anywhere. |
| Munteanu unsymmetric | Formulation 3 (off) | Non-symmetric (from the coupling-block mismatch *and* the gauge's own projection); best-conditioned per `tests/test_gauge_variants.cpp`'s kappa_D < kappa_A results, at the memory/factorization cost of Formulation 3 plus Method D's own fill-in (`M_ct*F^T` is generally denser than `M`'s own blocks). |
| Munteanu unsymmetric | Formulation 1 or 2 (on) | **Still non-symmetric.** The scaling only fixes the coupling-block mismatch; Method D's projection reintroduces asymmetry regardless. This combination pays Formulation 1/2's DC-degeneracy cost for zero symmetry benefit -- there is no reason to select it. |

**Solver-mode dispatch: a static lookup, not a per-solve numerical check
(revised Sept 2026).** Whether a given `(gauge, frequency_scaling)`
combination yields a symmetric matrix is not uncertain -- it is proven above,
algebraically, for all four combinations. So the solver front-end should
decide symmetric-vs-general dispatch with a cheap O(1) lookup on those two
enum values (the table above, encoded directly), not by computing
`||A - A^T||` from the actual assembled entries on every solve -- that would
be paying an O(nnz) cost (still cheap in absolute terms, but needless) to
re-derive something already established by proof.

What a numerical `||A - A^T|| / ||A||` check is actually useful for is
different: not verifying the *math* (settled above) but catching a future
*implementation bug* that silently violates it -- an edit to the assembly or
gauge code that gets a block wrong, breaking the proven invariant in practice
even though the static lookup still claims it holds. That risk is real enough
in a commercial numerical tool to guard against, but the right place for the
guard is the **test suite**, not the solve path: assert, once per code
change (e.g. in `tests/test_conditioning.cpp`), that each of the three valid
combinations actually produces a matrix with the symmetry this document
claims for it. That costs nothing at solve time, ever, while still catching
a regression the moment it's introduced -- rather than paying even a cheap
check on every production solve for something the proof already guarantees.

Why not just hard-code a threshold: a specific number like "1 Hz" is a property of
a particular mesh, material set, and solver -- not a universal constant of the
A-Phi formulation. Measuring it for your own problems is more correct, and it's
also the only way this project stays built entirely on public derivations and your
own results rather than on someone else's internal validation work.

## `estimate_condition_number`

A small dependency-free utility (power iteration + inverse power iteration on
`A^H A`) for exactly this kind of sweep on small-to-moderate dense systems. It is
not intended to scale to production mesh sizes -- once real meshes are in play
(post Phase 04), the iterative solver's own iteration count is the more practical
conditioning proxy, per Phase 05 of the roadmap.

## Relationship to diagonal equilibration

The two transforms above are physics-motivated: they change *which* algebraic
form of the A-Phi system you're solving. There's a separate, purely numerical
preprocessing step -- diagonal equilibration, to bring matrix entries into a
sane numerical range before factorization -- implemented in
`equilibration.hpp`/`.cpp`. The two compose (equilibrate whichever physics form
you've chosen); see `docs/LINEAR_SOLVER.md` for that module and for the sparse
solver build-vs-buy decision it feeds into.
