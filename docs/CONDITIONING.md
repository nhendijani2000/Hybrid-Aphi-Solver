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
DC singularity of Formulation 1 -- but the recovered `Phi` is now itself scaled by
`omega`, which can reintroduce numerical trouble in a different place (very small
`Phi'` values at very low frequency, if `Phi` itself doesn't vanish as `omega -> 0`).

## Choosing between them

There is no built-in default crossover frequency in this code, and there
shouldn't be one assumed from outside your own measurements. Use the Phase 05
frequency sweep from the project roadmap: for your actual mesh and materials, plot
`estimate_condition_number` (or, at real problem sizes, iterative-solver iteration
count) against frequency for both formulations, and read off where the curves
cross. `recommend_strategy(frequency_hz, crossover_hz)` takes that measured value
as an explicit argument -- it does not embed one.

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
