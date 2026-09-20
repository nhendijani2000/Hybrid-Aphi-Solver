# Scope addition: terahertz phased arrays as a primary target

## What changed (Sept 2026)

Following an independent technical/strategic review (conducted with Google
Gemini; saved alongside this repo in `Gemini/` —
`APhi_Solver_Proposal_Gemini.pdf` and `THz_Phased_Array_APhi_Evaluation.md`),
the project's target application set is sharpened to explicitly name
**terahertz (0.1–10 THz) phased array antennas** — specifically the multiscale
problem of sub-micron semiconductor feeds (photoconductive antenna gaps, HEMT
gates, via transitions) integrated with electrically large radiating
apertures — as a primary commercial beachhead alongside the general EDA
(signal/power-integrity) track from Phase 01.

This is a sharpening of scope, not a change of formulation: every governing
equation, gauge, and solver decision already in `docs/ROADMAP.md` carries over
unchanged. What's new is (a) two additional solver capabilities specific to
phased arrays, (b) an explicit acceleration-method choice for Phase 12, and
(c) a more specific commercial framing for Phase 11.

## Why terahertz is a natural fit for A-Φ, not a stretch

The same multiscale argument that motivates A-Φ over E-/H-field formulations
for EDA applies with more force at terahertz frequencies:

- At 1 THz, free-space wavelength λ₀ = 300 μm; a realistic array element pitch
  is ~λ₀/2 = 150 μm, so a 16×16 array spans roughly 8λ₀×8λ₀ — electrically
  large.
- Meanwhile the feed structures — photoconductive antenna gaps, HEMT gates,
  graphene junctions, via transitions — have critical dimensions from ~20 nm
  to ~500 nm.
- That's a spatial-scale ratio exceeding 10⁴:1 within a single mesh, which is
  exactly the regime where a single-field E-/H-formulation's edge elements
  can't separate the fast-varying capacitive (irrotational) behavior on the
  fine features from the inductive (solenoidal) behavior at the array scale —
  the same multiscale-breakdown mechanism already documented in Phase 01, just
  at a more extreme ratio than a typical EDA interconnect problem.

**Trade-off to keep honest:** A-Φ carries roughly 20–30% more DOFs than a pure
E-field solver on the same mesh (nodal Φ unknowns in addition to edge **A**
unknowns) — a real cost, not a free lunch. It buys numerical stability across
that 10⁴:1 scale ratio, which a same-size E-field system likely can't achieve
at all regardless of DOF count. Treat this percentage, and the illustrative
DOF-count comparison in the Gemini review (~15–25M vs. ~3–5M DOFs on a
hypothetical 16×16 array), as order-of-magnitude estimates, not benchmarked
numbers, until Phases 08/10 produce a real measurement on an actual meshed
array.

## New capabilities this adds (not a reformulation)

1. **Floquet periodic boundary conditions** on a single unit cell (infinite-
   array approximation) — added to Phase 07. This is the standard first
   modeling step for any phased array (scan blindness, central-element active
   impedance) and is far cheaper than the full finite-array FEM-BI hybrid, so
   it's sequenced *before* Phase 12's boundary-integral work as a fast way to
   validate the A-Φ core against COMSOL on a real array problem.
2. **Multi-port active impedance / active reflection coefficient extraction**,
   Γ_active,m(θ,φ) — added to Phase 06. Ordinary S-parameters (already in
   scope) are a two-port-at-a-time concept; a phased array needs the impedance
   seen by each element *as a function of scan angle*, with every other
   element simultaneously excited at its own scanned phase — a distinct
   post-processing step on top of the existing port machinery, not a
   replacement for it.
3. **Gauge validation on multi-port array topologies** — added to Phase 03.
   The existing mixed conductor/dielectric single-port benchmark generalizes
   to *many* ports sharing a common ground/substrate plane. Nothing suggests
   the two candidate gauges (Lagrange-multiplier Coulomb; Chew's generalized
   Lorenz) behave differently with more ports, but that should be checked on a
   realistic array-feed geometry before either gauge is treated as validated
   for this application.

## Locking in ACA over MLFMA for Phase 12

Phase 12 already named "fast multipole method or adaptive cross approximation"
as options for accelerating the dense boundary-integral blocks. For this
application specifically, **ACA is the recommended default**, not just one of
two equally-weighted options:

- MLFMA is asymptotically faster (O(N log N) vs. ACA's typically higher but
  still sub-quadratic cost) but requires an analytic multipole expansion of
  the Helmholtz kernel, and that expansion is well known to degrade at
  sub-wavelength scales — the spherical Hankel functions involved become
  ill-conditioned as their argument shrinks, the standard "low-frequency
  breakdown of MLFMA" — which is exactly the regime a full THz array with
  sub-micron feed detail lives in.
- ACA is purely algebraic (a low-rank cross-approximation computed directly on
  kernel-evaluated matrix blocks) and kernel-independent — it has no analytic
  low-frequency failure mode of its own, consistent with the all-frequency-
  stable design goal running through the rest of this solver.
- This isn't a new idea introduced here — it's the same all-frequency-
  stability reasoning that already drove the A-Φ-over-E/H and the
  potential-BEM-over-EFIE decisions earlier in this roadmap, now made explicit
  for the acceleration method too.

## Commercial framing (Phase 11)

The beachhead target sharpens from generic "EDA and scattering" to
specifically: **multiscale sub-THz/THz systems where sub-micron semiconductor
feeds must be co-designed with radiating apertures** — photoconductive THz
antenna arrays, on-chip antenna-circuit co-design, and automotive radar
packages are the named example verticals. This is a positioning refinement for
Phase 11, not a change to what gets built first; Phases 00–10 are unchanged.

## References for this addition

No new literature enters the project here — every technical claim above
traces back to papers already in `docs/REFERENCES.md` (Zhao & Fu 2017; Yan
2021; Sharma & Triverio 2021, 2022; Lee, Lee & Lee 2003; Ansari, Farquharson &
MacLachlan 2017; Chew, generalized Lorenz gauge). The Gemini review documents
in `Gemini/` are a secondary synthesis of this same public literature, not an
independent source, and were checked against it before folding into this
roadmap.
