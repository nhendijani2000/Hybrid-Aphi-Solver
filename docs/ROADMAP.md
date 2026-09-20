# A-Phi Solver Roadmap

*From formulation to a COMSOL-class full-wave solver*

A phased technical roadmap for an all-frequency-stable A-Φ (magnetic vector
potential / electric scalar potential) finite-element solver — valid from DC
through full-wave, not just low-frequency eddy-current problems — gauged first
by tree-cotree splitting and later also by Coulomb gauge, aimed at EDA
(signal/power-integrity) and scattering applications across that full frequency
range. An exact open boundary via a FEM–boundary-integral hybrid is planned as a
deferred, later phase once the core solver works end to end.

> **Updated (Sept 2026):** the target application set now explicitly names
> **terahertz phased array antennas** — sub-micron semiconductor feeds
> (photoconductive antenna gaps, HEMT gates, via transitions) integrated with
> electrically large radiating apertures — as a primary beachhead alongside
> general EDA. This sharpens scope; it changes nothing about the formulation,
> gauges, or solver strategy already below. See `docs/THZ_PHASED_ARRAY_SCOPE.md`
> for the full rationale, the new capabilities it adds (Floquet boundary
> conditions, multi-port active impedance, an explicit ACA-over-MLFMA choice
> for Phase 12), and an honest caveat on which numbers in it are illustrative
> vs. benchmarked.

| | |
|---|---|
| **Owner** | Nastaran Hendijani |
| **Language** | C++17, CMake |
| **Repo** | `APhi_Solver_Project_LowFrequency_EDA/APhi_Solver` |
| **Phases** | 14 (00–12 plus 03.5, last deferred) |
| **Engineering standards** | `docs/ENGINEERING_STANDARDS.md` (speed-first, then memory; parallelize where it genuinely helps; advanced sparse-matrix techniques from Phase 04 onward) |

Track legend used throughout: **Setup/Core** (foundational, not gauge- or
application-specific) · **Tree-cotree** · **Coulomb gauge** · **EDA** ·
**Scattering** · **Full-wave / FEM-BI**.

> **⚠ Formulation source boundary.** This roadmap and its equations are built
> only from the published papers in `APhi_Papers/` (Lee 2003, Munteanu,
> Li–Sun–Dai–Chew 2015, Yan 2021, and related). The
> `ProposalForA-Phi_Solver_Formulation.pdf` in that folder explicitly condenses
> separate "supplied technical material" — including a specific
> matrix-conditioning threshold and an implementation-history note that don't
> appear in the cited public papers. That reads as your former employer's
> internal know-how rather than published theory, so it is deliberately left
> out here. Keep this boundary in mind as the project grows: published papers
> and your own from-scratch derivations are safe to build on; anything
> condensed from internal slides, reports, or notes is not, regardless of who
> typed it up.

## Contents

00. [Environment, repository, and the IP boundary](#00--environment-repository-and-the-ip-boundary) — Setup
01. [Formulation lock-in and application scope](#01--formulation-lock-in-and-application-scope) — Core, Full-wave
02. [Mesh ingestion and basis functions](#02--mesh-ingestion-and-basis-functions) — Core
03. [Tree-cotree gauge implementation](#03--tree-cotree-gauge-implementation) — Tree-cotree
03.5. [Numerical foundations for assembly](#035--numerical-foundations-for-assembly) — Core
04. [Weak form and matrix assembly](#04--weak-form-and-matrix-assembly) — Core
05. [Linear solver and conditioning](#05--linear-solver-and-low-frequency-conditioning) — Core
06. [EDA track — sources, ports, conductors](#06--eda-track--sources-ports-conductors) — EDA
07. [Scattering track — excitation and open boundary](#07--scattering-track--excitation-and-open-boundary) — Scattering, Full-wave
08. [Verification, validation, and a COMSOL cross-check](#08--verification-validation-and-a-comsol-cross-check) — Core
09. [Coulomb gauge track](#09--coulomb-gauge-track) — Coulomb gauge
10. [Performance and scaling](#10--performance-and-scaling) — Core
11. [Productization — competing with COMSOL](#11--productization--competing-with-comsol) — Core
12. [Full-wave open boundary — FEM-BI hybrid (deferred)](#12--full-wave-open-boundary--fem-bi-hybrid-deferred) — Full-wave, Scattering

---

## 00 · Environment, repository, and the IP boundary

**Track:** Setup

**Goal:** A working C++ toolchain, a git repository with a clean skeleton, and
an explicit, written boundary on what source material this project may draw on.

1. **Repo skeleton — already created.** `APhi_Solver/` now exists next to your
   papers folder with `CMakeLists.txt`, `src/`, `include/aphi_solver/`,
   `tests/`, `docs/REFERENCES.md`, and a `.gitignore` tuned for Visual
   Studio/CMake builds.
2. **Install Visual Studio Community (currently 2026).** Go to
   `visualstudio.microsoft.com/downloads` directly — **not** the "Downloads
   for Visual Studio Subscribers" / older-versions page, which is a
   paid-subscription (MSDN) portal and isn't needed. The free **Community**
   edition is a "Free download" button right on the main downloads page, no
   sign-in wall. In the installer, select the C++ desktop-development
   workload (named "Desktop development with C++" in older versions; check
   the 2026 installer's exact wording) — this pulls in the MSVC compiler,
   the Windows SDK, and CMake tools for Visual Studio; confirm CMake tools
   are checked before installing.
3. **GitHub Copilot is now built into Visual Studio 2026** rather than a
   separate Extensions-marketplace install — look for the Copilot Chat
   window in the IDE and sign in with your GitHub account from there
   (create one at github.com if you don't have one — a free account is
   enough to start). If you end up on an older Visual Studio version
   instead, the original path still applies: `Extensions → Manage
   Extensions → search "GitHub Copilot" → Download`, restart, then sign in.
4. **Open the project.** `File → Open → CMake...` and point it at
   `APhi_Solver/CMakeLists.txt`. Visual Studio configures the CMake cache
   automatically; `Build → Build All` should succeed on the placeholder
   scaffold.
5. **Turn the folder into a git repository.** In a terminal inside
   `APhi_Solver/` (Visual Studio's `View → Terminal` works): `git init`, then
   `git add .`, then `git commit -m "Initial solver scaffold"`.
6. **Create the GitHub remote and push.** On github.com, create a new
   (private, to start) repository named e.g. `aphi-solver` — leave it empty,
   no README. Then locally: `git remote add origin <the URL GitHub gives
   you>`, `git branch -M main`, `git push -u origin main`.
7. **Write down the IP boundary once, in `docs/REFERENCES.md`.** Already
   started: only peer-reviewed papers and your own derivations go into the
   design. Add a line every time a new paper enters the project.
8. **Decide (don't implement yet): license and company structure.**
   Closed-source is the default assumption in this roadmap; revisit at Phase
   11 once there's something to protect.

**Ready for 01 when:** Visual Studio builds the scaffold, `git log` shows a
commit, and the repo exists on GitHub.

---

## 01 · Formulation lock-in and application scope

**Track:** Core, Full-wave

**Goal:** One precise statement of the governing equations you're implementing
first, and which unknowns/geometry classes the EDA and scattering tracks each
require. **Updated:** the project's target shifted from low-frequency-only to
an all-frequency-stable formulation, so the equations below are the full
form — not the reduced low-frequency system this phase originally specified.

> **Written up in full (Sept 2026):** `docs/FORMULATION.md` now carries this
> phase's complete deliverable — the equations below re-derived from
> Maxwell's equations rather than only asserted, the DC special case, the
> reduced low-frequency variant, the DOF layout (including the gauge-reduction
> forward pointer to Phase 03), and both target test problems (a coaxial
> via/power-plane analog for EDA, a PEC sphere Mie benchmark at k·a = 1 for
> scattering). Phase 01's "Ready for 02" checklist is satisfied by that
> document.

```
∇×(1/μ ∇×A) + (jωσ−ω²ε)A + (σ+jωε)∇Φ = J_imp
−∇·[(jωσ−ω²ε)A + (σ+jωε)∇Φ] = 0
```

1. **Keep the full form — retain the ε (displacement-current) terms.** Per
   Zhao & Fu (2017) and Yan (2021), both already in your literature folder,
   this same A-Φ potential system is stable from DC through
   full-wave/radiating regimes when the ε terms aren't dropped — that's the
   whole basis for staying with A-Φ rather than switching to an E-/H-field
   full-wave formulation. The gauge and DOF work in Phases 02–05 carries over
   unchanged; only which coefficients you zero out changes.
2. **Keep the low-frequency-reduced system available as a cheaper optional
   path.** For pure eddy-current/EDA problems where full displacement current
   isn't needed, the ε-dropped, Φ-confined-to-Ω꜀ system from the original
   plan is still valid and cheaper (fewer DOFs) — make it a runtime/config
   choice sharing the same operator blocks, not a second codebase.
3. **Note the ω = 0 special case.** Like most frequency-domain potential
   formulations, the frequency-dependent coefficients above vanish exactly at
   DC, so ω = 0 needs a dedicated magnetostatic solve rather than a naive
   limit of the frequency-domain system — decide now whether that's in scope
   for v1.
4. **Separate the boundary-condition question from the formulation
   question.** A simple zero-tangential-field (PEC) truncation wall cannot
   represent radiation to infinity, regardless of which field variables you
   use internally — that's fixed in Phase 07 (a first-order ABC now, exact
   FEM-BI coupling later in Phase 12), not here.
5. **Write out, for your own reference, why A-Φ is better conditioned here
   than E- or H-field formulations** (per Lee/Lee/Lee 2003, §II) — this
   becomes the one-paragraph pitch for the whole project.
6. **Enumerate EDA requirements separately from scattering requirements**: EDA
   needs stranded/massive conductor excitation, lossy dielectrics, port
   boundary conditions, and thin trace/via geometry. Scattering needs
   plane-wave/dipole excitation and a radiation or absorbing boundary — a
   materially different boundary-condition story from EDA.
7. **Pick the first two target problems**, one per track (e.g. a package
   power-plane / via structure for EDA; a radiating dipole or PEC/lossy
   sphere at a real operating frequency for scattering) — these become the
   running test cases for every later phase.
8. **Updated (Sept 2026): note the terahertz multiscale case as a third,
   sharper version of the same scattering/EDA-adjacent problem.** A phased
   array feed (photoconductive gap, HEMT gate, or via transition at tens to
   hundreds of nm) driving an electrically large radiating aperture pushes
   the spatial-scale ratio that motivates A-Φ over E-/H-field formulations
   past 10⁴:1 — more extreme than a typical EDA interconnect, but the same
   underlying breakdown mechanism from step 5 above. This doesn't add a third
   target problem yet (Phases 06/07 add the array-specific machinery first);
   it's a reason to keep the full-form equations and not defer to a
   quasi-static shortcut. See `docs/THZ_PHASED_ARRAY_SCOPE.md`.

**Grounds this phase:** Zhao & Fu (2017) and Yan (2021) (all-frequency-stable
A-Φ); Lee, Lee & Lee (2003); Dular et al. (2000)

**Ready for 02 when:** the full (all-frequency) equations, the reduced
low-frequency variant kept as an option, the DOF layout, and the two target
test problems (now spanning a real frequency range, not just low-frequency)
are written down in `docs/`.

---

## 02 · Mesh ingestion and basis functions

**Track:** Core

**Goal:** A tetrahedral mesh data structure with the incidence information
(node–edge, edge–facet) that every later phase depends on, plus lowest-order
Whitney basis functions.

1. **Mesh ingestion:** read an unstructured tetrahedral mesh (start with a
   simple format — Gmsh `.msh` is a reasonable first target) into nodes,
   edges, facets, tets.
2. **Build the incidence matrices** the whole gauge story depends on: the
   discrete curl operator **C** (edge–facet) and the discrete
   gradient/divergence operator **G** (node–edge), exactly as defined in the
   Munteanu paper — these are sparse, entry-per-orientation, and cheap to
   build once mesh connectivity exists.
3. **Basis functions — updated Sept 2026 (final), see `docs/FORMULATION.md`
   §5.1 for the full decision record and history.** Two earlier stances were
   superseded (kept for the record in FORMULATION.md): matched first order
   (Whitney **A** + Lagrange P1 Φ), then a brief matched-second-order plan
   for both fields. **Locked decision: mixed order.** Φ uses the verified
   10-node quadratic nodal element from Jin (2014), Ch. 5, Eq. (5.52)/Fig.
   5.3. **A** uses the classic first-order (lowest-order) Whitney/Nédélec
   tetrahedral edge element — 6 DOF per tet, one per edge — verified both
   against Bossavit (1998) and against the user's own confirmed prior
   implementation (`BasisFunction.docx` / `EdModel.cpp`, from her university
   "3dedyaphi" codebase; confirmed not Ansys-derived). This is a **deliberate
   mixed-order pairing**, not an oversight: `grad(P2 Phi)` is not exactly
   representable in the first-order **A** edge space (see the from-scratch
   proof in `docs/FORMULATION.md` §5.1), which caps the A-Phi coupling
   term's accuracy at first order — an accepted, documented trade-off, not a
   stability problem. Tree-cotree gauging is unaffected (Phase 03 below).
   The same source document also describes a separate 10-point-per-tet
   reconstruction of the solved **A** field (for combined plotting with Φ);
   that is a post-processing utility for later phases, not part of this
   phase's discretization.
4. **True higher-order/hierarchical bases and hp-adaptivity**, and any
   future matched-second-order upgrade to **A**, stay deferred to Phase 10
   — get a correct pipeline working end-to-end first. (Lee 2003's
   hierarchical hp scheme was considered as a Phase 10 reference, but note
   it is a 2-D triangular-element scheme, not 3-D tetrahedral — see
   `docs/FORMULATION.md` §5.1 — so it illustrates the general "match nodal
   order to edge order" principle only, not a directly reusable 3-D formula
   set. Graglia, Wilton & Peterson (1997) and García-Castillo et al. (2000)
   remain identified, unread candidates for that future upgrade.)
5. **Done (Sept 2026) — global edge-orientation signs for the **A** basis.**
   A tet's local vertex order comes from the mesh file and is arbitrary,
   while a global edge DOF is defined on the canonical low-index →
   high-index direction. Where the two disagree — measured at **~58% of
   (tet, local edge) pairs** on `meshes/cube_*.msh` — the Whitney function
   built from the local pair is the *negative* of the global edge's basis
   function, so scattering it into a global DOF without a sign silently
   negates more than half of all element contributions and breaks
   tangential (H(curl)) continuity between neighbouring tets. This is the
   `mEdgeSign` convention of the user's own prior 3dedyaphi implementation
   (`docs/FORMULATION.md` §5.1), implemented here as `Mesh::tet_edge_signs`
   (built in `build_topology`, alongside `tet_edges`) and applied by
   `whitney_edge_value_global` / `whitney_edge_curl_global` in
   `basis_functions.hpp` — **which are the forms Phase 04 assembly must
   call**; the local-oriented `whitney_edge_value` / `whitney_edge_curl`
   must not be scattered into a global DOF directly. Φ needs no sign: the
   P2 edge-midpoint shape function `4*L_v0*L_v1` is symmetric in its two
   vertices. Caught late because both single-tet test fixtures used
   vertex order `{0,1,2,3}` (ascending, so every sign is +1, making an
   orientation bug structurally invisible); the regression cover added with
   the fix is a two-tet mesh whose second tet lists vertices non-ascending,
   checking the global circulation identity (δ_jm along each edge's
   canonical direction) and tangential continuity across the shared face.
   Verified to genuinely catch the bug: with the sign neutralized the new
   checks fail on exactly the reversed edges (426/431), with it restored
   431/431.

**Grounds this phase:** Munteanu (tree-cotree condensation properties);
J.-M. Jin (2014), Ch. 5 (verified Φ-side second-order nodal element);
Bossavit (1998), Ch. 5 (verified first-order Whitney/Nédélec **A**); the
user's own prior implementation (`BasisFunction.docx` / `EdModel.cpp`);
Lee, Lee & Lee (2003), §IV–VI (background/motivation only — 2-D scope, see
`docs/FORMULATION.md` §5.1)

**Ready for 03 when:** **C** and **G** can be assembled for a test mesh and
**CG = 0** holds numerically (the discrete curl·grad = 0 identity) — this is
your first real unit test — and the global edge-orientation invariant of
step 5 holds on a mesh with non-ascending tet vertex ordering (circulation
δ_jm along each edge's canonical direction; tangential continuity across a
shared face). Both are now covered by `tests/`.

---

## 03 · Tree-cotree gauge implementation

**Track:** Tree-cotree

**Goal:** A working spanning-tree gauge so the assembled system is
non-singular, plus a deliberate choice among the known tree-cotree variants
rather than an arbitrary one. **Updated:** also decide, on real evidence
rather than assumption, what to do about ports that straddle a conductor and
a dielectric — a case tree-cotree doesn't handle cleanly, and where the
generalized (penalty) Coulomb gauge alone is not guaranteed to be reliable
either.

**Resolved (Sept 2026):** an earlier draft of this note flagged an open
question here — whether tree-cotree gauging applies to the whole element or
only a lowest-order subspace — raised while a matched second-order **A**
was the working plan. Phase 02's final decision kept **A** at plain
first-order (Whitney/Nédélec, mixed with a second-order Φ; see
`docs/FORMULATION.md` §5.1), so that question no longer applies here: the
classical tree-cotree treatment below goes through exactly as originally
written, with no dependency on the still-unread Graglia/Wilton/Peterson
(1997) or García-Castillo et al. (2000) candidates. Those remain of
interest only if a future matched-second-order **A** is revisited (Phase
10).

1. **Done (Sept 2026).** Implemented the spanning-tree search over the mesh
   graph, grounded in S.-C. Lee, J.-F. Lee, R. Lee (2003) Sec. V (fetched and
   checked against the actual paper before writing any code): node numbering
   with PEC/ground handling, then tree/cotree edge marking, fused into a
   single BFS pass (`include/aphi_solver/tree_cotree.hpp`,
   `src/tree_cotree.cpp`; 73/73 checks in `tests/test_tree_cotree.cpp`). A
   CSR `NodeAdjacency` and a `UnionFind` (path compression + union by rank)
   were used for this rather than the `std::map`-based lookups elsewhere in
   the codebase, per `docs/ENGINEERING_STANDARDS.md`.
2. **Done (Sept 2026).** Implemented two gauging variants from Munteanu's
   projection framework (fetched and verified against the actual paper --
   `docs/REFERENCES.md` has the full verification record and one disclosed
   caveat about a symbol collision in the OCR'd equations 19-20, resolved by
   deriving the reduced-matrix formula from her general framework instead of
   copying that specific equation text): **Albanese–Rubinacci** (simplest,
   keeps sparsity, but worse conditioning -- `build_albanese_rubinacci_gauge`)
   and **Munteanu unsymmetric** (denser, but better-conditioned per her
   numerical tests -- `build_munteanu_unsymmetric_gauge`), in
   `include/aphi_solver/gauge_variants.hpp` / `src/gauge_variants.cpp`. The
   "essential incidence matrix" F the unsymmetric variant needs is computed
   via an O(V) tree back-substitution (not a general dense matrix inversion),
   exploiting that the tree's own node-incidence structure lets each group's
   potential be written as a cumulative sum down from its reference root.
   The full step-by-step derivation of both reduced-matrix formulas (not
   just a summary of what was/wasn't trusted from the paper) is written up
   independently in `docs/TREE_COTREE_GAUGE.md`.
3. **Done (Sept 2026).** Instrumented condition-number reporting via a crude
   power-iteration estimate (`estimate_condition_number`, largest/smallest
   eigenvalue of A^T*A via power/inverse-power iteration -- correct for both
   variants uniformly, since Method D's reduced matrix is not symmetric and
   so its eigenvalues alone would not give a meaningful condition number).
   On both of this project's (tiny) test meshes, the unsymmetric variant
   came out better-conditioned than Albanese-Rubinacci (single tet: kappa_A
   = 4, kappa_D = 1; two tets sharing a face: kappa_A ≈ 13.1, kappa_D ≈ 2.0)
   -- the right *direction* relative to Munteanu's own kappa_D < kappa_A
   ordering, though these meshes are far too small to confirm her full
   ordering or exact values. `tests/test_gauge_variants.cpp` also checks a
   stronger, independent correctness property: for a manufactured,
   automatically-consistent source j = M*z, both gauge variants recover a
   vector potential whose curl matches curl(z) exactly, even though the two
   recovered potentials differ from each other and from z (different valid
   gauges, same physical field) -- 16/16 checks passing.
   **Updated (Sept 2026):** `estimate_condition_number` was rewritten to be
   sparse (matrix-free A^T*A matvecs, Conjugate Gradient inner solve instead
   of dense Gaussian elimination -- see `docs/ENGINEERING_STANDARDS.md` and
   `docs/REFERENCES.md`, "Matrix conditioning"), brought forward from its
   originally-planned Phase 04/05 slot so a gauge-comparison tool can run on
   realistic mesh sizes rather than only Phase 03's tiny validation meshes.
   Verified independently against diagonal test matrices with analytically-
   known condition numbers (identity, diag(1,2,4,100), and a 500x500
   diagonal with kappa = 1000 exactly), not just re-checked against the
   gauge variants -- 19/19 checks now passing in
   `tests/test_gauge_variants.cpp`.
4. **Handle multiple PEC bodies** explicitly (each grounded to its own
   reference node) — a common source of subtle bugs later. **Already done as
   part of step 1**: `build_tree_cotree` discovers PEC-tagged nodes'
   connected components automatically (via UnionFind) and gives each
   distinct body its own reference group rather than merging them, per
   `tests/test_tree_cotree.cpp`'s disconnected-PEC-bodies test case.
5. **Build a mixed conductor/dielectric port test case early, and benchmark
   two gauge alternatives against it.** The standard tree-cotree recipe
   grounds the spanning tree through whole PEC bodies; a port face that is
   only partly conductor doesn't fit that recipe, and per Rapetti, Alonso
   Rodríguez & De los Santos (2022), a tree gauge isn't a discretization of
   any orthogonality condition to begin with, so there's no variational
   fallback right at that boundary. The generalized (penalty) Coulomb gauge
   from Phase 09 isn't automatically safe here either: Ansari, Farquharson &
   MacLachlan (2017) show that an implicitly- or softly-enforced Coulomb
   gauge can produce non-unique A/Φ — different solvers converging to
   different potentials for the same physical field — when the normal
   component of **A** is discontinuous across a material interface, exactly
   the situation at a conductor/dielectric port; Demerdash & Wang (1990)
   documented a related Coulomb-gauge breakdown at high/low-permeability
   material contrasts. On this one synthetic geometry, benchmark: (a) an
   *explicitly*-constrained Coulomb gauge via a Lagrange multiplier that
   carries the interface-jump term directly (Ansari, Farquharson &
   MacLachlan's fix, rather than the softer penalty form already in Phase
   09), and (b) Chew's generalized Lorenz gauge for inhomogeneous media
   (∇·(εA) = −χ∂Φ/∂t, χ = αε²μ), which is built from the start for
   spatially-varying ε/σ and derives the correct conductor/dielectric
   interface conditions as part of the gauge itself rather than patching
   them on afterward. Decide the default only once you have numbers from
   this test case, not from the literature alone.
   **Updated (Sept 2026): implement Chew's generalized Lorenz gauge first.**
   It adds no new unknowns (no Lagrange-multiplier field, no saddle-point
   solver, no inf-sup pairing to get right) and its frequency-dependent
   gauge condition is a more natural fit for a solver whose whole premise is
   full-wave/radiating behavior, whereas Ansari's method comes from the
   galvanic/inductive geophysical literature. This is a prioritization for
   implementation order, not a decision made from the literature alone —
   Ansari's explicit Lagrange-multiplier Coulomb gauge stays in as the
   required benchmark comparison on the same synthetic geometry before
   either one is called the default, exactly as the sentence above already
   requires. Note also that, like the classical Lorenz gauge, Chew's
   generalized version is frequency-dependent and should be expected to
   degenerate at ω = 0 the same way — the dedicated DC/magnetostatic solve
   already planned in Phase 01, step 3 covers that case regardless of which
   AC gauge is chosen here.
6. **Updated (Sept 2026): extend the mixed-port benchmark to a multi-port
   array topology before calling either gauge validated for phased-array
   feeds.** The single-port test above generalizes to *many* ports sharing a
   common ground/substrate plane, which is the actual feed topology for a
   phased array. There's no specific reason in the literature to expect the
   Lagrange-multiplier Coulomb gauge or Chew's generalized Lorenz gauge to
   behave differently with more ports present, but that's an assumption, not
   a result — check it on a small synthetic array (e.g. a 2×2 or 4×4 feed
   grid) rather than assuming the single-port result generalizes for free.
   See `docs/THZ_PHASED_ARRAY_SCOPE.md`.

**Grounds this phase:** Munteanu, §IV–VI (variants A–E, condition-number
ordering κ_D < κ_C < κ_A < κ_E < κ_B); Lee, Lee & Lee (2003), §V; Rapetti,
Alonso Rodríguez & De los Santos (2022), tree gauge vs. Coulomb gauge
distinction; Ansari, Farquharson & MacLachlan (2017), explicit/Lagrange-
multiplier Coulomb gauge for mixed galvanic/inductive problems; Chew
(2014/2022), generalized Lorenz gauge for inhomogeneous anisotropic media;
Demerdash & Wang (1990), Coulomb-gauge failure at material contrasts.

**Ready for 04 when:** the gauged system is non-singular on your two test
meshes, the condition-number ordering you measure is consistent with
Munteanu's reported ordering, and the mixed conductor/dielectric port
benchmark has a documented result and a chosen gauge for that case (not left
as an open question for Phase 06).

---

## 03.5 · Numerical foundations for assembly

**Track:** Setup/Core

**Added Sept 2026, after a review at the 03 → 04 boundary.** Phase 04 was
stalling, and the reason turned out to be structural rather than
mathematical: as originally written it silently required four interlocking
prerequisites that did not exist yet — a complex sparse matrix type, mesh
region/boundary tags, a gauge reduction that can act on part of a larger
block system, and correct global edge orientation — each of which changes
the *signatures* the others are written against. Attempting them inside
Phase 04 meant every attempt at assembly code surfaced another gap, so the
work produced plans instead of matrices. Splitting them out into a phase
with its own exit gate is the fix; the numbering is fractional on purpose,
so the many `Phase 04`/`Phase 05` references in `docs/FORMULATION.md`,
`docs/CONDITIONING.md`, `docs/LINEAR_SOLVER.md` and
`docs/TREE_COTREE_GAUGE.md` stay valid.

**Goal:** Everything Phase 04's assembly loop needs to already exist and
have a settled type signature, so that phase is only about the weak form.

1. **Done (Sept 2026) — global edge-orientation signs.** Completed as Phase
   02, step 5 above; listed here because it is a hard prerequisite for any
   assembly and was the one item of the four that was an outright bug
   rather than a missing piece.
2. **Done (Sept 2026) — one templated sparse matrix type, replacing
   today's split.** Shipped as `Sparse<T>` in
   `include/aphi_solver/sparse_matrix.hpp` (aliases `SparseMatrixD`,
   `SparseMatrixZ`), with `incidence.hpp`'s old COO struct now an alias for
   the real instantiation and `APhiBlockSystem` holding sparse blocks with
   plain `std::vector<Complex>` right-hand sides. `conditioning.hpp` gains
   `assemble_sparse` (the production path) alongside `assemble_dense` (the
   ground truth, deliberately kept independent so their cross-check is not
   tautological); `equilibration.hpp` moved with them rather than forcing a
   dense round-trip mid-solve. The migration deleted more than it added:
   `incidence.cpp`'s hand-written coalesce/multiply/transpose and
   `gauge_variants.cpp`'s local matvec helpers all became members of the one
   type. Original text follows. The
   codebase currently has two incompatible halves: `SparseMatrix`
   (`incidence.hpp`) is real-valued COO with a `std::map`-based coalesce,
   and `ComplexMatrix` (`complex_matrix.hpp`) is complex but **dense** —
   and `APhiBlockSystem` is built from the dense one. Phase 04 needs
   complex *sparse*, which neither provides; a dense `K_AA` is already
   55 MB at `meshes/cube_6.msh` and impossible at any real EDA mesh size.
   Implement a single `Sparse<T>` (aliases `SparseMatrixD = Sparse<double>`
   for the incidence operators **C**/**G**, `SparseMatrixZ =
   Sparse<Complex>` for the system blocks): accumulate as triplets, then
   `compress()` once to CSR. **Keep a triplet export as well** — MUMPS's
   assembled centralized input format wants IRN/JCN/A arrays, not CSR, so
   both shapes are needed regardless (`docs/LINEAR_SOLVER.md`). Migrate
   `APhiBlockSystem` onto it; keep the dense `ComplexMatrix` only as the
   small-system ground truth for `solve_dense`, which Phase 05 wants
   anyway. This is the `docs/ENGINEERING_STANDARDS.md` item 3 requirement
   ("compressed sparse formats … rather than the COO/`std::map`-based
   `SparseMatrix`"), which that document already scheduled for Phase 04 but
   which nothing had been assigned to actually do.
3. **Mesh region and boundary tag ingestion.** `read_gmsh_msh` currently
   discards *every* element tag — including each tet's own physical-group
   id — and skips non-tet elements entirely, so the boundary triangles Gmsh
   uses to mark PEC walls, ports and the outer truncation never reach
   `Mesh` at all. Since DOF setup, material lookup and PEC handling all key
   off "which entity belongs to which named region," this is a hard
   prerequisite, not a later refinement. Add `Mesh::tet_tags` (one
   physical-group tag per tet, parallel to `Mesh::tets`) and a tagged
   boundary-face list capturing Gmsh elm-type 2 triangles with their first
   tag, remapped to 0-based node indices. `Mesh` carries integers only —
   what tag `7` *means* is the input file's job, keeping `gmsh_reader.cpp`
   format-agnostic and dependency-free.
4. **Done (Sept 2026) — generalize the gauge reduction to an explicit A-DOF
   index set.** Shipped as `GaugeIndexMap` +
   `build_albanese_rubinacci_index_map(a_dof_edge, num_phi_dofs, tc)` in
   `gauge_variants.hpp`, applied through `Sparse<T>::principal_submatrix`.
   The split that made this work: the index map holds *only indices* and no
   scalar type, so one map drives the reduction of a real or a complex
   matrix, while the reduction itself is one templated primitive. Method A
   is literally a principal submatrix, which is why it has no fill-in.
   `a_dof_edge[k]` gives the mesh edge that A-DOF `k` represents — the
   identity in the full-wave regime, a genuine mapping once PEC tangential
   edges leave the unknown set. `build_albanese_rubinacci_gauge` now routes
   through the same path with an identity map and zero Φ DOFs, so there is
   one implementation rather than two; its 19 existing checks pass
   unchanged and `compare_gauges` reports bit-identical numbers on
   `cube_4.msh`. `restrict_vector` / `expand_solution` handle the RHS and
   the solution, with eliminated tree entries expanding to exactly zero —
   which is the gauge condition `a_t = 0` itself, not padding. Original
   text follows. Both
   `build_albanese_rubinacci_gauge` and `build_munteanu_unsymmetric_gauge`
   currently assume the matrix's row/column index space *is* the mesh's
   global edge index space, 1:1 — true for their present callers
   (`test_gauge_variants.cpp` and `compare_gauges.cpp`, which only ever
   pass `M = CᵀC`), but wrong for Phase 04, where gauge reduction must act
   on only the A-DOF rows/columns of a larger coupled system while Φ rows
   and columns pass through untouched. Generalize the existing functions to
   take that index set rather than writing block-aware duplicates: the
   reduction math (test space W, trial space V — `docs/TREE_COTREE_GAUGE.md`)
   is identical either way, and duplicating it is exactly what
   `docs/ENGINEERING_STANDARDS.md` warns against. Combined with step 2's
   template, this also resolves the real-vs-complex question in one move
   rather than forcing a second complex-typed copy of the same three
   functions.
5. **Done (Sept 2026) — make Method D sparse.** Found while auditing what
   still densifies after step 2. `build_munteanu_unsymmetric_gauge`
   densified `M` (E x E) and built a dense `L^T`, and
   `compute_essential_incidence_matrix` built `F` dense via one O(V) tree
   solve per free group — O(V^2) time and memory, roughly 80 GB for `M` and
   10 GB for `F` at 10^5 edges, so Method D could not run on a real mesh at
   all. Method A was never affected (it is a pure principal submatrix),
   which is why this did not block the chosen path; but every fill-in
   number `compare_gauges` reported was only obtainable on toy meshes.
   The fix came from a structural observation rather than a data-structure
   swap: `G_t^{-1}[m, col]` is `+-1` exactly when group `col` is an
   ancestor of `m`, so **row `r` of F is exactly the fundamental cycle of
   cotree edge `r`** — a tree path bounded by the spanning tree's diameter,
   with everything above the endpoints' common ancestor cancelling. F is
   now built by walking both endpoints up to their common ancestor, and the
   reduction is a sparse-sparse product. `solve_tree_system` was deleted
   (the cycle walk replaced its purpose). Verified identical, not merely
   passing: `compare_gauges` reports the same kappa and nnz on `cube_4` and
   `cube_6` as the dense implementation, checked by building the previous
   commit side by side. Gauge-build time (excluding the condition-number
   estimate, which dominates wall time at these sizes): cube_4 22.9 -> 17.6
   ms, cube_6 87.3 -> 24.2 ms — the old path grows ~3.8x when the mesh
   grows 3x, the new one ~1.4x.
6. **Clean up two API collisions before they propagate.** `dense_solve`
   (real, `gauge_variants.hpp`) vs. `solve_dense` (complex,
   `complex_matrix.hpp`) differ only in word order, and
   `estimate_condition_number` exists twice with different types in
   different headers. Both are survivable today and actively confusing once
   step 2 makes the scalar type a template parameter.

**Ready for 04 when:** a complex sparse matrix can be built, compressed,
multiplied by a vector and exported as triplets; `APhiBlockSystem` holds
sparse blocks; a tagged `.msh` round-trips its tet and boundary-face tags
through `Mesh`; and **Albanese-Rubinacci** reduces a system in which the
A-DOFs are a strict subset of the index space, with the existing 19/19 gauge
checks still passing unchanged. Generalizing Munteanu unsymmetric the same
way is expected to fall out of the same index-set machinery (step 4) and
should be done alongside, but it does not gate Phase 04 — per Phase 04 step
5, Method A is the path to the first validated solve.

---

## 04 · Weak form and matrix assembly

**Track:** Core

**Goal:** A sparse assembly pipeline producing the coupled A-Φ block system,
with both the tree-cotree gauge choice and the frequency-scaling choice
available as solver-conditioning knobs — see `docs/CONDITIONING.md`'s
"Interaction with the tree-cotree gauge choice" (Sept 2026) for why these two
knobs aren't independent and what the input file should actually expose.

```
K_AA  = ∫ (∇×A)·(∇×A′) dΩ        M_AA  = ∫ A·A′ dΩ
K_AΦ  = ∫ ∇Φ·A′ dΩ                K_ΦΦ = ∫ ∇Φ·∇Φ′ dΩ
```

**Prerequisites:** Phase 03.5 above. This phase assumes a complex sparse
type, mesh tags, index-set-aware gauge reduction, and correct edge
orientation all already exist — it is about the weak form and nothing else.

1. **A `Problem` struct first, its file format second.** Define the
   in-memory problem description (material regions keyed by mesh tag,
   boundary-condition types, gauge choice, frequency-scaling choice,
   operating frequency, source) as a plain C++ struct, and let JSON parsing
   be a thin adapter written on top of it afterwards. The struct is the
   interface DOF setup actually needs; sequencing it first means DOF setup
   is not blocked on vendoring a JSON library. When the parser does land,
   have it cross-check that every tag the input file names exists in the
   loaded mesh — a typo'd region tag should fail loudly, not silently
   assemble a region with zero material properties.
2. **DOF numbering**, per `docs/FORMULATION.md` §5.2 and the mixed-order
   pairing locked in there (first-order Whitney **A**, P2 Φ — reaffirmed
   Sept 2026): `N_A = N_edges`, `N_Φ = N_vertices + N_edges`.
   **Φ lives everywhere in the full-wave regime (confirmed Sept 2026).**
   With displacement current on, Φ is defined over the whole domain —
   dielectric and free-space regions included — not confined to `Ω_c`.
   The `Ω_c` restriction belongs *only* to the reduced low-frequency
   variant (`docs/FORMULATION.md` §3), where `eps -> 0` removes the
   physical content Φ would carry outside the conductors. So the DOF map
   needs `Ω_Φ = whole domain` as its default path and `Ω_Φ = Ω_c` (derived
   from `sigma > 0`, not separately tagged) only when the reduced regime is
   selected — a domain-restriction input, not a branch in the mathematics.
   Derive `is_pec` for `build_tree_cotree` from the Phase 03.5
   boundary-face tags rather than the hand-built `std::vector<bool>`
   `compare_gauges.cpp` uses today. Independently testable ahead of any
   assembly: given a known small mesh and a known `Problem`, assert the
   exact DOF counts and index maps by hand, the way
   `tests/test_tree_cotree.cpp` already does.
3. **Element-level matrices**, tested on a single reference tet against
   hand-computed values, before any global assembly code runs. Two things
   worth fixing up front rather than discovering mid-implementation:
   - **Use `whitney_edge_value_global` / `whitney_edge_curl_global`**
     (Phase 02, step 5), never the local-oriented forms. More than half of
     all element contributions are wrong otherwise.
   - **A degree-2 (4-point) tetrahedral quadrature rule is exact for every
     integrand in this phase** — no guesswork needed. Whitney functions are
     linear in position, so `M_AA` and the `K_AΦ` coupling are quadratic;
     `curl(Whitney)` is constant per tet, so `K_AA`'s curl-curl term is
     exact at a single point; `grad(P2)` is linear, so `K_ΦΦ` is quadratic.
     Region-wise-constant material coefficients add no degree.
4. **Validate at DC before going complex.** At `omega = 0` the system
   decouples exactly (`docs/FORMULATION.md` §2): Φ solves DC conduction
   alone, **A** solves magnetostatics with the tree-cotree gauge. That
   makes a real-valued, uncoupled, analytically checkable first milestone
   (straight wire → `B = mu*I/(2*pi*r)`) which nonetheless exercises the
   whole pipeline end to end — tags → DOF map → oriented element matrices →
   global sparse assembly → gauge reduction → sparse solve. It is also
   already this phase's own "Ready for 05" criterion. Getting a checkable
   number out roughly half the pipeline sooner is the point; turn `omega`
   on only once it is green.
5. **Global sparse assembly** into the Phase 03.5 sparse type, respecting
   the tree-cotree reduction from Phase 03 (reduced-size system, not the
   full singular one). **Albanese-Rubinacci (Method A) is the first target
   (confirmed Sept 2026):** prescribe `a_t = 0` on every tree edge and
   solve the principal submatrix on the cotree DOFs. It is the simplest
   reduction available — pure row/column restriction, a genuine submatrix
   of the assembled system, so it preserves the sparsity pattern with no
   fill-in and is symmetry-neutral (`docs/CONDITIONING.md`, "Interaction
   with the tree-cotree gauge choice"). Get a correct end-to-end solve on
   Method A before wiring the gauge choice as a user-facing knob. Munteanu
   unsymmetric (Method D) stays implemented and tested, and the
   `(gauge, frequency_scaling)` selection in step 6 is still designed for
   both — but Method D's better conditioning is not worth paying for in
   debugging surface while the pipeline has never produced a validated
   field. Per-tet element work is embarrassingly parallel; the scatter into
   the global matrix needs thread-local accumulation plus a merge, or
   atomics, to stay correct (`docs/ENGINEERING_STANDARDS.md`).
6. **Implement the frequency-scaling choice** (natural/non-symmetric,
   row-scaling, or scaled-Φ — `docs/CONDITIONING.md` Formulations 1-3) as a
   build-time or run-time flag — keep all forms available since later phases
   compare them, and derive/validate the solver's symmetric-vs-general mode
   from the actual (gauge, scaling) combination rather than trusting it
   blindly (`docs/LINEAR_SOLVER.md`, "Planned: solver-mode selection").
7. **Right-hand-side assembly** for the source models you'll add properly in
   Phases 06–07; a simple uniform current density is enough here to exercise
   the pipeline.

**Grounds this phase:** general A-Φ weak form, standard in the cited
literature (Dular et al. 2000; Zhao & Fu 2017)

**Ready for 05 when:** the DC milestone (step 4) reproduces the straight-wire
field, the assembled sub-blocks are complex-symmetric where they should be
(K_AA, K_ΦΦ symmetric; K_AΦ and K_ΦA transposes of one another up to their
respective coefficients, per `docs/CONDITIONING.md`), and the same trivial
problem gives the right sign and rough magnitude of **B** at `omega != 0`
as well.

---

## 05 · Linear solver and low-frequency conditioning

**Track:** Core

**Goal:** A solver strategy that stays accurate as frequency → 0, built
around what the assembled matrix actually is: sparse, complex *symmetric*
(A = Aᵀ, not Hermitian — real Whitney/nodal basis functions times complex
coefficients), and typically ill-conditioned. Scaling and reordering happen
*before* factorization, not as an afterthought.

1. **Dense direct-solve baseline — done.** `solve_dense` in the conditioning
   module, ground truth for every sparse/iterative result on small problems.
2. **Diagonal equilibration — done.** Symmetric Ruiz-style scaling
   (`equilibration.hpp`): `diag(d) A diag(d)`, kept symmetric on purpose so a
   complex-symmetric A stays complex-symmetric — a numerical-stability step
   distinct from the physics-motivated conditioning transforms in Phase 04,
   and tested to compose correctly with them (`docs/LINEAR_SOLVER.md`).
3. **Sparse solver: recommendation is MUMPS, not a from-scratch
   factorization.** A robust complex-symmetric sparse direct solver with
   stable pivoting (Bunch–Kaufman-style) is a large, mature-library-sized
   undertaking on its own. MUMPS supports complex-symmetric indefinite
   systems directly with reordering/scaling built in, and is CeCILL-C
   licensed — permissive enough for closed-source commercial use.
   SuiteSparse's CHOLMOD (the module you'd actually want for performance) is
   GPL, which risks forcing this solver open-source if shipped closed-source
   — ruled out on licensing alone. Verify Intel MKL PARDISO's current
   commercial terms directly if considering it (they've changed hands). Full
   writeup: `docs/LINEAR_SOLVER.md`. Treat this as the leading candidate,
   confirmed once real problems exist to benchmark against.
4. **Sparse storage and fill-reducing reordering wait for Phase 04.**
   Committing to a CSR/COO layout before knowing what MUMPS wants as input is
   premature, and reordering (AMD/METIS) only means something relative to a
   *real* mesh connectivity graph — a synthetic test would validate the code
   runs, not that it reduces fill-in the way it's supposed to.
5. **Iterative solver track, in parallel**: COCG or COCR
   (conjugate-gradient-type methods that exploit complex-symmetric structure
   directly, cheaper per iteration than general non-symmetric GMRES) as the
   primary candidate, with GMRES as a fallback if COCG-type methods struggle
   on your specific matrices; basic preconditioner (ILU) to start.
6. **Sweep frequency on a fixed mesh** and plot condition number / iteration
   count vs. ω for each conditioning strategy from Phase 04 — this is where
   you empirically rediscover why low-frequency breakdown happens, what your
   gauge choice buys you, and where your own natural-vs-scaled crossover
   actually falls (see `docs/CONDITIONING.md`).

**Grounds this phase:** Yan (2021), §1–2 (low-frequency breakdown mechanism);
Li, Sun, Dai & Chew (2015), Table I (condition-number / iteration-count
reporting style to imitate); general sparse-direct-methods background (e.g.
Davis, *Direct Methods for Sparse Linear Systems*) and complex-symmetric
iterative methods (COCG: van der Vorst & Melissen 1990) — standard numerical
linear algebra, not specific to any one implementation.

**Ready for 06 when:** you have a frequency-vs-conditioning plot for your two
test problems, MUMPS (or the alternative you land on) actually linked and
factoring a real assembled system, and a default solver configuration
documented.

---

## 06 · EDA track — sources, ports, conductors

**Track:** EDA

**Goal:** The excitation and material machinery specific to
signal/power-integrity problems: stranded and massive conductors, lossy
dielectrics, and a usable port abstraction.

1. **Implement stranded-conductor and massive-conductor source models** as
   two distinct excitation types (they couple differently into the A-Φ
   system) — this is a well-studied distinction, not something to improvise.
2. **Add lossy-dielectric and finite-conductivity material regions**, and
   thin-conductor/via modeling appropriate for interconnect geometry.
3. **Define a port abstraction** (voltage/current or wave-port style)
   sufficient to extract S-parameters or impedance from the field solution —
   this is what makes the solver usable for actual EDA problems rather than
   only field plots.
4. **Run your Phase 01 EDA test problem end-to-end** and sanity-check against
   a hand or textbook estimate (e.g. DC resistance, or a simple RLC estimate
   for a via).
5. **Updated (Sept 2026): add multi-port active impedance / active
   reflection coefficient extraction, Γ_active,m(θ,φ).** This is distinct
   from the S-parameter/impedance extraction in step 3 above: ordinary
   S-parameters characterize ports one or two at a time, but a phased array
   needs the impedance *seen by each element while every other element is
   simultaneously excited at its own scanned phase* — the quantity that
   predicts scan blindness. Build this as a post-processing step on top of
   the existing port machinery, not a parallel solve path; it's what makes
   the beam-steering sweep in Phase 07/12 (factorize once, cheap
   forward/backward substitution per angle) actually useful for array work.
   See `docs/THZ_PHASED_ARRAY_SCOPE.md`.

**Grounds this phase:** Dular et al. (2000), massive/stranded inductor source
formulations

**Ready for 08 when:** the EDA test problem produces port-level quantities
(R, L, or S-parameters) within a sanity-check tolerance of a hand estimate.

---

## 07 · Scattering track — excitation and open boundary

**Track:** Scattering, Full-wave

**Goal:** Incident-field excitation and a boundary treatment that lets a
bounded mesh stand in for an open (unbounded) region — the part of this
project most different in character from EDA. **Updated:** since Phase 01 now
keeps the full (all-frequency) A-Φ system, this track targets real radiating
problems (a dipole, a sphere at an actual operating frequency), not only
quasi-static scattering — the open-boundary approximation below is an interim
step, upgraded to an exact treatment in Phase 12.

1. **Implement plane-wave and dipole incident-field excitation** across your
   target frequency range — the all-frequency-stable formulation from Phase
   01 is what makes this valid beyond low-frequency, not a limitation to work
   around anymore.
2. **Add an outer boundary treatment.** **Updated (Sept 2026): implement a
   first-order ABC (Silver-Müller/Sommerfeld type), not PML.** PML's known
   low-frequency degeneracy would put a second low-frequency failure mode
   right next to a formulation whose whole premise is that it doesn't have
   one — a direct conflict with the all-frequency-stability pitch driving
   every other choice in this project. PML is also the larger implementation
   lift here: it needs an added volumetric absorbing layer (mesh changes),
   anisotropic complex-tensor assembly inside it, and empirical
   thickness/grading tuning, on top of translating a technique whose
   published form is for E/H fields onto the coupled A-Φ system while
   preserving the chosen gauge inside the layer. A first-order ABC costs a
   short field-to-potential substitution instead (E = −∂A/∂t − ∇Φ applied to
   the standard Silver-Müller condition) and drops in as a single boundary
   term in the existing weak form, no new DOFs or mesh changes. Its
   trade-off — accuracy that degrades away from near-normal incidence — is a
   documented, bounded limitation, not a new failure regime. No published
   ABC or PML formulated directly in A and Φ was found in a literature
   check (see `docs/OPEN_BOUNDARY_ABC.md` for the search and the from-scratch
   derivation used instead); this is expected to be interim and
   documented, since Phase 12 replaces it with an exact boundary-integral
   coupling once the core solver is proven.
3. **Add near-to-far-field / RCS post-processing** on top of the field
   solution.
4. **Run your Phase 01 scattering test problem** (e.g. a dipole or PEC/lossy
   sphere at a real frequency) and check the induced-field pattern and
   recoverable cross-section against the known analytical solution (Mie
   series, or a canonical dipole result) — this is also your baseline for
   the Phase 12 comparison later.
5. **Updated (Sept 2026): add Floquet periodic boundary conditions on a
   single unit cell** (infinite-array approximation, applied to both **A**
   and Φ) — the standard first modeling step for any phased array (scan
   blindness estimation, central-element active impedance) and far cheaper
   than the finite-array FEM-BI hybrid in Phase 12. Sequence this *before*
   committing to the full boundary-integral work: it's a fast way to
   validate the A-Φ core against COMSOL on a real periodic-array problem,
   using the Phase 06 active-impedance extraction, while Phase 12's dense
   boundary-integral machinery is still being built. See
   `docs/THZ_PHASED_ARRAY_SCOPE.md`.

**Grounds this phase:** Zhao & Fu (2017); Yan (2021) (frequency range this
formulation now targets)

**Ready for 08 when:** the scattering test problem matches its analytical
solution within a documented tolerance across the frequency range you've
committed to, with the ABC choice and its known angle-dependent accuracy
limits written down.

---

## 08 · Verification, validation, and a COMSOL cross-check

**Track:** Core

**Goal:** A standing regression/verification suite you (and eventually
collaborators) trust, plus the first honest comparison against COMSOL.

1. **Build a small analytical-benchmark suite**: solenoid/coaxial inductor
   (closed-form L), conducting slab skin effect (closed-form field profile),
   low-frequency sphere scattering (Mie limit) — wire these into the
   `tests/` CMake target as automated regression tests.
2. **Cross-validate against COMSOL** on 2–3 canonical problems spanning both
   tracks; record relative error and runtime side by side. This is your
   first real competitive data point, not marketing — be honest about where
   you're behind.
3. **Turn every bug found here into a permanent regression test** before
   moving on.

**Ready for 09 when:** the automated test suite passes and you have a
written, numeric COMSOL comparison table for at least one EDA and one
scattering case.

---

## 09 · Coulomb gauge track

**Track:** Coulomb gauge

**Goal:** A second, independently selectable gauge, added once tree-cotree is
proven correct — not a replacement, a documented alternative with its own
trade-offs.

1. **Implement the generalized Coulomb gauge via divergence penalty** (Li,
   Sun, Dai & Chew 2015): expand the divergence-of-**A** term through nodal
   (Whitney-0) elements while keeping **A** itself in edge elements — this
   avoids decomposing **A** into scalar components, which was the
   traditional Coulomb-gauge failure mode.
2. **Reuse your Phase 04 assembly machinery**: the new gauge only adds the
   **G**_NN/**G**_NE divergence-coupling blocks (their eqs. 13–18), not a
   parallel solver stack.
3. **Evaluate the sparse-approximate-inverse (SAI) trick** for the
   **G**_NN⁻¹ term if the dense inverse becomes a bottleneck.
4. **Optionally, evaluate the auxiliary-potential all-frequency formulation**
   (Yan 2021): an implicit Coulomb gauge via a third potential ψ, stable from
   DC to microwave frequencies without tree-cotree splitting at all — worth
   a spike to see whether it simplifies your EDA and scattering tracks into
   one unified formulation.
5. **Re-run the Phase 08 benchmark suite under Coulomb gauge** and compare
   condition number, iteration count, and accuracy directly against the
   tree-cotree results.

**Grounds this phase:** Li, Sun, Dai & Chew (2015); Yan (2021), §2

**Ready for 10 when:** both gauges pass the same benchmark suite and you have
a written recommendation for which gauge is the default per problem class.

---

## 10 · Performance and scaling

**Track:** Core

**Goal:** The solver stops being a correctness prototype and becomes
something that can take on COMSOL-sized meshes in reasonable time.

1. **Profile before optimizing**: assembly vs. solve time breakdown on a
   realistic mesh size for each track.
2. **Parallelize assembly** (element-level parallelism is embarrassingly
   parallel) and evaluate a parallel/multithreaded sparse solver or
   preconditioner.
3. **Add true hierarchical/higher-order basis functions and hp-adaptivity.**
   Phase 02 settled on a mixed-order pairing (first-order **A**,
   second-order Φ; Sept 2026) rather than raising both fields together, so
   this phase now covers two possible upgrades, either or both: (a) a
   matched second-order **A** (Graglia, Wilton & Peterson 1997 and/or
   García-Castillo et al. 2000, still unread, identified as candidates), to
   remove the A-Phi coupling accuracy cap documented in
   `docs/FORMULATION.md` §5.1; and (b) true hierarchical/hp-adaptive bases
   beyond second order generally. The general principle (matching nodal and
   edge order improves accuracy per DOF and stability in the h-adaptive
   regime) follows Lee 2003, though that paper's own scheme is 2-D
   triangular and not directly reusable here — a 3-D hierarchical reference
   will need to be identified when this phase is reached.
4. **Evaluate GPU or matrix-free approaches** for the largest target
   problems, once CPU performance is well understood and profiled — don't
   start here.

**Grounds this phase:** Lee, Lee & Lee (2003), §III–IV (hierarchical bases,
hp-adaptivity)

**Ready for 11 when:** you have a documented scaling curve (time / memory vs.
problem size) for both tracks.

---

## 11 · Productization — competing with COMSOL

**Track:** Core

**Goal:** Turn a validated, performant solver into something a third party
(or a future customer) could actually use, and settle the licensing question
deferred from Phase 00.

1. **API / scripting layer**: a stable public interface (C++ and, likely,
   Python bindings) for defining geometry, materials, sources, and running
   solves — this is what a "product" needs that a research code doesn't.
2. **Standard geometry/mesh import** (STEP for CAD geometry, common mesh
   formats) so users aren't limited to hand-built test cases.
3. **Results post-processing and visualization**, or export to an existing
   viewer (ParaView via VTK output is a reasonable first step rather than
   building your own).
4. **Decide the license and packaging model** now that there's real IP to
   protect: closed-source commercial, source-available, or dual-licensed —
   this determines a lot about how you distribute and support it.
5. **Stand up a continuous COMSOL-comparison benchmark** (extending Phase
   08) that runs on every significant change — this becomes both your
   internal quality bar and, eventually, your sales material.
6. **Updated (Sept 2026): sharpen the beachhead target.** Rather than a
   generic "competes with COMSOL on EDA and scattering" pitch, the named
   target is multiscale sub-THz/THz systems where sub-micron semiconductor
   feeds must be co-designed with radiating apertures — photoconductive THz
   antenna arrays, on-chip antenna-circuit co-design, and automotive radar
   packages are the example verticals. This changes the pitch and the
   benchmark problems worth showcasing in step 5, not what gets built first;
   Phases 00–10 are unchanged. See `docs/THZ_PHASED_ARRAY_SCOPE.md`.

**Definition of done for the FEM-only solver (v1):** a documented,
reproducible benchmark set where the solver matches or beats COMSOL on
accuracy for your target problem classes, with a runtime/usability story you
can put in front of a first customer. Phase 12 below is the next major bet,
not a blocker for this one.

---

## 12 · Full-wave open boundary — FEM-BI hybrid (deferred)

**Track:** Full-wave, Scattering — **Deferred**

**Goal:** Replace the ABC approximation from Phase 07 with an exact,
non-approximate open boundary by coupling the interior FEM solve to a
boundary integral equation on the truncation surface. Deliberately deferred:
this is a substantial, separate undertaking on the scale of what
distinguishes the top commercial tier of full-wave solvers (HFSS's hybrid
regions, FEKO's MLFMM) — start it once Phases 00–11 give you a working,
validated FEM-only solver, not before. **Updated:** a focused literature
check turned up a specific, better-motivated coupling choice than a generic
EFIE/MFIE/CFIE pairing, and a plausible angle for a standalone paper
alongside the product work.

1. **Prefer a potential-based boundary formulation over classical
   field-based SIEs.** Most FEM-BI work couples to EFIE/MFIE (with CFIE as
   the standard fix for interior-resonance artifacts on closed surfaces) —
   but EFIE has its own well-known low-frequency breakdown, separate from
   the one A-Φ already solves on the FEM side, so a classical-SIE coupling
   can still degrade at low frequency even with a perfect interior solve.
   Pairing the all-frequency-stable A-Φ FEM with a *potential-based* BEM
   instead (Sharma & Triverio 2021, 2022) keeps the whole hybrid uniformly
   stable from DC through full-wave, and couples more naturally at the
   interface: both sides are already expressed in **A** and Φ, so there's no
   need to reconstruct E/H just to hand quantities across the boundary.
   Keep CFIE as a documented fallback if the potential-based route runs
   into trouble specific to your geometry classes.
2. **Triangulate the truncation boundary** — the surface enclosing all
   inhomogeneous/active regions — separately from the interior volume mesh,
   and implement the boundary-element matrices for whichever formulation
   you land on, including proper singular/near-singular integral quadrature
   (a well-known but fiddly numerical detail in its own right).
3. **Implement the FEM-BI coupling at the interface** (a
   Robin/Steklov–Poincaré-type condition linking the FEM's natural boundary
   condition to the BI operator) — most FEM-BI implementation bugs live
   here. Validate against a case with a known closed-form solution (a
   dielectric or PEC sphere) before trusting it on anything else.
4. **Address the dense BI matrix's cost before it becomes a bottleneck.** A
   direct dense solve is fine for small truncation surfaces, but real
   antenna/RCS-scale problems need a fast method to avoid O(N²)/O(N³)
   scaling. **Updated (Sept 2026): prefer Adaptive Cross Approximation (ACA)
   over MLFMA as the default**, not just one of two equally-weighted
   options — MLFMA's analytic multipole expansion of the Helmholtz kernel is
   well known to degrade at sub-wavelength scales (the spherical Hankel
   functions involved become ill-conditioned as their argument shrinks),
   which is exactly the regime a multiscale THz array with sub-micron feed
   detail lives in. ACA is purely algebraic (a low-rank cross-approximation
   on kernel-evaluated matrix blocks), kernel-independent, and has no
   analogous low-frequency failure mode — consistent with the
   all-frequency-stability reasoning behind every other formulation choice
   in this roadmap. This is itself a substantial subsystem; budget for it as
   such rather than as a footnote. See `docs/THZ_PHASED_ARRAY_SCOPE.md`.
5. **Re-run the Phase 08 benchmark suite** (extended with radiating/open-region
   analytical cases) under the FEM-BI hybrid and compare accuracy and cost
   against the Phase 07 ABC baseline.
6. **Treat this as a candidate paper, not only an implementation task.** A
   literature check (Sept 2026) found close but distinct prior work: a
   broadband A-Φ solver using discrete exterior calculus with no BI coupling
   (Zhang, Na, Jiao & Chew 2022); a full-wave, DC-to-high-frequency
   potential-based BEM that is pure BEM with no FEM coupling (Sharma &
   Triverio 2021, 2022); and older A-V/BEM–FEM couplings limited to
   magnetostatic/eddy-current problems. An all-frequency-stable A-Φ FEM
   coupled to a potential-based BEM, uniformly stable end-to-end, did not
   turn up in that search — treat that as a lead, not a confirmed gap, and
   do a fuller, citation-checked lit review before writing up a novelty
   claim.

**Grounds this phase:** Sharma & Triverio, potential-based BEM for lossy
conductors/dielectrics, DC to high frequency (2021, 2022); Zhang, Na, Jiao &
Chew, DEC-based broadband A-Φ solver (2022, for contrast — no BI coupling);
Jin & Volakis, hybrid finite-element–boundary-integral method for scattering
and radiation; general EFIE/MFIE/CFIE background (Peterson, Ray & Mittra,
*Computational Methods for Electromagnetics*) — general numerical-methods
literature, not specific to any one implementation.

**Ready to fold into the product when:** the FEM-BI hybrid matches or beats
the Phase 07 ABC/PML baseline in accuracy on radiating test cases, and the
dense-matrix cost is under control (FMM/ACA, or a documented size limit) for
your target problem sizes.

---

Built from the public literature in `APhi_Papers/` only — see
`docs/REFERENCES.md` for the running citation list. Revisit and re-order
phases as the two application tracks (EDA, scattering) reveal which needs
attention first. Phase 12 (FEM-BI hybrid) is a deliberate, later bet — it does
not block shipping the FEM-only solver first.
