# A-Phi Solver Roadmap

*From formulation to a COMSOL-class full-wave solver*

A phased technical roadmap for an all-frequency-stable A-Φ (magnetic vector
potential / electric scalar potential) finite-element solver — valid from DC
through full-wave, not just low-frequency eddy-current problems — gauged first
by tree-cotree splitting and later also by Coulomb gauge, aimed at EDA
(signal/power-integrity) and scattering applications across that full frequency
range. An exact open boundary via a FEM–boundary-integral hybrid is planned as a
deferred, later phase once the core solver works end to end.

| | |
|---|---|
| **Owner** | Nastaran Hendijani |
| **Language** | C++17, CMake |
| **Repo** | `APhi_Solver_Project_LowFrequency_EDA/APhi_Solver` |
| **Phases** | 13 (00–12, last deferred) |

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
2. **Install Visual Studio 2022 Community.** Download from
   `visualstudio.microsoft.com`. In the installer, select the **"Desktop
   development with C++"** workload — this pulls in the MSVC compiler, the
   Windows SDK, and CMake tools for Visual Studio.
3. **Install the GitHub Copilot extension.** Inside Visual Studio:
   `Extensions → Manage Extensions → search "GitHub Copilot" → Download`,
   restart Visual Studio to finish installing, then `Extensions → GitHub
   Copilot → Sign in` with your GitHub account (create one at github.com if
   you don't have one — a free account is enough to start).
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
   use internally — that's fixed in Phase 07 (ABC/PML now, exact FEM-BI
   coupling later in Phase 12), not here.
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
3. **Implement lowest-order Whitney edge elements** for **A** (tangential
   continuity) and nodal (Lagrange P1) elements for Φ (the `F¹_A` / `F⁰_Φ`
   spaces).
4. **Defer higher-order/hierarchical bases** (Lee 2003's hierarchical hp
   scheme) to Phase 10 — get a correct lowest-order pipeline working
   end-to-end first.

**Grounds this phase:** Munteanu (tree-cotree condensation properties); Lee,
Lee & Lee (2003), §IV–VI

**Ready for 03 when:** **C** and **G** can be assembled for a test mesh and
**CG = 0** holds numerically (the discrete curl·grad = 0 identity) — this is
your first real unit test.

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

1. **Implement the spanning-tree search** over the mesh graph (Lee 2003
   Algorithms 1–2: node numbering with PEC/ground handling, then tree/cotree
   edge marking) — this is graph theory on your edge/node incidence lists, no
   field solve required yet.
2. **Implement at least two gauging variants** from Munteanu's projection
   framework so you can compare them empirically rather than trust one
   blindly: **Albanese–Rubinacci** (simplest, keeps sparsity, but worse
   conditioning) and **Munteanu unsymmetric** (similarly simple, keeps the
   original condition number best per her numerical tests).
3. **Instrument condition-number reporting** from day one (even a crude
   power-iteration estimate) — this single number is what will guide every
   gauge decision for the rest of the project.
4. **Handle multiple PEC bodies** explicitly (each grounded to its own
   reference node) — a common source of subtle bugs later.
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

## 04 · Weak form and matrix assembly

**Track:** Core

**Goal:** A sparse assembly pipeline producing the coupled A-Φ block system,
with the symmetric row-scaling option available as a solver-conditioning knob.

```
K_AA  = ∫ (∇×A)·(∇×A′) dΩ        M_AA  = ∫ A·A′ dΩ
K_AΦ  = ∫ ∇Φ·A′ dΩ                K_ΦΦ = ∫ ∇Φ·∇Φ′ dΩ
```

1. **Element-level matrices first**, tested on a single reference tet
   against hand-computed values, before any global assembly code runs.
2. **Global sparse assembly** into a CSR (or similar) structure, respecting
   the tree-cotree reduction from Phase 03 (reduced-size system, not the full
   singular one).
3. **Implement the symmetric row-scaling option** (divide the Φ-row by jω) as
   a build-time or run-time flag — keep both the natural and symmetric forms
   available since later phases compare them.
4. **Right-hand-side assembly** for the source models you'll add properly in
   Phases 06–07; a simple uniform current density is enough here to exercise
   the pipeline.

**Grounds this phase:** general A-Φ weak form, standard in the cited
literature (Dular et al. 2000; Zhao & Fu 2017)

**Ready for 05 when:** the assembled sub-blocks are complex-symmetric where
they should be (K_AA, K_ΦΦ symmetric; K_AΦ and K_ΦA transposes of one another
up to their respective coefficients, per `docs/CONDITIONING.md`), and a
trivial problem (e.g. straight wire, known field) gives the right sign and
rough magnitude of **B**.

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
2. **Add an outer boundary treatment: ABC or PML, as an interim, approximate
   open boundary.** Either has its own error budget (an ABC is only accurate
   for a given incidence-angle range; PML adds volumetric DOFs and has known
   low-frequency degeneracy issues of its own) — document which you chose
   and its accuracy limits, since Phase 12 replaces this with an exact
   boundary-integral coupling once the core solver is proven.
3. **Add near-to-far-field / RCS post-processing** on top of the field
   solution.
4. **Run your Phase 01 scattering test problem** (e.g. a dipole or PEC/lossy
   sphere at a real frequency) and check the induced-field pattern and
   recoverable cross-section against the known analytical solution (Mie
   series, or a canonical dipole result) — this is also your baseline for
   the Phase 12 comparison later.

**Grounds this phase:** Zhao & Fu (2017); Yan (2021) (frequency range this
formulation now targets)

**Ready for 08 when:** the scattering test problem matches its analytical
solution within a documented tolerance across the frequency range you've
committed to, with the ABC/PML choice and its known limitations written down.

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
3. **Add hierarchical/higher-order basis functions and hp-adaptivity**
   (deferred from Phase 02) — per Lee 2003, this both improves accuracy per
   DOF and improves stability in the h-adaptive regime.
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

**Definition of done for the FEM-only solver (v1):** a documented,
reproducible benchmark set where the solver matches or beats COMSOL on
accuracy for your target problem classes, with a runtime/usability story you
can put in front of a first customer. Phase 12 below is the next major bet,
not a blocker for this one.

---

## 12 · Full-wave open boundary — FEM-BI hybrid (deferred)

**Track:** Full-wave, Scattering — **Deferred**

**Goal:** Replace the ABC/PML approximation from Phase 07 with an exact,
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
   antenna/RCS-scale problems need a fast method — the fast multipole method
   or adaptive cross approximation — to avoid O(N²)/O(N³) scaling. This is
   itself a substantial subsystem; budget for it as such rather than as a
   footnote.
5. **Re-run the Phase 08 benchmark suite** (extended with radiating/open-region
   analytical cases) under the FEM-BI hybrid and compare accuracy and cost
   against the Phase 07 ABC/PML baseline.
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
