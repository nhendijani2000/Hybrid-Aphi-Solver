# Engineering standards for this codebase

*Established Sept 2026, applies from Phase 03 onward (and to any future
rework of Phase 02).*

The project owner's explicit standing direction for all code
implementations from here on:

1. **Optimize for speed first, memory usage second.** These trade off
   against each other sometimes; when they do, prefer the faster option.
2. **Use parallelization wherever it is genuinely appropriate** — not
   reflexively everywhere. A loop that's already sub-millisecond, or an
   algorithm with a hard sequential data dependency, should stay
   single-threaded rather than gain a pragma that adds complexity and
   synchronization overhead for no measurable benefit. Each phase's code
   should say explicitly, in a comment, why a hot path was or wasn't
   parallelized.
3. **Use advanced sparse-matrix techniques and operations** rather than the
   simplest-possible approach, once the code is doing real linear algebra at
   production scale (Phase 04 assembly onward) — compressed sparse formats
   (CSR/CSC), fill-reducing orderings, sparse direct factorization
   (e.g. via SuiteSparse) or Krylov iterative solvers with a real sparse
   matrix-vector product, rather than the COO/`std::map`-based `SparseMatrix`
   in `incidence.cpp`. That type is explicitly a Phase-02-scale placeholder
   (its own header comment says so) — Phase 04/05 is where this standard
   actually bites.

## What this means concretely, phase by phase

- **Phase 02/03 (mesh, incidence, tree-cotree)** is graph-topology work, not
  numerical linear algebra: the relevant "speed" work is cache-friendly data
  structures (flat/CSR arrays instead of `std::map`-based lookups in hot
  traversal loops) and asymptotically optimal graph algorithms (BFS/DFS for
  spanning trees, Union-Find with path compression for connected-component
  grouping) — already asymptotically optimal, so there is no "more advanced"
  algorithm to reach for here, only a better-constant-factor implementation.
  Parallelizing a single BFS/Union-Find pass over a mesh graph has a hard
  sequential dependency between BFS levels and is not worth it at realistic
  FEM mesh sizes (thousands to low millions of nodes) — the traversal is
  already memory-bandwidth-bound and sub-millisecond to low-millisecond;
  true parallel BFS (level-synchronous, direction-optimizing) only pays for
  itself at graph-analytics scale (hundreds of millions to billions of
  edges), which this project's meshes are nowhere near. Where parallelism
  *does* pay off even at Phase 02/03 scale: independent per-element work
  (per-tet geometry, per-tet basis evaluation) is embarrassingly parallel
  and a natural `#pragma omp parallel for` candidate once there are enough
  tets that the loop overhead is no longer dominated by fixed cost.
- **Phase 04 (assembly)** is where per-tet element-matrix computation
  (embarrassingly parallel across tets) and scatter-into-global-matrix
  (needs either thread-local accumulation + merge, or atomics, to stay
  correct under parallelism) both matter.
- **Phase 05 (linear solve)** is where the sparse-matrix standard matters
  most: a real sparse format, a real factorization or iterative method, and
  parallel sparse matrix-vector products, replacing the placeholder
  `SparseMatrix`/`multiply` in `incidence.cpp` for anything beyond the
  `CG = 0` topology check it was built for.
- **Exception, made deliberately (Sept 2026):** `estimate_condition_number`
  (`src/gauge_variants.cpp`) was originally a dense implementation --
  correct, but O(n^2) memory and O(n^3) time per inverse-power-iteration
  step, fine only at Phase 03's tiny (single-digit-DOF) validation-mesh
  scale. It was rewritten to be sparse (matrix-free `A^T*A` matvecs, a
  Conjugate-Gradient inner solve instead of dense Gaussian elimination --
  see `docs/REFERENCES.md`, "Matrix conditioning") ahead of its originally-
  planned Phase 04/05 slot, specifically so a gauge-comparison tool (letting
  the user choose and compare the Albanese-Rubinacci vs. Munteanu-
  unsymmetric tree-cotree gauge on real mesh sizes) can be built without
  hitting that O(n^3) wall immediately. This is a narrow, motivated
  exception for one specific function, not a general "start using advanced
  sparse techniques now" -- the rest of Phase 03/04's own machinery
  (`SparseMatrix`/`multiply` in `incidence.cpp`, dense per-tet element
  matrices) is still exactly where the phase-by-phase guidance above says
  it should be.
- **Superseded (Sept 2026), for the two bullets above and item 3's
  "rather than the COO/`std::map`-based `SparseMatrix` in `incidence.cpp`":
  that placeholder no longer exists.** Phase 03.5 (`docs/ROADMAP.md`)
  replaced it with `Sparse<T>` in `include/aphi_solver/sparse_matrix.hpp`
  — triplet accumulation, one `compress()` into CSR, templated on the
  scalar so the real incidence operators and the complex assembled blocks
  share one implementation. `incidence.hpp`'s `SparseMatrix` is now an alias
  for `Sparse<double>`, `APhiBlockSystem` holds sparse blocks, and the
  hand-written `coalesce`/`multiply`/`transpose` and the local matvec
  helpers in `gauge_variants.cpp` were deleted rather than ported. So item 3
  is satisfied ahead of Phase 05, not still pending. What the text above
  still describes correctly is the *reasoning*: it was done when the code
  started doing real linear algebra, which turned out to be the Phase 03/04
  boundary rather than Phase 05. Dense per-tet element matrices remain
  appropriate and unchanged — those are 6x6/10x10, not a scaling concern.

## Style

- Prefer flat contiguous containers (`std::vector`, `std::array`,
  CSR-style offset/value arrays) over node-based containers
  (`std::map`, `std::set`) in any loop that runs once per mesh entity
  (node/edge/face/tet). `std::map` stays fine for one-off/rare lookups
  (e.g. `Mesh::find_edge`, called a handful of times per face during
  incidence-matrix assembly, not once per traversal step).
- New OpenMP-parallelized loops should be guarded by `#ifdef _OPENMP` so the
  code still builds correctly (serially) in a toolchain without OpenMP
  configured (e.g. this project's current MSVC/Visual Studio build has not
  yet had `/openmp` enabled in `CMakeLists.txt`).
- Every performance-motivated design choice (a data structure, an algorithm,
  a decision not to parallelize something) gets a short comment explaining
  the reasoning, the same way the rest of this codebase documents its
  formulation choices — so the next person reading the code (including a
  future Claude session) doesn't have to re-derive why.
