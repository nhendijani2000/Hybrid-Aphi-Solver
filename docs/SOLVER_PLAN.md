# Sparse linear solver — plan

Phase 05. The matrix exists and is cheap to rebuild (`docs/ASSEMBLY_PLAN.md`);
nothing solves it yet, so nothing in this project has been validated
numerically. Everything to date is structural or algebraic.

**Decided 26 Sept: in-house, no third-party library** (an optional, default-off
backend may be added later -- Sec. 12 -- but nothing in the build will ever
require one). This reverses the
"link MUMPS" recommendation in `docs/LINEAR_SOLVER.md`, which is marked
superseded there. The work is staged: **Stage 1** (§1–§9) is a correct,
deterministic, reference-quality direct solver; **Stage 2** (§10) is the
performance and robustness work, done in-house later *if measurements say it
is needed*. Nothing here depends on an external package at any stage.

---

## 1. What it has to do

| | |
|---|---|
| unsymmetric systems | `Conditioning::Natural` — LU |
| complex symmetric systems | `RowScaled`, `ScaledPhi` — `LDLᵀ` |
| frequency sweeps | one symbolic analysis, many numeric factorizations |
| diagonal scaling | `equilibration.hpp`, which already exists |
| fill-reducing reordering | new; the largest single lever on a direct solve |

**Complex symmetric, not Hermitian.** `A = Aᵀ` and `A ≠ Aᴴ`, and the system is
indefinite. Cholesky does not apply. `LDLᵀ` does, with `D` complex and no
positivity anywhere. This is asserted, not assumed: `test_assembly` checks
`A = Aᵀ` to the bit on both real meshes.

**The system is already reduced.** The gauge eliminated tree and Dirichlet
edges in the DOF map, so what the solver receives has no gradient null space to
work around. That is what makes a direct factorization meaningful at DC at all.

---

## 2. Architecture: symbolic once, numeric per frequency

The same split that paid for itself in assembly — 84.6 ms → 29.5 ms per
frequency, §13 there — and for the same reason: **the matrix's structure does
not depend on frequency, only its values do.** Measured, not assumed: across
1 kHz / 1 MHz / 10 GHz, `row_ptr` and `col_index` are identical while 1557013
of 1557522 values differ.

    Ordering        ->  permutation                    once per mesh
    analyze         ->  elimination tree, nnz(L)       once per mesh + ordering
    factorize       ->  numeric L, D (or L, U)         once per frequency
    solve           ->  triangular solves + refinement once per right-hand side

```cpp
enum class Ordering { Natural, ReverseCuthillMcKee, ApproximateMinimumDegree };

/// Everything about the factorization that does not depend on a value.
struct SolverAnalysis {
    Ordering ordering;
    std::vector<int> perm;       ///< new index -> old
    std::vector<int> iperm;      ///< old index -> new
    std::vector<int> parent;     ///< elimination tree, -1 at a root
    std::vector<int> col_count;  ///< nonzeros in each column of L
    SparsityPattern factor;      ///< L's pattern, in PERMUTED indices
    std::size_t predicted_nnz;   ///< sum of col_count; checked against reality
};

SolverAnalysis analyze(const SparsityPattern& a, Ordering ordering);
```

`analyze` takes the **symmetric pattern** in both cases: for `LDLᵀ` that is `A`
itself, for `LU` it is the pattern of `A + Aᵀ`. One ordering implementation
therefore serves both paths, which is worth the small amount of extra fill it
costs the unsymmetric case.

---

## 3. Reordering

The reason to do this at all: fill. Factorizing in the DOF map's natural order
(`[a | phi | V]`, which is edge and node numbering from the mesh file) produces
far more nonzeros in `L` than the matrix has, and both memory and flops scale
with that.

Three orderings, because a measurement needs a baseline:

| | what it is | why it is here |
|---|---|---|
| `Natural` | no permutation | **the control.** If AMD does not beat this substantially, the ordering is not working |
| `ReverseCuthillMcKee` | bandwidth reduction (Cuthill & McKee 1969, reversed by George) | ~100 lines, a sanity check that any ordering helps, and a fallback |
| `ApproximateMinimumDegree` | AMD (Amestoy, Davis & Duff 1996, 2004) | the standard fill-reducing ordering and the one to use |

All three are published algorithms. Nested dissection is Stage 2 (§10).

**The first deliverable of this whole plan is one table:** `nnz(L)` and
factorization time for all three orderings on `cylinder_box.msh` and
`loop_cut.msh`, symmetric and unsymmetric. That table decides where the direct
ceiling sits for this problem class, and it replaces the order-of-magnitude
estimate below with a number.

**Estimated, NOT measured:** for 3D FEM, `nnz(L)/n` is commonly 100–1000, so at
37368 unknowns expect 4–40 M nonzeros in `L`, i.e. 60–600 MB of complex values.
Feasible. At the 15–25 M DOFs that `docs/THZ_PHASED_ARRAY_SCOPE.md` §43 gives
as an order-of-magnitude figure for a 16×16 array, the same arithmetic gives
10¹⁰–10¹¹ — terabytes. **A direct solve does not reach that size, whoever
writes it.** That is an argument for the iterative phase (§11), not for a better
factorization.

**A memory budget, checked before allocating.** `analyze` knows
`predicted_nnz` before any value is touched, so the solver refuses a
factorization that would exceed a stated budget, naming the number, instead of
attempting a 40 GB allocation.

---

## 4. Scaling, and getting the un-scaling right

`compute_symmetric_equilibration` already exists: real positive `d`, Ruiz-style,
applied as `diag(d) · A · diag(d)`, which preserves `A = Aᵀ` — that is why it is
symmetric rather than independent row and column scalings.

The order is **scale, then order, then factorize**. Ordering is structural, so
scaling cannot affect it; scaling improves pivot quality, so it must precede
factorization.

The step that is easy to get wrong, written out because it will be tested:

    Ã = D A D,   solve  Ã y = D b,   then  x = D y

Not `x = y`, and not `D⁻¹ y`. A missing or inverted `D` on the way out gives a
plausible, smooth, completely wrong field. One negative control does exactly
that.

**This is equilibration's first real test.** A voltage port's constraint row
carries a unit diagonal beside matrix entries of ~3.4e11 (`ASSEMBLY_PLAN` §10) —
one badly scaled row, deliberately left for equilibration to fix, with no way to
measure the cost until now.

---

## 5. Pivoting: static, with refinement

**Decision: static pivoting plus iterative refinement, not dynamic pivoting.**

`LDLᵀ` with no pivoting can hit a zero pivot on a nonsingular complex symmetric
matrix — `[[0,1],[1,0]]` is the minimal example — so something must handle small
pivots. Two routes:

- **Dynamic** (Bunch-Kaufman, 1×1 and 2×2 pivots): always stable, and it
  **changes the fill pattern during numeric factorization**. That destroys the
  symbolic-once property of §2, so a 41-point sweep would pay for analysis 41
  times. It is also the single most intricate part of a sparse direct solver.
- **Static** (fixed ordering; perturb a pivot below a threshold; recover accuracy
  by refinement): keeps the pattern fixed and the sweep cheap, and is a few dozen
  lines. The cost is that accuracy now depends on refinement converging.

Static, therefore — and the honest consequence is that the solver **must report
what it did**:

```cpp
struct SolveReport {
    double residual_before_refinement;   ///< ||Ax - b|| / ||b||, original A
    double residual_after_refinement;
    int refinement_steps_used;
    int perturbed_pivots;                ///< 0 is the normal case
    double smallest_pivot_magnitude;
    std::size_t factor_nnz;              ///< actual; must equal predicted_nnz
    double analyze_ms, factorize_ms, solve_ms;
};
```

A perturbed pivot is not hidden and not fatal: it is counted, and the residual
says whether it mattered. Refinement that fails to converge is an error, not a
shrug.

**The residual is computed against the original `A` and `b`** — unscaled and
unpermuted. That makes it a true backward error rather than a statement about
the factorization's own arithmetic.

---

## 6. Factorization

**Up-looking, scalar.** The `LDLᵀ` factorization proceeds row by row using the
elimination tree's row subtrees — the approach of Davis, *Algorithm 849: A
concise sparse Cholesky factorization package* (ACM TOMS 2005), which is
compact enough to be read and checked in full. Correctness over speed is the
Stage 1 brief; blocking is Stage 2.

**The symmetric path factors the upper triangle in place of a transpose.**
`SparseSymmetric` stores the upper triangle, so the factorization is written as
`A = Uᵀ D U` with `U` unit upper triangular, rather than converting to lower
first. `find_slot` normalises index order, so no caller needs to know which
triangle is held.

**The unsymmetric path** is `LU` with the same ordering machinery on `A + Aᵀ`'s
pattern and the same static pivoting. Its fill pattern is computed from that
symmetric pattern, which over-estimates slightly and in exchange keeps one
symbolic implementation.

**At DC the matrix is block triangular**, and that is worth knowing rather than
exploiting: `α = 0` makes the entire `(Φ,A)` block exactly zero
(`test_assembly` asserts this), so

        [ K_AA   βC  ]
        [  0     βL  ]

Φ solves independently, then A follows. A general factorization handles this
without special-casing. It also gives a free independent check: solve the Φ
block alone and the answer must match the full solve's Φ entries.

---

## 7. The interface

```cpp
struct SolveOptions {
    Ordering ordering = Ordering::ApproximateMinimumDegree;
    int equilibration_iterations = 10;   ///< 0 disables scaling
    int max_refinement_steps = 3;
    double pivot_threshold = 1e-14;      ///< relative to the scaled matrix
    std::size_t memory_budget_bytes = 0; ///< 0 = no limit
};

/// Symmetric path. `a` holds the upper triangle.
SolveReport solve(const SparseSymmetricZ& a, const std::vector<Complex>& b,
                  std::vector<Complex>& x, const SolveOptions& options);

/// Unsymmetric path.
SolveReport solve(const SparseMatrixZ& a, const std::vector<Complex>& b,
                  std::vector<Complex>& x, const SolveOptions& options);
```

and, for a sweep, the split form that reuses the analysis — the whole point of
§2:

```cpp
SolverAnalysis analysis = analyze(pattern, options.ordering);
SymmetricFactor factor;
for (double f : problem.frequencies) {
    refill(system, ...);                                  // ASSEMBLY_PLAN §13
    factorize(system.matrix, analysis, factor, report);   // per frequency
    solve(factor, system.rhs, x, report);
}
```

`[solver]` in the input file gains `ordering = amd | rcm | natural` and
`refinement = N`, beside the existing `conditioning`. Which factorization runs
is **derived** from `conditioning`, not a separate field a user could set
inconsistently — `Natural` implies LU, the two symmetric forms imply `LDLᵀ`.
That rule already exists in `LINEAR_SOLVER.md` and still holds.

---

## 8. Verification

The point of Stage 1. In rough order of strength:

1. **Predicted `nnz(L)` == actual.** Symbolic and numeric agree or one of them
   is wrong. Free, and it catches most symbolic bugs immediately.
2. **Against a dense solve.** On small systems, the sparse solver and a dense
   LU must agree to machine precision. `complex_matrix.hpp` already provides the
   dense path.
3. **Residual, always.** `‖Ax − b‖/‖b‖` on the original system, reported in
   `SolveReport` on every solve, not only in tests.
4. **`LDLᵀ` and `LU` on the same matrix must agree.** The symmetric forms can be
   solved both ways; two independent factorizations reaching the same answer is
   the strongest check available here.
5. **The three conditionings must give the same physical answer.** They describe
   one system (`test_assembly` proves the matrices are related exactly), so
   after undoing the scaling their solutions must match — including `Φ = jω·Φ'`
   for `ScaledPhi`. This is an end-to-end check the assembly work already set up
   and it exercises the whole chain.
6. **The DC milestone.** `R = 0.1388 mΩ` exact against the meshed cross-section,
   `L = 0.3870 nH` within ~0.5 % — `examples/cylinder_box.aphi` states both and
   where they come from.
7. **Port current self-consistency.** The cylinder's P2 must report −1 A when P1
   drives +1 A. A voltage port's current is no longer a row residual
   (`ASSEMBLY_PLAN` §10), so extraction has to re-form the terminal's current
   balance — that lands with this phase.
8. **Ordering is a permutation.** `perm` and `iperm` are mutually inverse
   bijections, and `P A Pᵀ` has exactly as many nonzeros as `A`.
9. **Permuting is invisible.** Solving the permuted system and un-permuting must
   equal solving directly, to the bit where the arithmetic order allows.

### Negative controls

Each must be shown to fail:

- skip the `x = D y` un-scaling, or invert it
- permute the right-hand side but not the solution, or use `iperm` for `perm`
- transpose the permutation in `P A Pᵀ`
- drop iterative refinement while static pivoting perturbs a pivot
- use the natural ordering and assert it produces **more** fill than AMD (if it
  does not, AMD is not doing anything)
- drop the elimination-tree row subtree and factor a dense row instead
- forget that `ScaledPhi`'s unknown is `Φ'` and report it as `Φ`
- solve with the unsymmetric path on a symmetric matrix stored as one triangle
  (must be refused, not silently half-solved)

---

## 9. Implementation order

Each step verifiable before the next exists.

1. **RCM**, with the permutation invariants of §8.8. Cheapest possible start and
   it makes the ordering interface concrete.
2. **Elimination tree and symbolic factorization**, for a given permutation.
   Deliverable: `nnz(L)` for Natural and RCM on both real meshes.
3. **AMD.** Deliverable: the §3 table, which is where the direct ceiling stops
   being a guess.
4. **`LDLᵀ` numeric, no pivoting**, on a small symmetric matrix, against dense.
5. **Triangular solves and the residual.** First end-to-end solve.
6. **Equilibration wired in**, with the un-scaling test and its control.
7. **Static pivoting and refinement**, with `SolveReport` populated.
8. **`LU`** for the unsymmetric path, against the same dense reference.
9. **The DC milestone**, `R` and `L`, on `examples/cylinder_box.aphi`.
10. **The three-conditioning agreement check**, on both real meshes.
11. **Time it**, then decide about Stage 2 with numbers in hand — the same rule
    that governed threading in `ASSEMBLY_PLAN` §8 step 8, where the measurement
    said the obvious candidate was the wrong one.

---

## 10. Stage 2 — later, in-house, only if measured to be needed

Everything here is a **performance or robustness** upgrade, not a new
capability. None of it is required for a correct answer, and none of it needs a
third-party library. Each item has a trigger: do it when a measurement says so,
not on principle.

| | what it buys | trigger |
|---|---|---|
| **Nested dissection** ordering | better fill asymptotics than minimum degree on 3D meshes; raises the size ceiling perhaps 2–5× in `n` | AMD's `nnz(L)` becomes the binding memory constraint on a problem that must be solved directly |
| **Supernodal / BLAS3 blocking** | roughly 5–20× on factorization time, by turning scalar updates into dense block operations | factorization time dominates a sweep after the ordering work is done |
| **Bunch-Kaufman pivoting** | guaranteed-stable symmetric indefinite factorization | `SolveReport::perturbed_pivots` is regularly non-zero *and* refinement fails to recover the residual. Note the cost: it breaks the symbolic-once property of §2, so a sweep would re-analyse per frequency |
| **Multithreaded factorization** | wall-clock, via the elimination tree's independent subtrees | single-core factorization is the bottleneck and blocking is already done |
| **Out-of-core factors** | problems whose `L` exceeds RAM | rare; usually the signal to switch to iterative instead |

**Keep Stage 1's solver after Stage 2 exists.** Not as dead code — as the
reference. A scalar, single-threaded, statically pivoted factorization is
**bitwise deterministic**: same bits every run, every machine. A blocked,
threaded one generally is not, because the summation order changes. This project
asserts bitwise equality in several places already (the `ScatterMap` path,
`to_full`, symmetric storage), and a deterministic reference is what makes those
assertions possible. Stage 2 is then validated against Stage 1 the way the
pattern assembly was validated against the triplet path.

---

## 11. What this plan does NOT address

**The large-scale path is iterative, and it is a separate phase.** At the
15–25 M DOFs `THZ_PHASED_ARRAY_SCOPE` gives as an order-of-magnitude figure, no
direct factorization fits in memory — Stage 2 does not change that, it moves the
constant. The scalable route is COCG for the complex symmetric forms and
BiCGStab or GMRES for the unsymmetric one, with a **block preconditioner built
on the A-Φ structure**: auxiliary-space or multigrid on the curl-curl block,
multigrid on the Φ block, and the coupling handled approximately.

That is the genuinely differentiated linear-algebra work in this project, and it
is worth more than anything in Stage 2. What Stage 1 contributes to it:

- **an exact reference** to validate the iterative solver against, on problems
  small enough for both
- **equilibration**, which any iterative method wants too
- **an exact solve of small blocks**, usable inside a preconditioner
- **the residual and reporting machinery**, which is the same

What Stage 1 does *not* contribute: the ordering and symbolic work has no role
in an iterative solve. That is expected, and it is a reason to keep Stage 1
scoped tightly rather than gold-plated.

One thing not to repeat from the direct path: a frequency sweep cannot reuse a
factorization, and it cannot reuse an iterative solve either — **every value
changes with frequency** (1557013 of 1557522 on the cylinder). What it can reuse
is the analysis, the pattern, and a preconditioner's *structure*. Anything that
claims to reuse more than that is wrong.

---

## 12. A third-party backend — planned, after Stage 1, never relied on

**Committed 26 Sept: a library backend will be added once the in-house solver
works.** Not instead of it, and not as a dependency: behind a CMake option that
is **off by default**, so the whole suite stays runnable on a clean checkout
with nothing installed. That default is what separates *having* the option from
*relying* on one, and it is not negotiable.

**The choice: MUMPS.** It is the only strong candidate that supports
**complex symmetric indefinite** systems natively, which is exactly what
`RowScaled` and `ScaledPhi` produce, and it brings AMD plus METIS nested
dissection and Bunch-Kaufman pivoting -- the whole of Stage 2 (§10), already
written and tested by people who do this full time. Licensing is CeCILL-C,
which permits linking into closed-source software; the review is in
`LINEAR_SOLVER.md` and still stands. The alternatives and why not:

| | why not |
|---|---|
| SuiteSparse CHOLMOD / UMFPACK | the supernodal modules are GPL, which forces a choice between open-sourcing this solver and legal exposure |
| PARDISO | technically fine, but its commercial terms have changed hands; verify before relying, do not assume |
| SuperLU / SuperLU_DIST | BSD and its static pivoting matches §5's choice, but it has no symmetric-indefinite mode, so the half-memory benefit of `SparseSymmetric` is lost |
| PaStiX | a genuine alternative: complex symmetric, good performance, also CeCILL-C. Heavier build (wants Scotch) and a smaller community. The second choice if MUMPS disappoints |
| Eigen | MPL2, header-only, no Fortran -- by far the easiest to vendor, and the fallback if MUMPS's Fortran/BLAS toolchain proves painful on Windows. But `SimplicialLDLT` does not pivot, so it adds no robustness over Stage 1, and it is simplicial rather than supernodal so it is not much faster either |

**Two decisions, both involving the word MUMPS, and they are not the same one.**
`LINEAR_SOLVER.md`'s superseded recommendation was *link MUMPS **instead of**
writing a solver*. This section is *when a backend is added **alongside** ours,
which library*. The first was reversed; the second is answered MUMPS.

**It may turn out not to be needed at all**, which is why the decision waits for
step 3's number. If AMD puts the direct ceiling near 100-200 k unknowns and the
EDA problems sit under it, the in-house solver covers EDA outright. And the THz
case at 15-25 M DOFs fits in no direct solver's memory, MUMPS included -- that
needs §11's iterative path, which is a different choice of library or none. So
"in-house for EDA, iterative for THz, MUMPS never needed" is a live outcome.
Recording the choice now costs nothing; committing to the integration before
step 3 would be deciding without the fact that decides it.

Deliberately *not* part of Stage 1 or Stage 2. The decision in
§0 stands: nothing in this project's build requires a third-party package. But
having a library available as an **option** is a different thing from relying on
one, and it may be worth adding later.

**What it would be good for, and what it would not.**

As a *correctness* reference it would add little. `solve_dense` is already exact
on small systems, already present, and needs no dependency — and on a small
problem it is a stronger check than a library, because there is nothing to
configure wrongly. The `LDLᵀ`-vs-`LU`, three-conditioning and DC-milestone
checks in §8 are also all in-house.

As a *performance yardstick* and a *large-problem path* it would add real value:
it is how you find out whether Stage 2's supernodal work is worth doing, and it
would give a fast solve before Stage 2 exists.

**What it would cost.** MUMPS needs a Fortran toolchain plus BLAS/LAPACK, and
METIS on top if nested dissection is wanted — on Windows that is a real setup
burden. It would have to sit behind a CMake option that is **off by default**, so
the whole test suite stays runnable on a clean checkout with nothing installed.
That default is what separates *having* the option from *relying* on it, and it
is not negotiable. The other risk is organisational rather than technical: a fast
backend quietly becomes the only path anyone runs, and the in-house solver rots.
The answer to that is §10's — the in-house solver stays the bitwise-deterministic
reference, and Stage 2 is validated against it.

**When to decide.** After step 3 of §9 produces the ordering table. That
measurement says where the direct ceiling actually sits on these meshes, and
therefore whether a library is needed at all. Deciding before it is deciding
without the number.

**Shape, if it happens.** A backend interface both implementations satisfy —
`analyze` / `factorize` / `solve` with the §5 `SolveReport`, so the reporting is
identical either way — plus a comparison harness printing residual, `nnz(L)` and
time side by side. Stage 1 gets built first regardless, because an interface
designed against a working implementation is a better interface than one designed
in the abstract.

---

## 13. Step 1 done, 26 Sept: ordering, and what it can and cannot yet prove

`include/aphi_solver/ordering.hpp`, `src/ordering.cpp`,
`tests/test_ordering.cpp`. `Ordering`, `Permutation`, `compute_ordering`,
`permute_pattern`, `bandwidth_stats`. AMD **throws** rather than falling back to
another ordering, so it cannot report someone else's fill as its own.

RCM is Cuthill-McKee with a George-Liu pseudo-peripheral start, neighbours in
increasing degree with ties by index, components taken in order of their
lowest-numbered vertex — **deterministic**, which is what lets a factorization
built on it be compared bitwise later.

| | unknowns | bandwidth | mean bandwidth |
|---|---|---|---|
| shuffled 200-chain | 200 | 186 → **1** | — |
| `cylinder_box.msh` | 26907 | 26905 → 2753 | 15841 → **2112** |
| `loop_cut.msh` | 24342 | 24341 → 2781 | 14166 → **2100** |

The chain is the test that says whether it works at all: a path graph with
scrambled labels has an optimal bandwidth of 1, and RCM finds it.

Isolated unknowns and multiple components are handled rather than assumed away —
a voltage port's constraint row holds nothing but its diagonal (§10 of
`ASSEMBLY_PLAN`), so the adjacency graph genuinely has isolated vertices.

### The controls, including two that could not bite

| control | outcome |
|---|---|
| permute with `perm` where `iperm` belongs | **caught**, 5 checks |
| start at vertex 0 instead of a pseudo-peripheral one | **caught** — chain bandwidth 1 → 2, profile 200 → 210 |
| visit neighbours highest-degree first | **not an error in effect**, see below |
| skip the reversal (plain Cuthill-McKee) | **not detectable by these metrics**, by construction |

**Highest-degree-first changes almost nothing here.** Measured: mean bandwidth
2035.7 against the correct rule's 2111.5 on the cylinder — marginally *better* —
and 2149.5 against 2099.8 on the loop, marginally worse. ±4 % with no consistent
sign. Cuthill-McKee's degree rule matters on some graphs; on these FEM meshes it
does not, and asserting a direction would be fitting to noise. Recorded as a
measurement rather than papered over with a threshold.

**Skipping the reversal cannot be caught by any bandwidth or envelope metric on
a symmetric pattern**, and this is a theorem rather than a gap in the tests. For
a symmetric pattern the total lower envelope equals the total upper envelope —
each is the sum over rows of the distance to the furthest stored entry on one
side, and symmetry maps one sum onto the other — so reversal merely swaps them
and leaves every such total unchanged. The reversal's benefit is in the
**factor**, where elimination order is not symmetric. So it moves to step 2's
control list, to be re-run against `nnz(L)`.

That is the honest limit of step 1: it proves the permutation machinery is
correct and that RCM reduces what RCM reduces. Whether any of it reduces **fill**
is unmeasurable until the symbolic factorization of step 2 exists, and the
ordering table of §3 is not a table until then.

### A harness note

The first run of these controls reported "permute with `perm` instead of
`iperm`" as caught only by a crash, exit 9009. That was the control harness
reading a stale executable, not a result — re-run, it fails 5 checks cleanly.
Second time a harness artefact has been mistaken for a finding in this project
(`ASSEMBLY_PLAN` §13 has the first), so: a control result that says "crash" and
not "check" is worth re-running before it is believed.

---

## 14. Step 2 done, 26 Sept: the ordering table is a table

`include/aphi_solver/symbolic.hpp`, `src/symbolic.cpp`,
`tests/test_symbolic.cpp`. `elimination_tree`, `factor_row_counts`,
`predict_factor_size`, `analyze`. `SymmetricAdjacency` moved out of
`ordering.cpp` into `ordering.hpp`, so the ordering and the symbolic pass share
one definition of "the graph of the matrix" instead of two.

### §3's table, finally with numbers

| mesh | unknowns | nnz(A) | natural `nnz(L)` | RCM `nnz(L)` | |
|---|---|---|---|---|---|
| `cylinder_box.msh` | 26907 | 939575 | 202862505 — **3869 MB** | 27635653 — 527 MB | **7.34× less** |
| `loop_cut.msh` | 24342 | 763812 | 146936275 — 2802 MB | 25525947 — 486 MB | **5.76× less** |

**The natural ordering would ask for 3.9 GB of complex factor on a 27000-unknown
problem.** That is the entire case for reordering, and it is why the budget is
checked from the counts before the pattern is allocated — `predict_factor_size`
is O(n) memory whatever the fill, so it can report a ruinous ordering rather than
attempt it.

**What the table says about the direct ceiling.** RCM needs 527 MB at 27000
unknowns and about 1030 nonzeros per row of `L`. Fill per row grows with problem
size, so a straight-line reading is optimistic, but even so RCM runs out of a
workstation somewhere around 50–100 k unknowns. AMD typically improves on RCM by
2–5× on 3D problems, which would put the ceiling nearer 100–200 k. **That is the
number step 3 exists to produce**, and it is also the number that decides whether
the §12 backend is needed for EDA-scale work or only for the THz case.

### Verification

The strongest check is an **independent dense reference**: a naive symbolic
elimination that simulates fill by explicit set union in O(n³), sharing no code
with the elimination tree or the reachable-set walk. On 80 random patterns across
both orderings it agrees with `predict_factor_size` on `nnz(L)` **exactly, 80 of
80**. Everything clever in the real implementation — path compression, the
reachable-set walk, the marking — is a chance to be subtly wrong while still
producing a plausible number, and this is what closes that off.

Alongside it, three factors known by hand: diagonal (`nnz(L) = n`, n roots),
tridiagonal (`nnz(L) = 2n − 1`, tree a path), dense (`nnz(L) = n(n+1)/2`). Plus:
`parent[j] > j` or −1; `predicted_nnz == factor.nnz()` from two independent
passes; `L` lower triangular, each row ascending, diagonal last; `L` contains
`tril(A)`; an upper-triangle pattern gives the **same** tree and the same
`nnz(L)` as the full one.

### Controls, and the one step 1 could not test

| control | outcome |
|---|---|
| `row_reach` forgets to mark, so nodes are revisited | caught, 4 checks |
| the diagonal is not counted | caught, 7 checks |
| elimination tree without path compression | caught, 9 checks |
| **skip RCM's reversal** (plain Cuthill-McKee) | **caught — see below** |

§13 deferred the reversal control here because no bandwidth or envelope measure
can see it on a symmetric pattern. `nnz(L)` sees it clearly:

| | cylinder `nnz(L)` | per row | loop `nnz(L)` | per row |
|---|---|---|---|---|
| RCM, with the reversal | 27635653 | 1027 | 25525947 | 1049 |
| plain Cuthill-McKee | 55722203 | 2071 | 49098486 | 2017 |

**The reversal is worth 2.02× in fill** — about 535 MB on the cylinder. George's
reversal is not a cosmetic detail.

And the suite did **not** catch it at first: 67 checks passed with the reversal
removed, because nothing asserted a fill figure. So `mean_row < 1500` is now
pinned on both meshes — correct RCM gives 1027 and 1049, plain CM gives 2071 and
2017, and the threshold sits between with about 1.4× margin either side. A
measured threshold, with what it separates written beside it.

---

## 15. Step 3 done, 26 Sept: AMD, and the direct ceiling

`approximate_minimum_degree` in `src/ordering.cpp`, after Amestoy, Davis & Duff
(1996; Algorithm 837, TOMS 2004). Quotient graph with element absorption and
AMD's approximate external degree, selected through degree buckets.

**Deliberately left out**, per Stage 1's correctness-over-speed brief:
supervariables (indistinguishable-vertex detection by hashing), aggressive
absorption, and the paper's packed integer workspace with garbage collection —
replaced by per-vertex vectors, which cost memory and locality and are far
easier to be sure of. All three are speed, not quality.

**Why AMD was safe to write at this stage:** an ordering cannot be *wrong*, only
worse. Any permutation gives a correct factorization, so a defect shows up as
fill — and step 2 measures fill exactly. That is a rare luxury and it is why
this piece came before the numeric factorization rather than after.

### §3's table, complete

| ordering | cylinder `nnz(L)` | per row | loop `nnz(L)` | per row |
|---|---|---|---|---|
| natural | 202862505 — 3869 MB | 7539 | 146936275 — 2802 MB | 6036 |
| RCM | 27635653 — 527 MB | 1027 | 25525947 — 486 MB | 1049 |
| **AMD** | **12672169 — 241 MB** | **471** | **8668644 — 165 MB** | **356** |

**AMD is 16× better than no ordering and 2.2–2.9× better than RCM**, and it
orders in 230 ms / 161 ms. Quality is confirmed against an **exact** minimum
degree computed independently in the tests: worst `nnz(L)` ratio **1.02** over
25 random patterns, i.e. within 2 % of exact MD despite the approximate degree.

### The direct ceiling — the number this step existed to produce

At 27000 unknowns AMD needs **241 MB** and about 470 nonzeros per row of `L`.
Fill per row grows with problem size, so the following is an **estimate, not a
measurement**, and it is the first thing to re-measure on a bigger mesh:

| unknowns | estimated per row | estimated factor |
|---|---|---|
| 27 k (measured) | 471 | 241 MB |
| 100 k | ~700 | ~1.4 GB |
| 200 k | ~900 | ~3.6 GB |
| 500 k | ~1200 | ~12 GB |

So **the in-house direct solver should reach roughly 200–500 k unknowns** on a
workstation. Two consequences:

- **EDA-scale work is covered in-house.** If the problems of interest sit under
  a few hundred thousand unknowns, Stage 1 plus AMD is the whole answer and the
  §12 backend is not needed for them.
- **The THz case is not close.** 15–25 M DOFs is two orders beyond this, and no
  direct solver — MUMPS included — changes that. §11's iterative path is the
  only route there, which is what §12 already says.

### Controls: four run, and two of them were about the code, not the tests

| control | outcome |
|---|---|
| drop AMD's approximate degree, keep only the growth bound | **caught** — fill 12.7 M → 60.3 M, *worse than RCM* |
| elements never absorbed | **caught, by a time guard** — see below |
| neighbour lists not pruned of the new element's members | **caught, by a pinned figure** — see below |
| `min_degree` not reset after pruning | **not a defect** — the line was redundant and has been removed |

**Element absorption is worth 112× in time and nothing at all in fill.** Without
it, `nnz(L)` comes out *bit-identical* and the ordering takes **25567 ms instead
of 228**. That is exactly what the paper says absorption is for, and it means no
fill-based check can ever see it. So there is now a generous time bound —
under 5 s, against 230 ms measured and 26 s broken.

**Not pruning a neighbour's variable list costs 11 % of fill** — 471 → 523 per
row, about 27 MB here — and on small random patterns the same defect costs only
3 %, so the exact-MD ratio test does not see it either. Pinned at
`mean_row < 500` on both meshes, the same way RCM's figure is pinned in §14.

**`min_degree = 0` after the update loop was dead code I wrote.** Removing it
changed neither the fill nor the time, because the `std::min` inside the loop is
already the only place the minimum can fall — a variable not adjacent to the
pivot keeps its degree, and the scan had already established nothing was lower.
The control did not expose a missing test; it exposed a redundant line, which is
now gone with the reasoning in a comment.

---

## 16. Step 4 done, 26 Sept: LDLᵀ, and the first measured sight of the conditioning problem

`include/aphi_solver/factorization.hpp`, `src/factorization.cpp`,
`tests/test_factorization.cpp`. `permuted_lower`, `factorize_ldlt`,
`reconstruct`, plus `pattern_of` and `analyze`/`compute_ordering` overloads that
take a matrix instead of a pattern. `col_count` added to `SolverAnalysis`,
because an up-looking factorization lays out its storage by **column** while the
symbolic pass counts by row — different distributions, both summing to
`predicted_nnz`.

Up-looking and scalar, per Davis §4.8. `P A Pᵀ = L D Lᵀ` with `L` unit lower
triangular, its diagonal not stored. **No pivoting yet**: a pivot at or below the
floor stops the factorization and is reported, rather than being divided by.

### Verification: reconstruction

The test that matters forms `L D Lᵀ` again and compares against `P A Pᵀ` entry by
entry. **Worst relative error 1.4e-16 over 60 random cases** across two
orderings — machine precision. Every clever part of the algorithm (the tree walk,
the topological order, the column bookkeeping) is checked by that one number,
and all four negative controls were caught by it:

| control | worst reconstruction error |
|---|---|
| conjugate, as a Hermitian factorization would | 0.61 |
| process the row pattern in reverse topological order | 0.16 |
| `y[j]` not cleared after use | 0.63 |
| `permuted_lower` does not normalise below the diagonal | caught, 7 checks |

**One case exists solely to catch conjugation.** These matrices are complex
*symmetric* and not Hermitian, so a conjugating implementation still produces a
plausible `L` and `D`; only reconstruction against genuinely complex
off-diagonals exposes it. That control's error is 0.61 — not subtle once
measured, and completely invisible without it.

Also pinned: a diagonal matrix gives `L` empty and `D` the diagonal untouched; a
2×2 by hand (`L(1,0) = 3/2`, `D = [2, 1/2]`); the factor holds *exactly* the
nonzeros the symbolic pass predicted; `L`'s columns ascend strictly below the
diagonal; `permuted_lower` equals `A(perm[row], perm[col])` for every entry, with
all of `A` accounted for; and `[[0,1],[1,0]]` — the smallest complex symmetric
matrix that breaks an unpivoted `LDLᵀ` — is refused at column 1 rather than
divided by.

### The real matrix: it factorizes, and the pivots are alarming

`meshes/loop_cut.msh` at 1 MHz, `conditioning = row_scaled`:

| | |
|---|---|
| unknowns / stored nonzeros | 37064 / 796053 |
| `analyze` | 975 ms, `nnz(L)` = 28070495 — 535 MB |
| `factorize_ldlt` | **44955 ms**, all 37064 columns |
| **smallest \|D\|** | **8.1e-16** |
| **largest \|D\|** | **1.09e11** |
| **ratio** | **1.3e26** |

Two things to take from this.

**The factorization completed** — no zero pivot, `nnz` exactly as predicted. So
the machinery works on a real problem, not only on test fixtures.

**But a pivot ratio of 1.3e26 is ten orders beyond what double precision can
carry.** This is the first *measured* sight of the ill-conditioning that
`docs/CONDITIONING.md` and the whole formulation discussion were about, and it
arrives exactly where predicted: the `A` block sits near 1e11 while `RowScaled`
divides the Φ block by `jω` down to order 1, so the two blocks differ by eleven
orders and elimination lands on the difference. Note the loop has **no voltage
port**, so this is not the unit-diagonal constraint row of §10 — it is the
formulation's own scaling.

That is what steps 6 and 7 exist for, and the order is now clearly right:
equilibration first (`diag(d) A diag(d)`, which exists already and is exactly
aimed at this), then static pivoting with a floor, then refinement to recover
what the perturbation costs. **A residual is the only thing that can say whether
any of it worked**, and that is step 5.

**45 s is too slow for the test suite**, so the real-mesh factorization stays a
probe for now; step 5 turns it into a committed tool, since every step from here
needs it. Scalar and unblocked is Stage 1's brief, and §10 already names
supernodal blocking as the 5–20× that answers it.

---

## 17. Step 5 done, 26 Sept: the first solve

> **CORRECTION, added with Sec. 18.** This section originally concluded "at
> 100 MHz the solution is 6.9 % wrong. That is not an answer." **That was
> wrong**, and the fault was the measure, not the solve. The 6.9 % is
> `||Ax-b|| / ||b||` where `||b||` is 1.6e-9 with a single nonzero entry, against
> matrix entries of 4.4e12. The standard backward error for the same solve is
> 4.2e-25. The solve was accurate all along. Everything below about the residual
> being the only thing that can detect a bad answer still holds -- it is which
> residual that was wrong. Sec. 18 has the measurement and the fix.

`solve_in_place`, `solve`, `relative_residual` in `factorization.cpp`, and
`tools/solve_mesh.cpp` — a committed tool, because the suite cannot factorize a
real mesh in reasonable time and every step from here needs somewhere to watch
these numbers move.

Three sweeps: forward through `L` (unit lower, nothing to divide by), then `D`,
then backward through **`Lᵀ`** — the transpose, never the conjugate transpose.
The permutation is applied with `perm` **in both directions**: `b̃[i] = b[perm[i]]`
going in, `x[perm[i]] = x̃[i]` coming out, because the factorization is of
`P A Pᵀ`.

### Verification: three independent routes

1. **The residual** `‖Ax − b‖/‖b‖` on the original system.
2. **Three orderings must agree.** The solution cannot depend on the order the
   unknowns were eliminated in, so natural, RCM and AMD must produce the same
   `x`. This is what catches a `perm`/`iperm` confusion, which otherwise gives a
   self-consistent-looking field that solves nothing.
3. **Against `solve_dense`**, which shares no code with any of this.

Plus solves known without a solver (a diagonal, and `[[2,3],[3,5]]x = [5,8] →
[1,1]`), and a check of `Ax = b` against `A`'s **own entries** rather than through
`matvec`, so no shared code could absorb a conjugation on both sides.

**The residual is itself checked for being a real measurement**: perturbing one
entry of a correct solution by 0.01 must make it jump, and `x = 0` must give
exactly 1. A residual that cannot report a bad answer proves nothing.

All five controls caught:

| control | residual it produced |
|---|---|
| conjugate in the backward sweep (`Lᴴ` not `Lᵀ`) | 0.017 |
| permute the right-hand side with `iperm` | 0.50 |
| un-permute the solution with `iperm` | 0.38 |
| skip the division by `D` | 99.1 |
| backward sweep runs forwards | 2.8e-4 |

### The real mesh: it solves, and the answer is wrong

`examples/loop_sweep.aphi`, `row_scaled`, 37064 unknowns, AMD:

| frequency | factor | solve | **residual** | min \|D\| | max \|D\| |
|---|---|---|---|---|---|
| 10 kHz | 48153 ms | 97 ms | **1.8e-08** | 8.1e-16 | 1.08e11 |
| 100 MHz | 47394 ms | 92 ms | **0.069** | 8.1e-16 | 1.98e12 |

**At 100 MHz the solution is 6.9 % wrong. That is not an answer.** At 10 kHz,
1.8e-08 is eight orders worse than the 1e-16 the same code achieves on a
well-scaled test matrix.

This is the single most useful measurement so far, for three reasons.

**It vindicates reporting the residual on every solve** (§5). The factorization
returned `ok`, filled exactly the predicted nonzeros, and produced a smooth
field. Nothing but the residual distinguishes that from a correct answer.

**It confirms the conditioning problem is real and quantitative**, not a
theoretical worry: the pivot range spans 8.1e-16 to 1.98e12, a ratio of 2.4e27,
against the ~1e16 double precision carries. §16 saw this in the pivots; step 5
shows what it costs in the answer.

**It sets the acceptance criterion for steps 6 and 7.** They are no longer
speculative hardening — equilibration and static pivoting with refinement have a
number to move, and this table is what they have to beat. If they do not bring
the residual to ~1e-12 or better across the sweep, they have not worked.

Two things this does **not** show. The 45–48 s factorization is Stage 1's
declared scalar cost, not a surprise, and §10 names the fix. And nothing is
extracted from these solutions yet — no currents, no `R`, no `L` — so the
physics remains unvalidated. That is step 9, and it cannot usefully run until the
residual is small.

---

## 18. Step 6 done, 26 Sept: equilibration works, and §17's verdict was wrong

`compute_symmetric_equilibration` and `apply_symmetric_equilibration` gained
`SparseSymmetric` overloads; `solve_symmetric` in `factorization.cpp` ties the
whole pipeline together — **scale, order, factorize, solve, un-scale, measure** —
and `tools/solve_mesh.cpp` now runs every frequency **twice**, with and without
scaling, because a number only ever seen scaled cannot say whether the scaling
helped.

The scaling vector cannot be computed from the stored triangle as if it were a
general matrix: row `i`'s largest entry includes the mirrors of entries held in
earlier rows, so every stored entry contributes to **both** its indices. One test
checks the symmetric overload against the full-storage one, which sees both
triangles and cannot make that mistake.

### What equilibration did

`examples/loop_sweep.aphi`, `row_scaled`, 37064 unknowns, AMD:

| | before | after |
|---|---|---|
| pivot range | 8.1e-16 … 1.98e12 | **0.015 … 1.0** (10 kHz), 6.7e-5 … 1.0 (100 MHz) |
| pivot ratio | **2.4e27** | **66** and **1.5e4** |
| max \|L\| | — | **1.95** and **14.8** |
| scaling spread | — | 9.8e12 and 5.6e13 |

**Equilibration worked, and completely.** The pivot ratio fell by twenty-three
orders of magnitude, and `max |L|` — the multiplier growth that actually measures
whether an unpivoted factorization was stable — is **under 15**. A factorization
with those numbers is not merely adequate, it is clean.

### And that is what exposed the error in §17

The `||b||`-relative residual barely moved: 0.0685 → 0.0668. Taken with a pivot
ratio of 1.5e4 and `max |L|` of 15, that is a contradiction — such a
factorization cannot produce a 6.7 % error. So the measure was wrong, and
measuring the norms settled it:

| | |
|---|---|
| ‖b‖ | 1.59e-9, with **exactly 1 nonzero entry of 37064** |
| ‖x‖ | 57.8 |
| max \|A\| | 4.42e12 |
| ‖Ax − b‖ | 1.06e-10, absolute |
| ‖Ax−b‖ / ‖b‖ | **0.0668** — what §17 reported |
| ‖Ax−b‖ / (max\|A\|·‖x‖ + ‖b‖) | **4.2e-25** — the standard backward error |

`b` has one nonzero because `row_scaled` puts `r·I = I/(jω)` in the port row, and
at 100 MHz that is 1.6e-9. Dividing a rounding-level residual by it inflates
noise into a percentage. **The solve was accurate to 4e-25 all along.**

Three things follow, and the third is the uncomfortable one.

**`backward_error` is now the number to judge a solve by**, reported alongside
the `||b||`-relative one, and it is what `solve_mesh` thresholds against. Note
what the old measure also was: **not comparable between formulations of the same
problem**, since `||b||` scales with the conditioning choice — `Natural` would put
`I = 1` in that row where `RowScaled` puts `I/(jω)`.

**§17's acceptance criterion for steps 6 and 7 is void.** There was no 0.069 to
fix. Step 7's static pivoting and refinement are now insurance against matrices
that *are* hard, not a repair for this one — and their justification has to come
from a case that genuinely fails, not from this measurement.

**I stated a wrong conclusion confidently, in a commit message and a plan
section, on the strength of a metric I had written myself and not questioned.**
The tests around it were all passing and all correct; they simply never compared
the metric against a second opinion. The general form of this is the same failure
this project keeps finding — a check derived from the same assumption as the thing
it checks — and it is worth noting that it caught *me* rather than the code.

Two attempts to reproduce the discrepancy in a small fixture **failed**, and the
test file says so rather than pretending otherwise: large matrix entries alone do
not do it, because `||x||` then shrinks with `||b||`; nor did a hand-made
near-singular block. The separation needs `||x||` to stay large while `||b||` is
tiny, which is a property of that assembled system. The unit test therefore
asserts the invariant that always holds — `backward_error ≤ relative_residual`,
same numerator and a larger denominator — plus that neither measure can be fooled
by a wrong answer, and the real evidence stays where it was measured.

### On the fixture where equilibration barely helps

The badly-scaled 24×24 test matrix shows 6.3e-11 unscaled against 5.0e-11
equilibrated. That is not a defect either: the relative residual is a *backward*
error, and on a matrix spanning twelve orders it already sits near what the
conditioning allows. What equilibration buys is accuracy in `x`, which a backward
error does not see. The threshold there is set from the measurement, with the
reason written beside it.
