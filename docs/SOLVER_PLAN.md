# Sparse linear solver — plan

Phase 05. The matrix exists and is cheap to rebuild (`docs/ASSEMBLY_PLAN.md`);
nothing solves it yet, so nothing in this project has been validated
numerically. Everything to date is structural or algebraic.

**Decided 26 Sept: in-house, no third-party library.** This reverses the
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
