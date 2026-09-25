# Element matrices and global assembly — plan

*Opus 5, 25 Sept 2026. Plan for review; nothing here is implemented yet.
`docs/ROADMAP.md` Phase 04 steps 3 and 5, and the step after the DOF map in
`docs/INPUT_FILE_PLAN.md` §6.*

**This replaces an earlier draft** that returned a 4.7 KB `ElementMatrices`
struct per tetrahedron and left assembly for later. That design was wrong in
three ways, all pointed out in review:

- it **materialised** every element matrix, when each one is wanted once and
  then thrown away;
- it computed `C_PhiA` independently when `C_PhiA = jω·C_APhiᵀ` **exactly**,
  so one of the two couplings is free;
- it deferred assembly, which hides the decision that actually governs the
  cost — **when the global sparse structure is built**.

The revision keeps element matrices on the stack, scatters them immediately,
and builds the sparsity pattern in a separate pass beforehand.

---

## 1. Shape of the thing

```
  symbolic pass          numeric pass
  -------------          ------------
  DofMap + Mesh          for each body:            <- material coefficients
        |                  for each of its tets:      once per body
        v                    compute tet geometry
  SparsityPattern            kernel_AA    -> 6x6  \
  (row_ptr, col_idx)         kernel_APhi  -> 6x10  |  stack buffers,
        |                    kernel_PhiPhi-> 10x10 /   ~2.6 KB total
        v                    scatter all four blocks into the values array
  Sparse<Complex>            (PhiA is the transpose of APhi, times j*omega)
  values all zero
```

Two passes, because the pattern depends only on connectivity and the values
only on materials and frequency. A frequency sweep therefore builds the
pattern **once** and refills the values 41 times — which the single-pass
triplet approach cannot do at all.

---

## 2. Kernels

With `ν = 1/μ`, `α = jωσ − ω²ε`, `β = σ + jωε` (`FORMULATION.md` §1.4):

| kernel | out | integrand |
|---|---|---|
| `kernel_AA` | 6×6 | `ν (curl Wᵢ)·(curl Wⱼ) + α Wᵢ·Wⱼ` |
| `kernel_APhi` | 6×10 | `β Wᵢ·∇S_b` |
| `kernel_PhiPhi` | 10×10 | `β ∇S_a·∇S_b` |
| *(PhiA)* | 10×6 | **not computed** — it is `jω · APhiᵀ` |

**Three kernels, not five.** The curl-curl and mass terms are both A–A and
*add* in the assembled matrix, so computing them together means one block
and one scatter instead of two. And `α = jω·β` makes `C_PhiA = jω·C_APhiᵀ`
an identity, not a coincidence — it is where the system's asymmetry comes
from, and what `CONDITIONING.md` Formulation 1 removes by dividing the Φ
rows by `jω`. Under that scaling the two coupling blocks become exact
transposes and the matrix is complex-symmetric.

### All three formulations, from two scalars

All three of `CONDITIONING.md`'s formulations are assembled, selected by a
flag, per `ROADMAP.md` Phase 04 step 6 ("keep all forms available since
later phases compare them"). They differ only in two numbers.

Write the **coefficient-free** element integrals as

```
K = ∫ (curl W)·(curl W)     M = ∫ W·W     C = ∫ W·∇S     L = ∫ ∇S·∇S
```

and introduce a **row scale `r`** (multiplying every Φ-row contribution and
the Φ/port right-hand side) and a **column scale `c`** (the unknown
substitution `Φ = c·Φ'`). Then every formulation is:

| block | assembled as |
|---|---|
| (A,A) | `νK + αM` |
| (A,Φ) | `c · βC` |
| (Φ,A) | `r · αCᵀ` |
| (Φ,Φ) | `r·c · βL` |
| Φ and port RHS | `r ·` (port currents) |
| prescribed Φ or V | `÷ c` before use |
| recovered Φ or V | `× c` after the solve |

with

| | `r` | `c` | (A,Φ) | (Φ,A) | symmetric |
|---|---|---|---|---|---|
| **F3 — natural** | `1` | `1` | `βC` | `αCᵀ` | no |
| **F1 — row scaling** | `1/(jω)` | `1` | `βC` | `βCᵀ` | **yes** |
| **F2 — scaled Φ** | `1` | `jω` | `αC` | `αCᵀ` | **yes** |

Both symmetric cases work because `α = jω·β`. In F1 the `jω` on the Φ row
is divided out; in F2 it is instead supplied to the Φ *column*. Verified
against `CONDITIONING.md`: F1 turns `(jωσ − ω²ε)` into `(σ + jωε)`; F2 gives
`K_APhi_new = jω·K_APhi` and `K_PhiPhi_new = jω·K_PhiPhi` with `K_AA`,
`K_PhiA` and both right-hand sides untouched. Both match.

**The symmetry condition is exactly `c == r·jω`** — one runtime assertion
rather than a claim in a comment.

The kernels are untouched by the choice; it only alters what multiplies a
block on its way into the matrix. And **the sparsity pattern is identical
for all three**, so it is still built once and shared between modes and
across a sweep.

#### Ports, which F1 and F2 treat differently

`V_k` is Φ-like — its row *is* the port-current equation — so it takes both
scales. `PORT_AND_DOF_PLAN.md` §8 records this; getting it wrong makes the
read-back current wrong by exactly `jω`, which reads as a physics error
rather than a bookkeeping one.

| | F1 | F2 |
|---|---|---|
| unknown | `V_k` unchanged | `V'_k = V_k/(jω)` |
| the port row | ÷ `jω` | unchanged |
| current port RHS | `I_k/(jω)` | `I_k` |
| voltage port: prescribe | `V_k` | `V_k/(jω)` |
| read back `I_k` | `jω ×` residual | residual |
| after the solve | nothing | `V_k = jω·V'_k`, `Φ = jω·Φ'` |

#### Three more consequences, each needing a test

1. **F1 requires ω > 0**; `1/(jω)` is undefined at DC. **F2 requires it
   too**, but for the prescribed values rather than the matrix: a voltage
   port would need `V_k/(jω)`. Both must *refuse* DC, not divide by zero.
   No loss — at ω = 0 the system decouples (`FORMULATION.md` §2) and there
   is no coupling block to symmetrise. F3 works at any frequency.
2. **They degenerate at low frequency in opposite directions**, which is
   why keeping all three is worth the flag rather than picking one: F1's Φ
   block is `βL/(jω)` = O(σ/ω) and blows up; F2's is `αL` = O(ω) and
   vanishes; F3's is `βL` = O(σ) and stays put. (`CONDITIONING.md`'s
   Formulation 2 paragraph had this backwards and was corrected 25 Sept;
   the derivation is recorded there.)
3. **Symmetry needs gauge Method A.** Per `CONDITIONING.md`'s own table,
   the Munteanu projection reintroduces asymmetry regardless, so Method D +
   F1 or F2 pays the DC-degeneracy cost for no symmetry at all. Reject that
   combination rather than produce a non-symmetric matrix while claiming
   otherwise.

**Deferred: storing only one triangle.** A complex-symmetric matrix needs
half the storage, and MUMPS accepts that form. Doing it now would mean
guessing the solver's expected layout before Phase 05 has chosen one, and
it complicates the cross-check in §7.4. Assemble full storage in all three
modes; revisit when the solver is actually linked.

Signature — **the caller owns the memory**, which is what keeps element
matrices off the heap and makes the kernels testable in isolation:

```cpp
void kernel_AA     (const TetGeometry& g, const ElementCoefficients& c,
                    std::complex<double>* out);   //  36 =  6x6
void kernel_APhi   (const TetGeometry& g, const ElementCoefficients& c,
                    std::complex<double>* out);   //  60 =  6x10
void kernel_PhiPhi (const TetGeometry& g, const ElementCoefficients& c,
                    std::complex<double>* out);   // 100 = 10x10
```

**All three are pure functions of the geometry** -- no `Mesh`, no tet index,
no global state. They use the LOCAL Whitney basis, and the global edge
orientation is applied by the scatter through `DofMap`.s `coeff` (§9.1).
That is exact, because `W_global = s * W_local` gives
`A_global[p][q] = s_p * s_q * A_local[p][q]`, which is what
`cp * cq * block[p][q]` computes; for the (A,Phi) block `grad(S)` carries no
orientation, so `cq = 1` and only `s_p` applies.

Being pure is what makes them testable against hand-computed numbers with a
`TetGeometry` built by hand and no mesh at all.

### Quadrature

A degree-2 (4-point) rule is **exact** for all of these, not an
approximation:

| integrand | degree |
|---|---|
| `(curl W)·(curl W)` | **0** — `curl W = 2 ∇Lᵢ × ∇Lⱼ`, constant per tet |
| `W·W`, `W·∇S`, `∇S·∇S` | 2 — `W` and `∇(P2)` are linear |

Region-wise-constant materials add no degree. The curl-curl part of
`kernel_AA` needs no quadrature loop at all: it is `ν·V·(curl·curl)`.

---

## 3. The symbolic pass

```cpp
struct SparsityPattern {
    int rows = 0;
    std::vector<int> row_ptr;    // size rows + 1
    std::vector<int> col_index;  // ascending within each row
    std::size_t nnz() const { return col_index.size(); }
};

SparsityPattern build_sparsity(const DofMap& dofs, const BoundProblem& bound,
                               const Mesh& mesh);
```

For each tet, gather its live global indices from `DofMap::local_dofs` (at
most 16, fewer once Dirichlet edges and absent Φ nodes drop out) and mark
every ordered pair. Build by counting per row, prefix-summing, filling, then
sorting each row — no `std::set`, no per-row `std::vector`.

**Why this is worth a separate pass.** The alternative, `Sparse<T>::add()`
into a triplet list followed by `compress()`, is what the existing type
does. On the cylinder that means roughly **2.3 M triplets (~55 MB) and a
sort of all of them**, against a CSR of ~1.5 M nonzeros (~30 MB) with no
sort. It also cannot be reused across a frequency sweep.

`Sparse<T>` needs one addition to accept this:

```cpp
static Sparse<T> from_pattern(const SparsityPattern& p, int cols);
int find_slot(int row, int col) const;   // binary search within the row
```

The existing triplet path stays — it is right for the small hand-built
matrices the tests use.

---

## 4. The numeric pass

```cpp
void assemble(const BoundProblem& bound, const Mesh& mesh, const DofMap& dofs,
              const SparsityPattern& pattern, double omega,
              Sparse<Complex>& matrix, std::vector<Complex>& rhs);
```

```
for each body                       <- coefficients computed ONCE per body
    for each tet of that body
        g = compute_tet_geometry(mesh, tet)
        d = dofs.local_dofs(tet, mesh, bound)
        kernel_AA     -> aa[36]     \
        kernel_APhi   -> ap[60]      |  stack, ~2.6 KB, reused every tet
        kernel_PhiPhi -> pp[100]    /
        scatter
```

Looping bodies-then-tets rather than tets-directly is what makes the
material coefficients a per-body cost instead of a per-tet one. Bodies
partition the tets, so every tet is still visited exactly once, and
`BoundBody::tets` is ascending so the mesh access stays ordered.

### The scatter

For local DOFs `p`, `q` with global indices `i`, `j` and coefficients
`cᵖ`, `cᑫ`:

```
matrix[i][j] += cp * cq * block[p][q]
```

A **prescribed** DOF (Dirichlet edge, grounded cut side, pinned node, fixed
port voltage) has no column, so its contribution moves to the right-hand
side instead: `rhs[i] -= cp * block[p][q] * value(q)`. That is the one place
`TetDofs::phi_fixed` is used.

The four blocks are scattered as: `aa` into (A,A), `ap` into (A,Φ), `jω·apᵀ`
into (Φ,A), `pp` into (Φ,Φ).

---

## 5. Parallelism

**The per-tet work is embarrassingly parallel; the scatter is not** — two
tets sharing an edge write to the same row. `ENGINEERING_STANDARDS.md`
already records this.

**Proposed order, and I want to measure before choosing:**

1. **Single-threaded first, and time it.** The kernels are ~300 flops per
   tet; on 16 040 tets I expect assembly in the tens of milliseconds. If it
   is, threading it is not where the time goes and the DOF map/binding
   measurements said the same about their own hot spots. *Measure, do not
   assume* — this project has twice now found the cost somewhere other than
   the obvious place.
2. **Then, if it matters**, the two viable schemes:

| scheme | how | cost |
|---|---|---|
| **atomic adds** | `#pragma omp parallel for` over tets, atomic `+=` on each value slot | simple; contention on shared rows. `std::complex` is not atomic, so the values array must be viewed as `double[2·nnz]` and updated with two atomics — a real wrinkle, not a detail |
| **graph colouring** | colour tets so no two of a colour share a DOF; a parallel loop per colour needs no atomics at all | a colouring pass, and fewer tets per parallel region |

MSVC's `/openmp` is OpenMP 2.0, which does support `#pragma omp atomic` on
scalar `+=`. `/openmp:llvm` gives a newer runtime if needed.

**What to build now so either stays possible:** keep the kernels free of
shared state (they already are — caller-owned output buffers), keep the
per-tet geometry local, and keep the scatter a separate inlineable step. No
OpenMP pragmas until there is a measurement saying where they would go.

### One optimisation deliberately deferred

Precomputing each local pair's position in the values array (`≈ 2.3 M ints,
~9 MB`) would remove the binary search from the inner scatter. It is the
obvious next step *if* the scatter turns out to dominate — but a search in a
row of ~60 entries is ~6 comparisons, and guessing that this matters is
exactly the mistake §5.1 warns against.

---

## 6. Memory, with numbers

**Measured on the cylinder**, 25 Sept, once `build_sparsity` existed --
26 907 unknowns, 16 040 tets, **941 665 nonzeros**, 35.0 per row:

| | |
|---|---|
| CSR `col_index` | 3.6 MB |
| CSR `values` (complex) | 14.4 MB |
| **total** | **17.96 MB** |
| element matrices | **~2.6 KB, on the stack, reused** |
| *(triplet alternative)* | *~55 MB plus a 2.3 M-element sort* |

This replaces an estimate of ~1.5 M nonzeros and ~30 MB, which was high by
about 60 %: the guess assumed 50-70 entries per row and the real figure is
35. Worth recording that the guess was wrong in the safe direction, but
wrong.

---

## 7. Testing

**The kernels are the testable unit**, which is the other reason they take a
caller-owned buffer: a test calls them with a local array and checks numbers,
with no assembly involved.

### 7.1 Hand-computed, on the unit tet

`(0,0,0), (1,0,0), (0,1,0), (0,0,1)`: `V = 1/6`, `∇L₀ = (−1,−1,−1)`,
`∇L₁ = (1,0,0)`, `∇L₂ = (0,1,0)`, `∇L₃ = (0,0,1)`.

```
curl W₀ = 2 ∇L₀ × ∇L₁ = 2 (−1,−1,−1) × (1,0,0) = (0, −2, 2)
K_AA[0][0] (curl part) = ν V |curl W₀|² = ν (1/6)(8) = 4ν/3
```

The whole curl-curl 6×6 is hand-derivable this way and goes in the test as
literals — never as a re-implementation of the code under test.

### 7.2 Properties true on *any* tet

| property | why | a failure means |
|---|---|---|
| `kernel_AA`, `kernel_PhiPhi` symmetric | real basis, symmetric form | an index transposed |
| `K_AA·(G p) = 0` for any nodal `p` | `curl(grad) = 0`, and the discrete gradient of a P1 field lies exactly in the Whitney space | the curl or the sign convention is wrong |
| `rank(curl-curl part) = 3` | 6 edges minus a 3-dimensional gradient subspace | same |
| `kernel_PhiPhi` row and column sums `= 0` | `Σ S_a = 1` so `Σ ∇S_a = 0` | the P2 gradients are not a partition of unity |
| mass part positive definite for real `α > 0` | it is a mass matrix | a sign or a weight |
| at `ω = 0`: `α = 0`, so no mass and no Φ–A coupling | the system decouples (`FORMULATION.md` §2) | the DC path is not what §2 says it is |
| scaling by `s`: curl-curl `→ 1/s`, mass `→ s`, coupling `→ s`, `∇S·∇S → s` | dimensions: `W ~ 1/L`, `curl W ~ 1/L²`, `∇S ~ 1/L`, `dV ~ L³` | a volume factor dropped or applied twice |

The curl-curl and mass terms scale in **opposite** directions (this table
first had both wrong, corrected 25 Sept when the test was written). Refining
a mesh grows the curl-curl term and shrinks the mass term, so their ratio
goes as `1/s²` -- which is the low-frequency conditioning problem in
miniature, showing up in a single element.

### 7.3 The three formulations

- **F1 and F2 really are symmetric**: `M == Mᵀ` entry for entry after
  assembly, not merely in the two coupling blocks. F3 must NOT be -- a
  negative control, since a bug that symmetrised everything would otherwise
  look like success.
- **`c == r·jω` holds** for F1 and F2 and fails for F3. Asserted at run
  time, not just claimed.
- **All three describe the same problem.** Each is the natural system with
  its Φ rows scaled by `r` and Φ columns by `c`, so the exact relation is
  checkable without a solver: taking the F1 or F2 matrix, multiplying every
  Φ and port row by `1/r` and every Φ and port column by `1/c`, must
  reproduce F3 entry for entry, RHS included.
- **F1 and F2 refuse DC**, rather than dividing by zero. F3 does not.
- **F1 and F2 refuse gauge Method D**, since its projection would break the
  symmetry the caller asked for.
- **The pattern is formulation-independent**: `build_sparsity` output is
  identical for all three, which is what lets a sweep share it.
- **Port handling**: a voltage port under F2 must have its prescribed value
  divided by `jω`, and its recovered `V_k` multiplied back. A round trip
  through prescribe-then-recover must return the value written in the input
  file.

### 7.4 Assembly-level

- **Every pattern slot the scatter writes to exists.** `find_slot` returning
  −1 during assembly is a symbolic/numeric mismatch and must abort, not
  silently drop a term.
- **Symmetry of the assembled blocks**: `K_AA` and `K_PhiPhi` symmetric,
  and `C_PhiA = jω·C_APhiᵀ` after assembly, not just per element.
- **Against the existing triplet path**: assemble a small mesh both ways —
  `from_pattern` + scatter, and `add()` + `compress()` — and require them
  equal entry for entry. Two independent routes to one answer.
- **Row count sanity**: a row belonging to an absent or prescribed DOF must
  not exist at all.

### 7.5 Negative controls

Each must be shown to fail: drop the global edge sign; transpose one block;
use a 1-point rule for the quadratic terms; drop the volume factor; skip the
`jω` on the F3 PhiA scatter; forget to scale the port rows under F1; forget
to divide a prescribed voltage by `jω` under F2; swap `r` and `c`.

---

## 8. Implementation order

1. **Quadrature rule**, tested by integrating degree 0/1/2 monomials over
   the reference tet against exact values — before any physics uses it.
2. **`ElementCoefficients`**, with the `ω = 0` degeneracy asserted.
3. **`kernel_AA`** — the constant curl-curl part first, hand-checked, then
   the mass part.
4. **`kernel_PhiPhi`**, **`kernel_APhi`**.
5. **`SparsityPattern` + `Sparse<T>::from_pattern`**, tested on a small mesh
   against the triplet path.
6. **The scatter and RHS**, including prescribed-DOF elimination, in the
   unsymmetric mode first.
7. **F1 and F2**, which are two scalars and their consequences (§2),
   each checked against F3 by the exact row/column-scaling identity in
   §7.3, plus the port round trip.
8. **Time it**, then decide about threading with a number in hand.

---

## 9. Decisions to confirm

1. **Where the edge orientation sign is applied. Decided 25 Sept: in the
   DOF map, via `coeff`; the kernels use the LOCAL basis.**

   An earlier draft of this plan recommended the opposite -- kernels using
   the `_global` accessors, with `DofMap`.s `coeff` reduced to 1. Writing
   out the kernel signatures changed that: taking the sign out of the
   kernels makes all three **pure functions of `TetGeometry`**, with no
   `Mesh` argument and no tet index, which is what lets them be tested
   against hand-computed numbers without a mesh existing at all. The sign
   also then stays in the one place that already implements and tests it
   (`test_local_dof_map` asserts `coeff == tet_edge_signs`), and no
   committed code has to change -- nothing outside tests uses the `_global`
   accessors yet.

   What must NOT happen is both, which would apply the sign twice and
   cancel it, and which no current test would catch. One test asserts it is
   applied exactly once.
2. **Complex throughout, including DC.** *Recommend yes*, per
   `FORMULATION.md` §5.2 — one assembly and solve path, at a factor of two
   in arithmetic on the DC milestone.
3. **Assemble the full matrix, or only the gauge-reduced one?** The DOF map
   already omits tree and Dirichlet edges, so what is assembled *is* the
   reduced system. *Recommend that*, with `principal_submatrix` retained as
   the independent check: assembling ungauged then reducing must equal
   assembling reduced directly.
4. **Frequency sweep**: build the pattern once and refill the values per
   frequency. *Recommend yes* — it is most of the reason for the two-pass
   split.
5. **All three formulations, from one flag** (§2). *Agreed 25 Sept.* The
   default for a solve is the open question, not whether to have all three:
   F1 and F2 buy `LDLᵀ` over `LU` in Phase 05 -- roughly half the
   factorisation memory -- but only above DC and only with gauge Method A,
   and they degenerate in opposite directions as omega falls.
   *Recommend: F3 as the default until a validated field exists to compare
   against, then choose between them by the Phase 05 frequency sweep, which
   is what `CONDITIONING.md` says to do and what `recommend_strategy`
   already takes measured arguments for.* Getting one correct answer out
   first has been the right call at every previous fork here.
