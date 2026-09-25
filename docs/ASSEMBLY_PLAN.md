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

### Symmetric and unsymmetric, from one code path

Both forms are assembled, selected by a flag, per `ROADMAP.md` Phase 04
step 6 ("keep all forms available since later phases compare them").

The symmetric form is `CONDITIONING.md` **Formulation 1** — divide the Φ
equation rows by `jω`. Because `α = jω·β`, that division does not *cancel*
a `jω`; it means one was never needed:

```
C_PhiA / (jω)  =  (jω · C_APhiᵀ) / (jω)  =  C_APhiᵀ      exactly
```

So one scalar covers both modes, applied to every Φ-row contribution:

| | `phi_row_scale` | (Φ,A) scatter | (Φ,Φ) scatter | Φ and port RHS |
|---|---|---|---|---|
| unsymmetric | `1` | `jω · C_APhiᵀ` | `K_PhiPhi` | unscaled |
| symmetric | `1/(jω)` | `C_APhiᵀ` | `K_PhiPhi/(jω)` | `÷ jω` |

The kernels themselves are unchanged — the flag only alters what multiplies
a Φ-row block on its way into the matrix — and **the sparsity pattern is
identical**, so it is still built once and shared between modes and across
a sweep.

Four consequences, each of which needs a test:

1. **Port rows scale too.** `V_k` is Φ-like: its row *is* the port-current
   equation. `PORT_AND_DOF_PLAN.md` §8 records that Formulation 1 must
   divide the `V_k` rows **and their `I_k` right-hand sides** by `jω`.
   Missing it makes the read-back current wrong by exactly `jω`, which
   reads as a physics error rather than a bookkeeping one.
2. **Symmetric requires ω > 0.** `1/(jω)` is undefined at DC. Not a real
   limitation — at ω = 0 the system decouples (`FORMULATION.md` §2) and
   there is no coupling block to symmetrise — but the code must *refuse*
   symmetric + DC, not divide by zero.
3. **No recovery step.** Row scaling multiplies an equation by a constant
   and leaves the unknowns alone, so Φ comes out unchanged. (Formulation 2,
   `Φ = jωΦ'`, *does* change the unknown and needs
   `recover_scaled_scalar_potential`. Keeping the two straight matters.)
4. **Symmetry needs gauge Method A.** `CONDITIONING.md`'s own table: the
   Munteanu projection reintroduces asymmetry regardless, so Method D +
   Formulation 1 pays the DC-degeneracy cost for no symmetry at all. Reject
   that combination rather than produce a non-symmetric matrix while
   claiming otherwise.

**Deferred: storing only one triangle.** A complex-symmetric matrix needs
half the storage, and MUMPS accepts that form. Doing it now would mean
guessing the solver's expected layout before Phase 05 has chosen one, and
it complicates the cross-check in §7.3. Assemble full storage in both modes;
revisit when the solver is actually linked.

Signature — **the caller owns the memory**, which is what keeps element
matrices off the heap and makes the kernels testable in isolation:

```cpp
void kernel_AA(const TetGeometry& g, const Mesh& mesh, int tet,
               const ElementCoefficients& c, std::complex<double>* out);   // 36
void kernel_APhi(const TetGeometry& g, const Mesh& mesh, int tet,
                 const ElementCoefficients& c, std::complex<double>* out); // 60
void kernel_PhiPhi(const TetGeometry& g,
                   const ElementCoefficients& c, std::complex<double>* out); // 100
```

`kernel_PhiPhi` needs no mesh or tet: P2 shape functions carry no
orientation, only the barycentric gradients, which are in `TetGeometry`.
`kernel_AA` and `kernel_APhi` do, because the Whitney functions need
`tet_edge_signs`.

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

On the cylinder (26 907 unknowns, 16 040 tets, ~1.5 M nonzeros expected):

| | |
|---|---|
| CSR `col_index` | ~6 MB |
| CSR `values` (complex) | ~24 MB |
| element matrices | **~2.6 KB, on the stack, reused** |
| *(triplet alternative)* | *~55 MB plus a 2.3 M-element sort* |

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
| scaling by `s`: curl-curl `→ s`, mass `→ s³`, coupling `→ s²`, `∇S·∇S → s` | dimensions | a volume factor dropped or applied twice |

### 7.3 Symmetric vs unsymmetric

- **The symmetric assembly really is symmetric**: `M == Mᵀ` entry for entry
  after assembly, not merely in the two coupling blocks.
- **Both modes describe the same problem.** Scaling rows does not change the
  solution, so for the same mesh and frequency the two assembled systems
  must have the same solution. Until a solver exists, check the weaker but
  still exact statement: multiplying every Φ and port row of the symmetric
  matrix by `jω` reproduces the unsymmetric one entry for entry, RHS
  included.
- **Symmetric + DC is refused**, not divided by zero.
- **Symmetric + gauge Method D is refused**, since the projection would
  break the symmetry the caller asked for.
- **The pattern is mode-independent**: `build_sparsity` output is identical
  for both, which is what lets a sweep share it.

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
`jω` on the unsymmetric PhiA scatter; forget to scale the port rows in
symmetric mode.

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
7. **The symmetric mode**, which is one scalar and its four consequences
   (§2), checked against the unsymmetric one by the row-scaling identity in
   §7.3.
8. **Time it**, then decide about threading with a number in hand.

---

## 9. Decisions to confirm

1. **Where the edge orientation sign is applied.** The DOF map already
   carries `tet_edge_signs` in its `coeff`. If the kernels *also* use the
   `_global` basis, the sign is applied **twice and cancels** — and no
   current test would catch that. *Recommend: kernels use the `_global`
   basis; `DofMap`'s edge `coeff` becomes 1, with a test asserting it is
   applied exactly once.* Either choice is correct; both is a silent bug.
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
5. **Symmetric and unsymmetric both, from one flag** (§2). *Agreed 25 Sept.*
   The default for a solve is the open question, not whether to have both:
   symmetric buys `LDLᵀ` over `LU` in Phase 05 — roughly half the
   factorisation memory — but only above DC and only with gauge Method A.
   *Recommend: unsymmetric as the default until a real solve exists to
   compare against, then switch once symmetric is shown to give the same
   answer.* Getting one validated field out first has been the right call
   at every previous fork in this project.
