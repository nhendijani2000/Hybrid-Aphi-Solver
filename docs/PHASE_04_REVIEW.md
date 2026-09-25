# Phase 03 → 04 review and revised plan

*Opus 5 review, Sept 2026. Findings are from reading the actual code and
running the actual build, not from the design docs' own account of
themselves. The plan in Section 4 is mirrored into `docs/ROADMAP.md` (Phase
02 step 5, the new Phase 03.5, and the rewritten Phase 04) — the roadmap is
the source of truth; this document is the reasoning behind that edit.*

**Supersedes** `docs/DOF_AND_INPUT_FILE_PLAN.md`. That draft's
mesh-tag gap (its §2) and gauge-reduction index-space gap (its §5) were both
real and are carried forward; its open questions 1 and 2 are answered below;
its implementation sequence is replaced by Section 4.

---

## 1. Diagnosis

Phase 04 was not stalling for a mathematical reason. As written it silently
required four interlocking prerequisites that did not exist yet:

1. a complex **sparse** matrix type,
2. mesh region/boundary tag ingestion,
3. a gauge reduction that can act on part of a larger block system,
4. correct global edge orientation for the **A** basis.

Each of these changes the type signatures the others are written against.
That is why every attempt at assembly code surfaced "one more gap" and the
work produced plans instead of matrices. The fix is structural: split them
into their own phase with their own exit gate, so Phase 04 is about the weak
form and nothing else.

Of the four, item 4 was an outright bug rather than a missing piece.

---

## 2. Findings

### 2.1 Critical — missing global edge-orientation sign (fixed, commit `919cce5`)

`whitney_edge_value(g, local_edge, L)` returned the Whitney function for the
**local** vertex order `tv[vi] -> tv[vj]`. A global edge DOF is defined on
the canonical low-index → high-index direction (`mesh.cpp`, `edge_key`).
Where the two disagree, the returned function is the **negative** of the
global edge's basis function. No sign existed anywhere in the codebase.

Measured directly from the project's own meshes:

| mesh | nodes | tets | edges | local edges needing a sign flip |
|---|---|---|---|---|
| `cube_2.msh` | 27 | 48 | 98 | 168 / 288 (58%) |
| `cube_4.msh` | 125 | 384 | 604 | 1344 / 2304 (58%) |
| `cube_6.msh` | 343 | 1296 | 1854 | 4536 / 7776 (58%) |

Concretely: `cube_4.msh` tet #2 is `1 7 6 32` → local verts `[0,6,5,31]`;
local edge 3 spans global nodes (6,5), so assembly would have added `−N`
under the DOF for canonical edge (5,6).

**Why it survived to Phase 04.** Both single-tet test fixtures used vertex
order `{0,1,2,3}` — ascending, so every sign is +1 and an orientation error
is structurally invisible. The convention was *documented* (it is the
`mEdgeSign` of the prior 3dedyaphi implementation, cited in
`docs/FORMULATION.md` §5.1) but never implemented.

**Fix.** `Mesh::tet_edge_signs` (built in `build_topology` alongside
`tet_edges`) plus `whitney_edge_value_global` / `whitney_edge_curl_global`,
which are the forms assembly must call. Φ needs no sign — its P2
edge-midpoint shape function `4*L_v0*L_v1` is symmetric in its two vertices.

**Regression cover, verified to actually catch it.** A two-tet mesh whose
second tet lists vertices non-ascending, checking (a) the global circulation
identity `δ_jm` along each edge's canonical direction and (b) tangential
continuity across the shared face — the H(curl) conformity the sign exists
to protect. With the sign neutralized the new checks fail on exactly the
reversed edges (426/431); restored, 431/431. Full suite: 564 checks green.

### 2.2 Architectural — there is no complex sparse matrix type

The codebase has two incompatible halves:

- `SparseMatrix` (`incidence.hpp`) — real-valued COO, `std::map`-based
  coalesce, naive multiply. Its own header calls it a Phase-02-scale
  placeholder.
- `ComplexMatrix` (`complex_matrix.hpp`) — complex, but **dense**.

`APhiBlockSystem` is built from the dense one. Phase 04 needs complex
*sparse*, which neither provides. A dense `K_AA` is already ~55 MB at
`cube_6.msh` and impossible at any real EDA mesh size.

`docs/ENGINEERING_STANDARDS.md` item 3 already requires compressed sparse
formats "from Phase 04 assembly onward" and explicitly names this
`SparseMatrix` as the thing to replace — but nothing had been assigned to
do it. The earlier draft plan framed this as "template the gauge functions
or write a complex copy," which treats the symptom; the blocker is that
`APhiBlockSystem` itself is the wrong type.

### 2.3 Mesh ingestion discards every tag

Carried forward from the earlier draft, which had this right.
`read_gmsh_msh` reads and throws away every element tag — including each
tet's own physical-group id — and skips non-tet elements entirely, so the
boundary triangles Gmsh uses to mark PEC walls, ports and the outer
truncation never reach `Mesh`. DOF setup, material lookup and PEC handling
all key off region membership, so this is a hard prerequisite.

### 2.4 Gauge reduction assumes a bare edge-indexed matrix

Also carried forward. Both gauge builders assume the matrix index space *is*
the mesh's global edge index space, 1:1 — true for their current callers
(which only ever pass `M = CᵀC`), wrong for a coupled block system where Φ
rows and columns must pass through untouched.

### 2.5 Minor — two API collisions

`dense_solve` (real, `gauge_variants.hpp`) vs. `solve_dense` (complex,
`complex_matrix.hpp`) differ only in word order, and
`estimate_condition_number` exists twice with different types in different
headers. Survivable now; actively confusing once the scalar type becomes a
template parameter.

---

## 3. Decisions taken

| Question | Decision | Rationale |
|---|---|---|
| Basis pairing | **Keep mixed order** (1st-order Whitney **A**, P2 Φ) as locked in `docs/FORMULATION.md` §5.1 | Owner's call. Reviewed and re-affirmed; see the cost note below. |
| Φ support in full-wave | **Everywhere**, whole domain | With displacement current on, Φ carries physical content outside conductors. `Ω_c` restriction belongs only to the reduced low-frequency variant (`FORMULATION.md` §3). |
| First gauge | **Albanese–Rubinacci** (`a_t = 0` on tree edges) | Simplest reduction available: pure row/column restriction, a genuine submatrix, no fill-in, symmetry-neutral. Get a correct solve before making the gauge a user-facing knob. Munteanu unsymmetric stays implemented and tested but does not gate the first validated field. |
| Real vs. complex gauge reduction | **One template**, `Sparse<T>` | Answers the earlier draft's open questions 1 and 2 in a single move rather than maintaining a second complex-typed copy. |
| Input file sequencing | **`Problem` struct first, JSON second** | The struct is the interface DOF setup needs; JSON is a thin adapter. Decoupling means DOF setup is not blocked on vendoring a library. |

**Cost note on the retained mixed order, recorded for Phase 10.** Because
`N_Φ = N_vertices + N_edges` and `N_edges ≈ 5–7 × N_nodes`, P2 Φ costs ~1.85×
the total unknowns of a matched-P1 pairing (measured: 1.83× on `cube_4`,
1.84× on `cube_6`), which is roughly 3× the direct-factorization time, while
`B = curl A` stays first-order accurate either way. It also breaks
`grad(V_Φ) ⊆ V_A`, so equation (II) is no longer the exact discrete
divergence of equation (I). Neither is a correctness problem and the
decision stands; both are worth re-measuring against a working baseline at
Phase 08/10 rather than re-arguing from theory.

---

## 4. Plan

Mirrored into `docs/ROADMAP.md`. Summary only here.

### Phase 03.5 — numerical foundations (new, own exit gate)

1. ~~Global edge-orientation signs~~ — **done**, commit `919cce5`.
2. **One templated sparse type.** `Sparse<T>`, with `SparseMatrixD` for the
   real incidence operators **C**/**G** and `SparseMatrixZ` for the complex
   system blocks. Accumulate as triplets, `compress()` once to CSR, and
   **keep a triplet export** — MUMPS's assembled centralized format wants
   IRN/JCN/A arrays, not CSR, so both shapes are needed regardless. Migrate
   `APhiBlockSystem` onto it; keep dense `ComplexMatrix` only as the
   small-system ground truth for `solve_dense`.
3. **Mesh region and boundary tag ingestion.** `Mesh::tet_tags` plus a
   tagged boundary-face list from Gmsh elm-type 2 triangles. `Mesh` carries
   integers only — what a tag *means* is the input file's job.
4. **Generalize gauge reduction to an explicit A-DOF index set**, rather
   than writing block-aware duplicates of the same reduction math.
5. **Resolve the two API collisions** before the template multiplies them.

*Exit gate:* complex sparse matrices build/compress/matvec/export;
`APhiBlockSystem` holds sparse blocks; a tagged `.msh` round-trips its tags;
Albanese–Rubinacci reduces a system whose A-DOFs are a strict subset of the
index space, with the existing 19/19 gauge checks unchanged.

### Phase 04 — weak form and assembly

1. `Problem` struct (materials by tag, BCs, gauge, scaling, frequency,
   source); JSON adapter afterwards, cross-checking every tag against the
   loaded mesh.
2. DOF numbering — `N_A = N_edges`, `N_Φ = N_vertices + N_edges`, Φ
   everywhere in the full regime; `is_pec` derived from the new boundary
   tags. Testable before any assembly exists.
3. Element matrices, via the `*_global` basis forms. **A degree-2 (4-point)
   tetrahedral quadrature rule is exact for every integrand in this phase**
   — Whitney is linear so the mass and coupling terms are quadratic,
   `curl(Whitney)` is constant so curl-curl is exact at one point,
   `grad(P2)` is linear so `K_ΦΦ` is quadratic. Region-wise-constant
   materials add no degree.
4. **DC milestone first.** At `ω = 0` the system decouples exactly
   (`FORMULATION.md` §2) — real-valued, uncoupled, analytically checkable
   (straight wire, `B = μI/2πr`), while still exercising tags → DOF map →
   oriented element matrices → global assembly → gauge reduction → sparse
   solve. Turn `ω` on only once it is green.
5. Global sparse assembly, Albanese–Rubinacci first.
6. Frequency-scaling choice (`CONDITIONING.md` Formulations 1–3).
7. RHS assembly — uniform current density is enough here.

*Exit gate:* the DC milestone reproduces the straight-wire field; the
assembled sub-blocks are complex-symmetric where `CONDITIONING.md` says they
should be; and the same problem gives the right sign and magnitude of **B**
at `ω ≠ 0`.

---

## 5. Sequencing rationale

Two reorderings do most of the work. **Phase 03.5 step 2** removes the type
churn that made every assembly attempt reopen the same questions. **Phase 04
step 4** produces a checkable physical number roughly half the pipeline
sooner than a full complex coupled solve would. Writing assembly against a
dense block system *and* a full complex coupled solve simultaneously is what
made Phase 04 feel unbounded.
