# Field post-processing: exactly what `compute_fields` does

How a solution vector becomes `A`, `B`, `H`, `E` and `Φ`. Every formula here is
the one in `src/postprocess.cpp`; where the code makes a choice, this says which
choice and why. The design rationale is `POSTPROCESSING_PLAN.md` §2; this is the
implementation.

Companion documents: `FORMULATION.md` (the equations being solved),
`POSTPROCESSING_PLAN.md` §9 (what the cylinder validation measured).

---

## 1. What the solution vector holds

`DofMap` lays it out as three blocks:

```
x = [ a  (free edges) | phi  (free P2 nodes) | V  (one per port) ]
      0 .. num_a         num_a ..              num_a + num_phi ..
```

`Solution` carries the vector together with the `Conditioning` it was obtained
under and `omega`, because under `Conditioning::ScaledPhi` **the stored unknown
is not Φ**. The substitution is `Φ' = Φ / (jω)`, so recovering Φ needs

```
    phi_scale = jω    under ScaledPhi
              = 1     under Natural and RowScaled
```

`RowScaled` scales the Φ *rows* — it changes the equations, not the unknown — so
it needs no factor. This is applied in exactly one place, `Solution::phi_scale()`,
rather than at each call site; omitting it is listed as a negative control in
`SOLVER_PLAN.md` §8 and is the likeliest way to produce a wrong field that still
looks smooth.

**A prescribed value is never scaled.** It was written into the right-hand side
as a physical volt and never went through the substitution, so it must not come
back through it either.

---

## 2. The reference element

Barycentric coordinates `L_0 … L_3` with `Σ L_i = 1`. `TetGeometry` holds
`grad_L[i]` (constant over the tet) and the signed `volume`.

Local edges, `Mesh::kTetLocalEdgeVerts`:

| local edge | vertices | | local edge | vertices |
|---|---|---|---|---|
| 0 | (0, 1) | | 3 | (1, 2) |
| 1 | (0, 2) | | 4 | (2, 3) |
| 2 | (0, 3) | | 5 | (1, 3) |

The 10 local P2 nodes: `0…3` are the tet's vertices, and `4+e` is the midpoint
of local edge `e` — the same ordering, so local P2 node `4+e` sits on local
edge `e`.

---

## 3. Basis functions

### `A`: first-order Whitney (Nédélec) edge functions, 6 per tet

For local edge `e = (i, j)`:

```
    W_e(x)      = L_i ∇L_j  -  L_j ∇L_i
    curl W_e    = 2 ∇L_i × ∇L_j                     ← constant over the tet
```

The defining property is `∫_{e'} W_e · dl = δ_{e e'}`: the unknown `a_e` **is**
the line integral of `A` along edge `e`. That is what makes the uniform-`B` test
in §9 possible.

### `Φ`: second-order (P2) nodal functions, 10 per tet

```
    vertex i      N_i      = L_i (2 L_i - 1)
                  ∇N_i     = (4 L_i - 1) ∇L_i

    midpoint of
    local edge e  N_{4+e}  = 4 L_{v0} L_{v1}
    = (v0, v1)    ∇N_{4+e} = 4 ( L_{v1} ∇L_{v0} + L_{v0} ∇L_{v1} )
```

---

## 4. Orientation: the sign, applied exactly once

A tet's local vertex order comes from the mesh file and is arbitrary, so the
local pair `(v_i, v_j)` runs *against* the global edge's canonical low→high
direction for **roughly 58 %** of (tet, local edge) pairs on a real mesh. The
global-oriented basis is

```
    W^glob_e = sign(t, e) · W^loc_e          sign = Mesh::tet_edge_signs
```

Assembly and post-processing apply this sign in **different but equivalent**
places, and the difference is worth stating because getting it wrong is silent:

| | basis used | sign applied |
|---|---|---|
| assembly | local `W^loc` | at scatter, as `DofEntry::coeff` |
| `compute_fields` | global `W^glob` | inside the basis function |

So post-processing uses the **raw** global unknown and does *not* multiply by
`coeff` again. Applying the sign twice, or not at all, negates `A` on 58 % of
edges — and **no residual would reveal it**, because the residual is computed
against the same matrix that was assembled. A test (§9) pins the two conventions
together.

`Φ` has no such issue: its `coeff` is always `1`.

---

## 5. Resolving one tet's coefficients

`DofMap::local_dofs(t, mesh, bound)` returns the 6 edge and 10 Φ entries for tet
`t`. It is consulted **per tet**, not per node, for one reason: a node on an
internal cut reads the port's unknown or a grounded zero depending on **which
side this tet is on**, so the mapping is not a property of the node alone.

```
    edge L:   a_L  = x[ edge[L].index ]         if index ≥ 0
                   = 0                          otherwise

    P2 node b: φ_b = phi_scale · coeff · x[ phi[b].index ]   if index ≥ 0
                   = phi_fixed[b]                            otherwise
```

An eliminated **edge** is zero for one of two reasons, and the distinction
matters even though the value is the same: `n × A = 0` on a Dirichlet surface is
*physics*; a tree edge is the *gauge*, an arbitrary choice. Pick a different
tree and you get a different `A` with the same `B` (§8).

An eliminated **Φ** node carries `phi_fixed[b]`, already physical — see §1.

---

## 6. The fields

Under the `e^{+jωt}` convention of `FORMULATION.md` §1:

```
    A(x)  =  Σ_L  a_L  W^glob_L(x)                       linear in the tet

    B     =  curl A  =  Σ_L  a_L  curl W^glob_L          CONSTANT per tet

    H     =  B / μ ,     μ = μ_r · μ_0                   per the tet's body

    ∇Φ(x) =  Σ_b  φ_b  ∇N_b(x)                           linear in the tet

    E(x)  =  -jω A(x)  -  ∇Φ(x)                          linear in the tet
```

`μ_0` and `ε_0` come from `aphi_solver/constants.hpp`, shared with assembly.
They were file-local in `element_matrix.cpp` until this work: `H = B/μ` computed
with a different `μ_0` than the one behind `ν = 1/μ` is wrong in a way no
residual would show.

`B` is *genuinely* piecewise constant — not approximated as such — because
`curl W_e` is constant. Any smoothness in a plotted `B` is something
post-processing added.

### Where each is evaluated

| output | evaluated at | why |
|---|---|---|
| `b_tet` | — | constant; no point needed |
| `a_tet`, `e_tet` | centroid `L = (¼,¼,¼,¼)` | for a linear field this is its mean over the tet |
| vertex contribution | that corner, `L = e_c` | keeps the linear variation a centroid value would flatten |

Evaluating at the **corner** rather than reusing the centroid value is a
deliberate accuracy choice: `A` and `E` vary linearly within the tet, so a
centroid value smeared to all four vertices would discard exactly that
variation.

---

## 7. From per-tet to per-node

All four field arrays are sized `num_p2_nodes` and indexed **identically to
`phi_node`** — `DofMap`'s numbering, vertex `v` at index `v` and the midpoint of
edge `e` at `num_nodes + e` — so a field and the potential are read at the same
index with no second convention.

### At a vertex: volume-weighted over the incident tets

```
    X_v  =  ( Σ_{t ∋ v}  |vol_t| · X_t(v) )  /  ( Σ_{t ∋ v}  |vol_t| )
```

The absolute value is not cosmetic: a signed volume leaking through would let
contributions cancel. The test in §9 checks that the accumulated weight is
exactly `4 ×` the mesh volume, which is what catches that.

A vertex with no incident tet keeps weight `0`, is left at zero, and is counted
in `num_orphan_vertices` rather than silently dividing by zero.

### At an edge midpoint: the mean of the two endpoints

```
    X_m  =  ½ ( X_{v0} + X_{v1} )          for A, B, H, E  —  and ONLY these
```

`A`, `B`, `H` and `E` have **no mid-edge degree of freedom**. There is nothing
solved there to read, so the endpoint mean is the standard reconstruction.

---

## 8. Φ is copied, never averaged

**This is the one place the treatment of Φ and of the fields deliberately
differs, and it is the easiest thing in this document to "fix" by mistake.**

`Φ` lives in a **P2** space. The midpoint of every edge carries a genuine
unknown that the solve determined. Averaging the two endpoint values would
replace the quadratic with its linear interpolant and discard precisely the term
that makes the space second order.

In the code, `phi_node` is assigned exactly once:

```cpp
    out.phi_node = phi.value;       // from potential_at_nodes -- a copy
```

and the mid-edge loop writes only `a_node`, `b_node`, `h_node`, `e_node`. Φ is
never an operand of an average anywhere in `compute_fields`.

Solvers that *do* average Φ at mid-edge nodes do so because their Φ is P1 and
the mid-edge node exists only to draw a curved element. Ours is not that.

| | mid-edge degree of freedom? | mid-edge value |
|---|---|---|
| `Φ` | **yes**, a solved P2 unknown | **copied, exact** |
| `A`, `B`, `H`, `E` | no | mean of the two endpoints |

---

## 9. What is gauge-dependent, and what is not

The discrete system is **exactly gauge-invariant**. Under `A → A + ∇ψ`,
`Φ → Φ - jωψ` the coupling term becomes

```
    αA + β∇Φ + (α - jωβ) ∇ψ
```

and `α = jωβ` is an identity here, while `curl ν curl ∇ψ = 0`. The discrete
freedom is exactly `ψ ∈ P1`: gradients of P1 hat functions lie in the Whitney
space, and there are `#vertices - 1` of them modulo a constant — precisely the
number of tree edges that tree-cotree pins.

Consequences, which govern how these outputs may be read:

| quantity | gauge-dependent? |
|---|---|
| `Φ`, `A` | **yes** — the tree picks the representative |
| `B`, `H`, `E` | no |

Because the shift is `-jωψ` with `ψ` real to leading order, it is **almost
purely imaginary**: `Re Φ` is nearly clean while `Im Φ` carries an arbitrary P1
function. Measured on the 50 Hz cylinder (`POSTPROCESSING_PLAN.md` §9): `Re Φ`
matches the exact `z/l` to 4.02e-07, while `Im Φ` violates the problem's own
mirror antisymmetry by twice its own peak.

**So `phi_imag` and `a_*` are not results on their own. `E`, `B` and `H` are.**

---

## 10. What the tests pin down

`tests/test_postprocess.cpp`. The exact-reproduction tests run against a
hand-built **all-free** DOF map, because the real map pins a spanning tree and
the Dirichlet edges, so an arbitrary `A` is simply not in the constrained space
— that is a property of the constraint, not a correctness question, and the two
must not be conflated.

| test | what it would catch |
|---|---|
| global basis = sign × local, over every (tet, edge) | the §4 conventions drifting apart |
| …and the sign is negative somewhere (`flipped > 0`) | the above passing vacuously on a mesh with no flips |
| uniform `B` from `A = ½ B × r`, exact in every tet | any error in the edge reconstruction or orientation |
| linear `Φ`, `A = 0` ⟹ `E = -∇Φ` exactly; `B` exactly 0 | the P2 gradient, and the sign of `E` |
| vertex weight = exactly 4 × mesh volume | a signed volume cancelling |
| `b_node` at a vertex recomputed independently from `b_tet` | the weighting |
| `phi_node` equals `potential_at_nodes` bit for bit under `ScaledPhi` | the `jω` applied twice or not at all |
| mid-edge fields are endpoint means **and** mid-edge Φ is not | §8 being "fixed" |

The mid-edge test is driven by a genuinely **quadratic** Φ (`x² + 2yz`). With a
linear Φ the endpoint mean and the solved value coincide, so the test would pass
whichever rule were applied and would prove nothing.

**Negative control**, run and recorded: replacing `whitney_edge_curl_global`
with the local form fails exactly the two `B` checks and exits 1.

---

## 11. What is not here yet

Stated so nobody looks for it:

- **`J = σE`** is not produced. It is `E` times the tet's body `sigma`, and is
  meaningful only in conductors.
- **`sample_fields` / `probe_line`** (`POSTPROCESSING_PLAN.md` §3) are not
  implemented — both need point location, a separate concern.
- **No writer carries these fields yet.** `write_vtk` still writes Φ only. The
  plan §3 specified a Gmsh `$NodeData` writer on the grounds that Gmsh was
  already installed and no second tool would be needed; that reasoning has since
  expired, because ParaView is installed and `write_vtk` works, so extending the
  existing VTK writer is the better option. This is a deviation from the plan and
  is flagged rather than made silently.
- **`R` and `L` extraction** (plan §4) depends on the above.
