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

---

## 12. Material interfaces: where the nodal average is meaningless

Added after the first field run, because the defect it describes is invisible in
a plot and was found only by checking a number.

`E`'s **normal** component genuinely jumps across a conductor/insulator
interface -- `J_n = 0` at a free conductor surface requires it -- and `H`'s
tangential component jumps across a change of `mu`. A vertex on such a surface
has incident tets on **both** sides, so the volume-weighted average of §7 mixes
two physically different fields and produces a value that is neither.

Measured on the 50 Hz cylinder, `E` in the wire against the exact `-1000 V/m`:

| | nodes | mean `Re(E_z)` | worst deviation | worst transverse |
|---|---|---|---|---|
| **per-tet** (`e_tet`, no averaging) | 7535 tets | **-999.9993** | 8.6e-04 | 1.2e-06 |
| nodal, `material_interface = 0` | 5440 | -999.9992 | 8.7e-04 | 9.3e-07 |
| nodal, `material_interface = 1` | 6598 | mixed | **126** | **2705** |

A transverse field of 2705 V/m on an axisymmetric problem, nearly three times
the axial field, is not a small error -- and a surface plot of `E_magnitude`
shows it as a bright rim that looks like a skin effect. It is not one. The skin
depth here is 9.35 mm against a 0.2 mm radius; there is no skin effect at all.

`FieldOutput::on_material_interface` is 1 at every such node (a mid-edge node
inherits it from either endpoint), and `write_vtk` emits it as
`material_interface`. **Threshold it to 0 before reading a nodal `E` or `H`.**
The per-tet arrays have no interface to straddle and need no such care.

This is a limitation of one value per node, not a bug to be fixed in place: the
field really is two-valued there. Resolving it properly would mean a value per
(node, body), which is a larger change than the outputs currently justify.

---

## 13. What `write_vtk` writes

With a `FieldOutput` supplied. The legacy format has no complex type, so each
phasor is a pair of real `VECTORS`.

**Point data**, one per P2 node, indexed as `phi_node`:

| array | |
|---|---|
| `phi_real`, `phi_imag`, `phi_magnitude` | Φ, exact P2 values |
| `phi_present` | 0 where Φ does not live at all |
| `A_real`, `A_imag` | **gauge-dependent** (§9) |
| `B_real`, `B_imag`, `B_magnitude` | |
| `H_real`, `H_imag` | |
| `E_real`, `E_imag`, `E_magnitude` | |
| `material_interface` | 1 where the nodal average is meaningless (§12) |

`*_magnitude` is the phasor amplitude `sqrt(|Xx|² + |Xy|² + |Xz|²)`, which is
**not** the length of either the real or the imaginary vector: where the field
is elliptically polarised neither of those is the physical peak.

**Cell data**, one per tet, none of it averaged:

| array | |
|---|---|
| `body_tag` | the Physical Volume tag |
| `sigma` | the tet's body conductivity |
| `B_cell_real`, `B_cell_imag` | the **exact** piecewise-constant `B` |
| `J_real`, `J_imag` | `sigma * E_tet`, zero in an insulator |

Carrying both `B_cell` and the nodal `B` is deliberate: the difference between
them is the smoothing that §7 added, and nothing else makes it visible.

---

## 14. Material interfaces: each field takes its own side

A node shared between two materials has no single value for `E`'s normal
component or `H`'s tangential one. One side has to be chosen, and the rule is:

| field | side chosen |
|---|---|
| `E`, `J` | the **higher conductivity** |
| `B`, `H` | the **higher permeability** |
| `A` | neither -- it carries no material, so it averages over everything |

The two selections are tracked separately because they need not coincide. At a
copper/air surface the electrical side is the copper, while the magnetic one is
a tie (`mu_r = 1` on both) and no choice is made at all.

Each field is then averaged **over its own side only**, with its own weight:

    e_node = sum over tets with sigma = sigma_max   vol * E_t(node)  /  sum vol
    b_node = sum over tets with mu    = mu_max      vol * B_t        /  sum vol

and the two constitutive relations are applied at the node, exactly as written,
with the selected material:

    j_node = sigma_max * e_node          h_node = b_node / mu_max

So `J = sigma E` and `B = mu H` hold to the digit at every node, and neither
`sigma` nor `mu` is ever interpolated. Where a property is uniform across a
node's tets there is no ambiguity and nothing is excluded.

**Restricting the FIELD matters as much as choosing the material.** `sigma` from
the copper multiplied by an `E` averaged over copper *and* air still carries the
air's contribution. Measured on `01_OneCylinder` at 50 Hz, nodal `J` at the
conductor surface:

                                   |Jz| / axis      |Jr| / |Jz|
    sigma_max x E over all tets          1.060           0.1417
    sigma_max x E over copper only       1.062           0.0006
    exact (Bessel, a/delta = 1.07)       1.062           0

`Jr` must be zero at a free surface -- no current leaves the conductor -- and
restricting `E` to the conducting side takes it from 14 % of `Jz` down to
0.06 %, a factor of 236, while bringing `Jz` onto the exact value.

### Mid-edge nodes are accumulated, not averaged from endpoints

This is a change from the rule in Sec. 7, and the interface selection forces it.
An edge running from a conductor-surface vertex out into the air has one
endpoint carrying the copper-side field and one carrying the air's; their mean
is neither, and no choice of material at the midpoint can repair it.

Accumulating at the midpoint instead gets it right for free: a conducting tet's
four vertices all lie in the conductor, so such an edge belongs to **no**
conducting tet, its `sigma` range is uniformly the air's, and `J` there is zero
-- correctly, since that midpoint is inside the insulator.

The cost is that a mid-edge field no longer equals the mean of its endpoints.
Against a field whose answer is known -- a linear `Phi` with `A = 0`, so `E` is
constant -- the direct evaluation is exact, and a test pins that.

`Phi` is untouched by all of this: it is a solved P2 unknown at every node,
copied, never averaged, never assigned a side.

### Which to integrate

`j_tet` and `b_tet` -- exact, no averaging, no interface to straddle.
`tools/pv_extract_rl.py` reads `j_tet`. The nodal forms exist because they plot
smoothly.

### Measured: change the gauge, watch Phi move and E not

Asserted several times in this document, and eventually doubted, so it was run.
The spanning tree was changed two ways on the **identical mesh, identical node
numbering, identical everything else**, and the case re-solved:

| gauge change | worst change in Phi | worst change in E |
|---|---|---|
| BFS neighbours reversed | 0.0605 V | 3.5e-11 V/m (4e-12 of peak) |
| **breadth-first -> depth-first** | **0.7548 V** | 8.8e-11 V/m (1e-11 of peak) |

`E` is gauge-invariant to round-off. `Phi` moves by three quarters of the
applied volt.

And `Phi`'s SHAPE is not merely shifted. Along the conductor surface:

    z/L        0.15    0.35    0.55    0.75    0.95
    BFS       0.026   0.068   0.609   0.966   0.992      rises 0 -> 1
    DFS       0.547   0.566   0.649   0.691   0.797      never reaches either

Its azimuthal scatter at mid height changes too, 0.31 to 0.11. **What a plot of
Phi looks like is a choice, not a result.**

One intermediate result is worth keeping because it nearly misled. Reversing
the BFS neighbour order left the azimuthal scatter *identical* to three decimals
(0.3098 against 0.3095), which looked like evidence that the scatter was
physical rather than gauge. It was not: that perturbation changed only 250 of a
million stored nonzeros -- the same tree, essentially -- so it sampled almost
none of the gauge space. A weak perturbation showing no change is not evidence
of invariance.

### When `Phi = z/l` is the right target, and when it is not

The objection that finally settled this: tree-cotree exists to give a unique
`A` and `Phi`, and on a simple cylinder `Phi` should come out `z/l`.

Both halves are right, with one distinction each.

Tree-cotree makes the **matrix** non-singular and gives a unique solution *for a
given tree*. It does not give a tree-independent one: `A_tree = 0` is a
condition whose meaning depends on which tree, which is why the DFS variant
above moved `Phi` by 0.75 V.

And `Phi = z/l` is the **resistive-limit** answer. Measured on the same mesh,
same solver, same gauge, changing only the frequency:

    f       omega*L/R    worst |Re(Phi) - z/L|    profile at z/L = .15 .35 .55 .75 .95
    1 Hz      0.059            0.0021            0.150  0.352  0.553  0.754  0.961
    50 Hz     2.96             0.5420            0.021  0.048  0.633  0.973  0.994

At 1 Hz `Phi` is `z/l` to 0.2 %. At 50 Hz the same geometry is not a resistor:
of the 1 V applied, 0.33 V is resistive drop and 0.94 V is inductive EMF, in
quadrature. The potential distribution in an inductance-dominated structure is
not linear, and the split between `grad(Phi)` and `j*omega*A` is exactly what
the gauge fixes.

The boundary conditions are exact either way: both caps read 0 and 1 to the
last bit, with zero imaginary part.

**So `Phi = z/l` is a valid validation target only when `omega*L/R << 1`.**
Above that, a `Phi` that fails it is not a defect to chase.

There is a trap in this for the geometry design. Enlarging the conductor to make
the skin effect visible also makes the structure inductive -- both scale with
`omega*mu*sigma*a^2` -- so the case built to show one necessarily loses the
other. `examples/cylinder_50hz.aphi` (a = 0.2 mm) and
`regression_tests/01_OneCylinder/cylinder_1hz.aphi` are the resistive controls;
the 50 Hz cylinder is the inductive one.

### The tree-cotree gauge inflates |A| by an order of magnitude

Measured after a user pointed out that Ansys Maxwell, using the same tree-cotree
method and the same flux-tangential boundary, shows a continuous potential over
a cylinder surface.

    case                  mean |A| in conductor   mu0*I/(4*pi)   ratio
    0.2 mm wire, 50 Hz          4.84e-03 Wb/m       7.21e-04      6.71
    10 mm rod,  50 Hz           1.54e-01            1.45e-02     10.66

A Coulomb-gauge solution would sit near the physical scale. Ours is 7-11x above
it in **both** cases, so the inflation is a property of the gauge as
implemented, not something that appeared with the larger geometry. What differs
between the two is only how `omega*A*l` compares with the drive:

    0.2 mm wire   omega*A = 1.52 V/m over 1 mm   = 1.5e-03 V   0.15 % of 1 V
    10 mm rod     omega*A = 48 V/m over 40 mm    = 1.9 V        190 % of 1 V

which is why one case has a clean `Phi` and the other does not.

**Using tree-cotree does not fix which tree.** Another implementation of the
same method can produce a much smaller `|A|`, and therefore a much smaller gauge
term in `Phi`, without being a different formulation.

The Ansys A-Phi technical documentation is consistent with everything measured
here. It states the same gauge ("The degrees of freedoms of the edge elements on
the spanning tree of the finite element mesh are set to zero") and the same
invariance property ("the potentials at the terminal will be floating so they
are not unique. However, B, J, H, and all other quantities will be unique since
they depend on the derivative of the potential"). It also describes its voltage
plot as showing "the total potential which includes the ohmic electric potential
and the contribution from the eddy effect" -- an assembled quantity rather than
the bare nodal unknown, which may be why its plots look smooth.

**Open work.** A Coulomb-gauge projection as a post-process -- solve for `psi`
with `div(A + grad psi) = 0`, then `A' = A + grad psi`, `Phi' = Phi - j*omega*psi`
-- would make `Phi` unique, shrink `|A|` toward the physical scale, and leave
every gauge-invariant quantity untouched. It is testable: the BFS and DFS trees
must then give the same `Phi`.

### What sets Phi's smoothness under tree-cotree, and what does not

The gauge function is `psi(v) = -integral A.dl` along the tree path from the
root, so for two ADJACENT nodes joined by a cotree edge

    psi(v) - psi(u) = -closed integral A.dl round the fundamental cycle
                    = -flux enclosed by that cycle

**`psi` jumps between neighbours by exactly the flux their fundamental cycle
encloses.** Short cycles enclose little and `psi` is smooth; a cycle that runs
the length of the conductor encloses the whole flux linkage and `psi` jumps
there by the EMF.

That predicted the seam: the outer boundary is one group, BFS advances inward
from both caps, the fronts meet at `z = L/2`, and a node there has tree paths
to OPPOSITE caps. Observed jump across mid height ~0.6 V against an EMF of
`omega L I = 0.94 V` -- right size, right place.

**The prediction was wrong about the cure.** Re-entering the boundary nodes in
order of increasing `z`, so the sweep runs from one cap, changed nothing:

    worst |Re(Phi) - z/L|     0.5420  ->  0.5420
    azimuthal std at mid       0.3243 ->  0.3219

Two reasons, both structural:

1. **BFS from a SET gives every node its nearest boundary point** whatever order
   the set is queued in. Reordering changes which neighbour is discovered first
   within a level, not the depth structure, so the fronts still meet where they
   met before.
2. **`psi` constant on the whole boundary is forced, not chosen.** `n x A = 0`
   on every boundary edge means walking along the boundary accumulates no
   `psi`, so the entire box surface is one equipotential for it -- for any tree.

Since the conductor spans from one boundary face to the other, loops with both
ends on the boundary enclose the full flux linkage. `psi` must therefore vary by
~the EMF somewhere inside, and **no choice of spanning tree removes that.** What
a tree can do is move where the variation sits and whether it is concentrated or
spread; it cannot make it go away.

So on this geometry, with `flux_tangential` on the whole box and a large flux
linkage, tree-cotree alone cannot deliver a smooth `Phi`. On the resistive cases
it can, because the EMF that has to be absorbed is negligible: 1.5e-03 V on the
0.2 mm wire against a 1 V drive.

### The regime is set by `omega*L/R`, not by frequency

Written down after comparing against an Ansys Maxwell voltage plot that shows a
perfectly continuous potential over three conductors -- at 50 Hz, the same
frequency as the case here that does not.

There is no contradiction, because the parameter is not frequency. For a round
conductor `omega*L/R` scales as `omega * a^2`, so radius and frequency enter the
same way. From the Maxwell case's own reported terminal impedance,
`R = 3.4575e-04 ohm` and `omega*L = 1.2395e-05 ohm` at 1 A:

    pi a^2 = l / (sigma R) = l / 20053     ->   a ~ 1.26 mm for l = 100 mm

                    a        a/delta @50Hz    omega*L/R
    Maxwell rods   ~1.26 mm      0.13           0.036
    this case      10 mm         1.07           2.96

Scaling by `a^2`: `(10/1.26)^2 = 63`, and `0.036 * 63 = 2.3` against the measured
2.96. Both cases sit where their geometry puts them.

**`a/delta` and `omega*L/R` are the same parameter.** Making the skin effect
visible (`a/delta >~ 1`) necessarily makes the structure inductive
(`omega*L/R >~ 1`), which is exactly when `Phi` stops being readable. A case
cannot show both a strong skin effect and a clean potential. The Maxwell plot is
continuous *because* those rods have no skin effect to show.

Verified on this solver, same mesh and same tree-cotree gauge, at
`omega*L/R = 0.059`: `Phi` is a smooth gradient, `z/l` to 0.2 %. The
0.2 mm wire at 50 Hz (`omega*L/R = 8.8e-04`) gives `z/l` to 4e-07.
