# Ports and DOF assignment — plan

*Opus 5, Sept 2026. Draft for review — nothing here is implemented yet. Once
approved it folds into `docs/ROADMAP.md` Phase 04 (steps 1–2) and a new
`docs/FORMULATION.md` section on port and boundary conditions.*

---

## 0. Summary

1. **One mechanism covers all four port types.** Each port owns one extra
   global unknown, its voltage `V_k`. The equation row belonging to `V_k` is,
   identically, the port current `I_k` (derived in §2). So:
   - **current port** → `V_k` is unknown, `I_k` goes on the right-hand side;
     the solve returns `V_k`.
   - **voltage port** → `V_k` is fixed; the solve returns `I_k` from its row.

   Boundary vs. internal changes only *how Φ couples to `V_k`*: a boundary
   terminal **slaves** its Φ nodes to `V_k`; an internal cut makes Φ **jump**
   by `V_k` across the surface.
2. **One mechanism covers all DOF bookkeeping.** Every element-local DOF maps
   to a short list of `(global index, coefficient)` plus a fixed value. Edge
   orientation signs, tree edges, PEC walls, ground, terminal slaving, cut
   jumps and fixed voltages are all just different entries in that list.
   Assembly never special-cases any of them (§5).
3. **A required change, measured (§6).** The Phase 03 tree-cotree gives each
   PEC body its *own* root. With two or more separate `n×A = 0` surfaces —
   which one port terminal plus one ground wall already creates — the
   gauged system is **singular**. Fix: collapse each Dirichlet component to a
   single super-node and build one spanning tree.
4. **Two corrections to the DC milestone (§7).** At ω = 0 Φ must live only on
   conductors (its equation is empty in dielectrics). And the DC inductance
   to compare against is `(μ₀/2π)(ln(b/a) + ¼)·ℓ`; the `L'` in FORMULATION
   §7.1 is the external (high-frequency) inductance only, and would show a
   false **15.5 %** error at DC.

---

## 1. Terminals, ports, and sign conventions

- A **terminal** is a tagged surface in the mesh.
- A **port** is a terminal + a reference + an excitation.

| | terminal faces | how Φ couples | condition on **A** | you specify | solve returns |
|---|---|---|---|---|---|
| **boundary current** | boundary (1 adjacent tet) | Φ = `V_k` on the terminal | `n×A = 0` | `I_k` | `V_k` |
| **boundary voltage** | boundary | Φ = `V_k` on the terminal | `n×A = 0` | `V_k` | `I_k` |
| **internal current** | interior (2 adjacent tets), forming a cut | Φ jumps by `V_k` across the cut | none | `I_k` | `V_k` |
| **internal voltage** | interior, forming a cut | Φ jumps by `V_k` | none | `V_k` | `I_k` |

Boundary vs. interior is derivable — `Mesh::face_tets` already says how many
tets touch each face. The input file states it anyway (`kind`, §9), as a
declaration of intent that is *checked* against the mesh: a surface you
believed was on the boundary but is embedded in the volume is exactly the
mistake worth catching loudly.

**Conventions (one set for all four types):**
- `V_k = Φ(terminal) − Φ(reference)` for boundary ports; `V_k = Φ⁺ − Φ⁻` for
  a cut.
- `I_k` = current **entering** the structure through the terminal, or
  crossing the cut **from the − side to the + side**.
- Then `P_k = ½ Re(V_k I_k*)` is the power the port delivers into the
  structure, and `Z_k = V_k / I_k` is the impedance it sees. Phasors are
  **peak** amplitudes (that is what the ½ assumes).

---

## 2. Why the port row *is* the port current

Write the total (conduction + displacement) current density as
`J_tot = (σ + jωε)E = −(jωσ − ω²ε)A − (σ + jωε)∇Φ`. Equation (II) of
FORMULATION §1.4 is `∇·J_tot = 0`. Its weak form with a scalar test function
`Φ'`:

```
∫_Ω F·∇Φ' dΩ = −∮_∂Ω (n·J_tot) Φ' dS,     where F = −J_tot
```

The left side is what assembly builds; the right side is boundary data.

**Boundary terminal.** Φ on terminal `T_k` equals `V_k`, so testing with
`V_k`'s own basis means testing with `v_k = Σ_{i∈T_k} N_i` — the sum of the
P2 basis functions of every terminal node. `v_k = 1` on `T_k`. Elsewhere on
the boundary it is nonzero only on faces touching the terminal's rim — and
those carry no current data: validation (§10) forbids a terminal touching
ground, and a port sheet is natural for Φ (`n·J_tot = 0`). So the assembled
row is

```
(K x)_{V_k} = −∫_{T_k} n·J_tot dS = I_k          (n outward ⇒ −flux = inflow)
```

**Internal cut.** Φ = Φ_c + `V_k`·v_k with v_k = Σ_{cut nodes} N_i restricted
to the tets on the + side. Integrating by parts over that one-element layer,
v_k = 1 on the cut and 0 on the rest of the layer's boundary:

```
(K x)_{V_k} = ∫_cut n_cut·J_tot dS = I_k          (current crossing − → +)
```

Same identity, same sign. Consequences:

- **Current port:** put `I_k` in the right-hand side of the `V_k` row.
- **Voltage port:** `V_k` is known and moves to the RHS; after the solve,
  `I_k` is the residual of that row. No separate surface integral is needed —
  and this *weak* (reaction) current is more accurate than integrating
  `−σ∇Φ·n` over the terminal, because it uses the equation rather than a
  differentiated field.
- **Open circuit = current port with `I = 0`. Short circuit = voltage port with
  `V = 0`. A floating conductor = a current port with `I = 0`.** No extra
  machinery for any of them.
- **Z-matrix:** column *j* = excite port *j* with `I_j = 1`, every other port
  `I = 0`, read all `V_i`. One factorization, N right-hand sides — the
  "factorize once, cheap back-substitution per excitation" pattern Phase
  06/07 already plan for active impedance.

---

## 3. Attaching each type to the mesh

### 3.1 Boundary terminal (electrode)

- **Φ:** every P2 node on the terminal — each vertex *and* each edge midpoint
  of every terminal face — is slaved to `V_k`. The terminal is an
  equipotential.
- **A:** `n×A = 0` on every edge of every terminal face.

Why the *essential* condition on A and not the natural one, using the coax
itself: current enters the via's end face along z, so **B** circulates in the
plane of that face (`B_φ`). The natural condition `n×H = 0` would force
`H_φ = 0` on the face — wrong. `n×A = 0` implies `B·n = B_z = 0` and leaves
`H_φ` free — right. The assumption behind it, stated: current crosses the
terminal roughly normally, so **B** is tangential to it. That holds for any
terminal that is a conductor cross-section, which is what terminals are.

### 3.2 Port sheet (recommended for boundary ports)

A boundary port normally also has a surface *bridging* its terminal to its
reference — for the coax, the dielectric annulus at z = 0 (tag 21). Give it
`n×A = 0` with Φ **natural**:

- On the sheet `E_t = −∇_tΦ` (A has no tangential part), so
  `V_k = Φ_T − Φ_ref` equals `∫E·dl` along **any** path within the sheet.
  The port voltage is then gauge-invariant and path-independent, **even at
  full-wave**.
- Without a sheet (terminal and a distant ground only), `V_k` is still a
  Φ difference, but Φ differences are not gauge-invariant: with the Method A
  tree gauge, it equals `∫E·dl` along whichever tree path happens to join
  them. Exact at DC, increasingly arbitrary as the port grows electrically.

So a sheet is what turns a boundary port into a proper lumped port. Recommend
**required for ω > 0, optional at DC** (decision 2 in §13).

### 3.3 Internal cut

- **Mesh:** a surface embedded in the volume (Gmsh: include it in
  `BooleanFragments`) and tagged as a physical surface. Its faces have two
  adjacent tets.
- **Φ:** `Φ = Φ_c + V_k·v_k`, where `Φ_c` is continuous and
  `v_k = Σ_{i ∈ cut nodes} N_i` restricted to tets on the + side. Φ therefore
  jumps by exactly `V_k` across the cut. There is no DOF duplication; in the
  §5 map, a cut node seen from a + side tet simply has two entries.
- **A:** no condition. **A** is continuous across the cut.
- **Which side is +:** the face normal. That needs the triangle's original
  winding, which **the reader currently discards** (`TaggedFace` sorts its
  nodes) — step 1 of §12 keeps it. The input file may also give a
  `direction` vector, checked against every face.
- **Which tets are on the + side:** a tet that touches the cut only at a node
  or an edge (not across a face) still has to be assigned a side. Robust
  rule: flood-fill the star of tets around each cut node, seeded from the tet
  across a cut face on the + side, never crossing a cut face. Works for
  curved cuts, needs no geometric tolerance.
- **Gauge invariance:** under a gauge change Φ → Φ − jωψ with ψ continuous,
  the jump is unchanged. And since **A** is continuous and finite, `∫E·dl`
  across an infinitesimally thin cut equals the jump. So the cut voltage is
  well-defined at full-wave, like a port with a sheet.
- **Why cuts are necessary, not just convenient:** they are the only way to
  drive a **closed conductor loop** (a ring, a coil turn, a via-plane-via
  loop). A loop has no boundary face to attach an electrode to.

**Rim policy.** The cut's rim (edges belonging to exactly one cut face) must
lie on the boundary of Φ's support:
- At DC, Φ lives only on conductors (§7), so a cut across a conductor cross-
  section always qualifies.
- At ω > 0 in the full-wave regime Φ also lives in the dielectric, so a cut
  through a conductor only has its rim *inside* Φ's domain. That is a
  **delta-gap source**: restrict `v_k` to conductor tets, accept that Φ is
  discontinuous along the rim, and the field there is singular. That is the
  standard, well-known idealization — but it should be a deliberate choice
  (decision 3 in §13).

### 3.4 Internal electrode

A terminal patch on an *internal* surface — e.g. where a lead attaches to a
conductor inside an air box. Slaving (§3.1) does not care how many tets touch
a face, so this is the boundary-terminal code path with no changes. Supported
for free.

### 3.5 Reference

- **`"ground"`:** the grounded PEC surfaces (Φ = 0).
- **Another terminal (differential port):** terminal nodes map to
  `[(V_ref,1), (V_k,1)]` and reference nodes to `[(V_ref,1)]`, with one extra
  common-mode unknown `V_ref`. The `V_ref` row is `I_T + I_R = 0`: whatever
  enters one terminal leaves the other. Cheap in the §5 map, but it adds a
  validation case (decision 4).

---

## 4. Boundary conditions, per field

A and Φ need **independent** conditions — the coax port sheet needs
`n×A = 0` with Φ free, which no single "PEC / natural" switch can express.

| type | **A** | **Φ** | typical use |
|---|---|---|---|
| `pec` (grounded) | `n×A = 0` | Φ = 0 | shields, ground planes, shorts |
| `pec_floating` | `n×A = 0` | Φ = unknown constant, net current 0 | isolated conductor surface (internally: a current port with I = 0) |
| `flux_tangent` | `n×A = 0` | natural (`n·J_tot = 0`) | port sheets; planes that **B** is tangent to |
| `natural` (default) | natural (`n×H = 0`) | natural | insulating walls |
| `abc` | Phase 07 | Phase 07 | open boundary |

Terminal, sheet and cut roles are assigned through `ports`, not here.

**Coax assignment** (mesh tags from `tools/coax_via.geo`):

| tag | surface | assignment |
|---|---|---|
| 10 `shield_pec` | r = b wall | `pec` |
| 11 `short_end` | z = L caps | `pec` |
| 20 `port_conductor` | via end, z = 0 | **terminal** of P1 |
| 21 `port_return` | annulus, z = 0 | **sheet** of P1 (`flux_tangent`) |

With tag 21 left `natural` instead, `H_φ` would be forced to zero on the
feed plane and L would come out wrong. Tag 21 is misnamed for this role
(it carries no return current); rename to `port_sheet` (decision 6).

**Precedence where surfaces meet**, per node/edge:

| meeting | result |
|---|---|
| terminal ∩ ground | **error** — the port is short-circuited |
| terminal ∩ another terminal | **error** |
| terminal ∩ sheet | terminal wins (the sheet's inner rim *is* the terminal) |
| sheet ∩ ground | ground wins (the sheet's outer rim *is* the reference) |
| cut rim ∩ ground | **error** — Φ⁺ = Φ⁻ + V and Φ = 0 cannot both hold |
| `n×A = 0` from any source | union — an edge is Dirichlet if any rule says so |

---

## 5. DOF assignment

### 5.1 Unknown vector

```
x = [ a  (free edges)  |  φ  (free Φ nodes)  |  V  (port and floating unknowns) ]
```

Port unknowns go last: each couples to every DOF of its terminal or cut
layer, so their rows and columns are comparatively dense. MUMPS reorders
regardless, but last keeps them out of the way of the A/Φ structure.

### 5.2 Classification

**Edges:** `FREE` | `TREE` (a = 0, the gauge) | `DIRICHLET` (a = 0, `n×A = 0`).

**Φ nodes** (vertices and edge midpoints, P2):
`FREE` | `FIXED` (ground, Φ = 0) | `SLAVE(k)` (terminal of port k) |
`CUT(k)` (free, with a + side jump) | `ABSENT` (outside Φ's support — §7).

A voltage-port terminal is `SLAVE(k)` with `V_k` known, not `FIXED(V)`. The
distinction is what keeps its current readable as the `V_k` row.

### 5.3 The local → global map

For tet *t* and local DOF ℓ (6 edges + 10 Φ nodes), store a short list of
`(global index, coefficient)` and a fixed value:

| local DOF | entries | fixed value |
|---|---|---|
| free edge | `(a_e, s_e)`, `s_e` = `Mesh::tet_edge_signs` | 0 |
| tree or Dirichlet edge | — | 0 |
| free Φ node | `(φ_i, 1)` | 0 |
| ground Φ node | — | 0 |
| terminal node of port k | `(V_k, 1)` | 0 |
| cut node, tet on + side | `(φ_i, 1)`, `(V_k, 1)` | 0 |
| cut node, tet on − side | `(φ_i, 1)` | 0 |
| absent Φ node | — (the local DOF is dropped) | — |

Assembly scatters `c_a · c_b · K_e[ℓ_a, ℓ_b]` for every pair of entries.
Fixed values contribute to the RHS. Known `V_k` (voltage ports) are moved to
the RHS the same way.

The point of doing it this way: **the cut function, terminal slaving, and the
edge-sign fix from `919cce5` are not separate code — they are entries in this
table.** The map is also exactly the thing to unit-test before any element
matrix exists.

### 5.4 Gauge by elimination at assembly

Method A is "tree edges map to nothing": the gauge is applied while
assembling, rather than assembling everything and then taking
`principal_submatrix`. `principal_submatrix` stays as the independent check —
*assemble ungauged, then reduce* must equal *assemble directly reduced*,
entry for entry. That is a free, exact test.

Method D cannot be written this way: its test space (plain cotree selector)
differs from its trial space, so it stays a post-assembly reduction. Not on
the critical path, since Method A goes first.

### 5.5 Predicted counts on the coax

With the §4 assignment, the **entire** boundary is `n×A = 0` (shield, short,
terminal, sheet). The boundary is a closed genus-0 surface with 1,710
triangles, so it has 2,565 edges and 857 nodes (χ = 2). It is one connected
Dirichlet component, so the tree has `5,099 − 857 = 4,242` edges, and

```
free A DOFs = 34,539 − 2,565 − 4,242 = 27,732
```

That exact number becomes a DOF-map test.

---

## 6. Required change: the tree-cotree gauge with Dirichlet surfaces

Measured on `cube_4` with `M = CᵀC`, two opposite faces playing ground wall
and port terminal. Nullity of the reduced matrix by dense elimination with
complete pivoting; dropped pivots ≈ 1e−14, kept ≥ 0.05, so the counts are
unambiguous.

| case | tree edges | Dirichlet edges | nullity |
|---|---|---|---|
| (i) no PEC | 124 | 0 | **0** |
| (ii) 1 PEC face grouped, `n×A = 0` not imposed | 100 | 0 | **24** |
| (iii) 1 PEC face, `n×A = 0` imposed | 100 | 56 | **0** |
| (iv) 2 Dirichlet faces, each its own root *(Phase 03 design)* | 75 | 112 | **1** |
| (v) 2 Dirichlet faces collapsed, one spanning tree *(proposed)* | 76 | 112 | **0** |

**Why (iv) fails.** With *k* separate roots, take ψ = c_j, a different
constant on each tree. Then `a = Gψ` vanishes on every tree edge and every
Dirichlet edge but not on the cotree edges between trees — a gradient, hence
in curl-curl's null space, that survives the gauge. *k* − 1 null directions.
The tree in (iv) is exactly one edge short (75 vs 76).

**Why it happened.** Phase 03 step 4 merged two different questions:
- *Which conductors are at which potential?* — a question about **Φ**. Each
  body can be at its own potential; that is ground / floating / port data.
- *Where is the A-tree rooted?* — a question about the **gauge**. It must be
  one tree per connected mesh component, with each Dirichlet component
  collapsed to a single super-node.

They separate cleanly: the tree handles the second, the §5 map handles the
first.

**Case (ii)** shows the other half: grouping PEC nodes without eliminating
their edges leaves (#body nodes − 1) null directions. The grouping only makes
sense together with `n×A = 0` elimination, which Phase 04 adds — so both land
together.

**A second hazard, in how bodies are identified today.** `build_tree_cotree`
unions two PEC nodes if *an edge* joins them. In a thin substrate, a single
tet edge can span from a trace to the ground plane beneath it — and today's
rule would merge trace and ground into one body: a **short circuit**, in
exactly the EDA geometry this solver targets. The correct rule unions only
along edges that lie **on** a tagged Dirichlet face.

**Change:** `build_tree_cotree(mesh, dirichlet_edge_mask)` — union-find over
Dirichlet edges only, breadth-first search over the contracted graph, one
root per connected component. Invariant: tree edges = #super-nodes −
#components. The existing "disconnected PEC bodies → 2 reference groups"
test asserts the old semantics and will change; Method D's essential
incidence matrix, which is built from the groups, needs the same update. The
probe above becomes a permanent test.

---

## 7. DC (ω = 0)

- **Φ's support.** At ω = 0 equation (II) becomes `−∇·(σ∇Φ) = 0`, which is
  *identically empty* where σ = 0 — every dielectric Φ row would be zero, and
  the matrix singular. So Φ's support depends on the solve, not just the
  regime:

  | | Φ lives on |
  |---|---|
  | full-wave, ω > 0 | whole domain |
  | reduced regime, any ω | conductors (Ω_c) |
  | any regime, ω = 0 | conductors (Ω_c) |

- **Two stages.** Stage 1: Φ alone (electrokinetic), carrying every port.
  Stage 2: A alone (magnetostatic), sourced by `J = −σ∇Φ`, with the collapsed
  tree gauge. At DC the ports live entirely in stage 1.
- **Stage 2 is solvable.** Its source must be orthogonal to the gauge's null
  space; with one Dirichlet component that reduces to "net current through
  the Dirichlet boundary is zero", which stage 1 guarantees. The mixed-order
  pairing does not break this: P1 ⊂ P2, so the P2 current balance implies the
  P1 one. A discrete divergence check on the source becomes part of the
  milestone test.
- **Floating conductors.** A conductor component touching no ground and no
  port has Φ determined only up to a constant → singular. Detect with
  union-find over conductor tets; error, or treat as `pec_floating`.
- **DC current paths.** Every current port's terminal and reference must lie
  in the same conductor component (through conductors and grounded
  surfaces). The coax with its z = L end opened has no DC path through the
  via — that should be a load-time error, not a singular matrix.
- **Coax targets** (σ = 5.8e7, a = 0.1 mm, b = 0.5 mm, ℓ = 2 mm):
  - `R_DC = ℓ/(σπa²) = 1.098 mΩ`. The mesh's via cross-section is 0.985 of
    πa², so expect ≈ **1.114 mΩ** (≈ 1.5 % high). That is geometry, not solver
    error.
  - `L_DC = (μ₀/2π)(ln(b/a) + ¼)·ℓ = 0.744 nH` — external 0.644 nH plus the
    via's internal 0.100 nH, because at DC the current is uniform across the
    via. FORMULATION §7.1's `L' = (μ/2π)ln(b/a)` is the external part only:
    correct at high frequency where skin effect pushes current to the
    surface, but a DC solve compared against it would show a false 15.5 %
    error.

---

## 8. Interaction with conditioning

Port unknowns are Φ-like: the `V_k` column is a sum of Φ columns (it *is* a
Φ value) and the `V_k` row is a sum of Φ rows. So:

- **Formulation 1** (Φ rows ÷ jω) must divide the `V_k` rows **and** their
  `I_k` right-hand sides by jω too — otherwise symmetry is not restored and
  the read-back current is off by a factor jω.
- **Formulation 2** (Φ = jωΦ′) must treat `V_k` the same way: a fixed voltage
  becomes `V′_k = V_k/(jω)`, and the returned voltage is `V_k = jωV′_k`.
- The (gauge, scaling) → solver-mode lookup in `CONDITIONING.md` is
  unchanged: ports preserve the symmetry class of the Φ block.

---

## 9. Input file

Surfaces can be named by tag or by physical name; names resolve through
`$PhysicalNames` by (dimension, name), which is what the
`Mesh::physical_name` work was for.

```json
{
  "mesh": "meshes/coax_via.msh",
  "frequency_hz": 0.0,
  "regions": {
    "via_conductor": { "sigma": 5.8e7, "eps_r": 1.0, "mu_r": 1.0 },
    "dielectric":    { "sigma": 0.0,   "eps_r": 4.3, "mu_r": 1.0, "loss_tangent": 0.02 }
  },
  "boundary_conditions": {
    "shield_pec": { "type": "pec" },
    "short_end":  { "type": "pec" }
  },
  "ports": [
    {
      "name": "P1",
      "kind": "boundary",
      "excitation": "current",
      "amplitude": [1.0, 0.0],
      "terminal": "port_conductor",
      "reference": "ground",
      "sheet": "port_return"
    }
  ]
}
```

An internal voltage port — a delta-gap across the via at z = L/2, needing a
`via_cut` surface in the geometry:

```json
{ "name": "P2", "kind": "internal", "excitation": "voltage",
  "amplitude": [1.0, 0.0], "terminal": "via_cut", "direction": [0, 0, 1] }
```

A Z-matrix sweep over all ports: `"solve": { "mode": "z_matrix" }` — each
port driven in turn with `I = 1`, the others held at `I = 0`.

| field | values | notes |
|---|---|---|
| `kind` | `boundary`, `internal` | checked against the faces' tet count |
| `excitation` | `current`, `voltage` | |
| `amplitude` | `[re, im]` | peak; amperes or volts |
| `terminal` | surface name or tag | |
| `reference` | `ground`, or a surface (differential) | boundary ports |
| `sheet` | surface name or tag | boundary ports; see §3.2 |
| `direction` | 3-vector | cuts: which side is +; optional override of face winding |

---

## 10. Validation — all at load time, before assembly

**Errors:**
1. A named or tagged surface does not exist, or has the wrong dimension.
2. A `boundary` port's faces are not all boundary faces; an `internal` port's
   are not all interior.
3. A terminal shares a node with ground, with another terminal, or with a
   PEC surface (the port is shorted).
4. A cut's faces disagree on orientation, or disagree with `direction`.
5. A cut's rim touches a grounded surface.
6. No reference anywhere — no ground, no voltage port — so Φ is defined only
   up to a constant.
7. At DC: a current port with no conduction path from terminal to reference.
8. At DC: a floating conductor component (unless declared `pec_floating`).
9. Two ports use the same terminal.

**Warnings:**
- A boundary port without a sheet at ω > 0 — its voltage is path-dependent
  (§3.2).
- A cut whose rim lies inside Φ's support — delta-gap (§3.3).
- A terminal on faces with no conductor behind them — at DC it can carry no
  current.

---

## 11. Outputs

Per port: `V_k`, `I_k`, `Z_k = V_k/I_k`, `P_k = ½ Re(V_k I_k*)`. For a sweep:
the Z-matrix. S-parameters are Phase 06, from Z and a reference impedance.

---

## 12. Implementation sequence

Each step is independently buildable and testable.

1. **Reader.** Keep each triangle's original winding (needed for cut
   orientation). Rename `tagged_boundary_faces` → `tagged_faces`: cuts are
   interior. Rename coax tag 21 → `port_sheet` and regenerate the mesh.
2. **Tree-cotree rework (§6).** Collapse Dirichlet components, one tree per
   component, bodies identified along face edges. The probe becomes a test
   (nullity 0 with several Dirichlet bodies). Update Method D's group
   handling.
3. **`Problem` struct.** Materials, per-field boundary conditions, ports,
   regime, frequency, gauge, scaling. C++ only; JSON is a later adapter.
4. **Port and BC resolution + validation (§10).** Tags → faces → edges and
   nodes; boundary/interior classification from `face_tets`; + side
   flood-fill for cuts; precedence (§4); DC connectivity.
5. **`DofMap`** (§5). Tests: the coax's predicted **27,732** free A DOFs;
   edge-sign entries; two-entry cut entries; *direct reduced assembly ==
   ungauged assembly + `principal_submatrix`*.
6. **Element matrices.** Phase 04 step 3 as planned — degree-2 quadrature.
7. **Port rows, RHS, and V/I read-back.** Exact tests on a unit-conductivity
   cube: a linear potential lies in P2, so `R` is exact to round-off on any
   mesh. Boundary voltage port → `I` exact; boundary current port → `V`
   exact; the same through a cut at mid-length; an open port (`I = 0`)
   carries nothing.
8. **Coax DC milestone.** `R_DC ≈ 1.114 mΩ`; `L_DC ≈ 0.744 nH`.
9. **ω > 0.**

---

## 13. Decisions needed from you

1. **Tree-cotree rework (§6).** Changes the semantics of Phase 03 step 4 and
   one existing test. *Recommend: yes — the current design is singular with
   two Dirichlet bodies.*
2. **Port sheets.** *Recommend: required for boundary ports at ω > 0,
   optional at DC.*
3. **Cut rims inside Φ's support at full-wave.** *Recommend: allow as a
   delta-gap, with a warning.*
4. **Differential ports** (reference = another terminal). *Recommend: after
   ground-referenced ports work; the map supports it cheaply.*
5. **PEC default.** *Recommend: grounded (Φ = 0), with `pec_floating` for
   isolated conductors.*
6. **Coax tag names.** *Recommend: `port_return` → `port_sheet`, and
   `port_conductor` → `port_terminal`; regenerate the mesh.*
7. **Amplitude convention.** *Recommend: peak (matches ½ Re(VI\*)).*
