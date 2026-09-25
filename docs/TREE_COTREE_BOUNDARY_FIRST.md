# Tree-cotree construction: boundary-first, no node collapsing

*Opus 5, Sept 2026. Proposal for review, not yet implemented. Replaces the
"collapse each Dirichlet component into a super-node" wording in
`port_and_dof_plan.md` §6 with an equivalent construction that merges no nodes
at all. Both probes below were run on `APhi_Solver/meshes/cube_4.msh` with the
test matrix `M = CᵀC`.*

---

## 1. The proposal in one paragraph

Build the spanning tree over **every node of the mesh, individually** —
conductor interiors included, nothing merged. Build it in a specific order:
**first a spanning tree of each Dirichlet (`n×A = 0`) surface, using only
edges that lie on that surface; then grow into the interior from one root
surface, entering every other Dirichlet surface exactly once.** Separately, give
each equipotential surface **one Φ unknown (or none, if grounded)**, however many
nodes it has. The first part is the A gauge; the second is the Φ bookkeeping.
The two are independent.

---

## 2. Background: what the tree is for

The A-Φ system is unchanged by `A → A + ∇ψ, Φ → Φ − jωψ`: E and B do not move.
In the discrete system that relabelling freedom is an exact null space, one
direction per node. The Albanese-Rubinacci gauge removes it by setting A = 0 on
the edges of a spanning tree.

With `n×A = 0` imposed on a boundary Γ_D, the edges lying on Γ_D are already
zero. The gradients still to be removed are `∇ψ` with ψ **constant on Γ_D** —
one per node not on Γ_D. So a correct gauge adds exactly

```
N_interior tree edges = N − V_b
```

constraints, where `N` is the number of nodes and `V_b` the number of nodes on
Γ_D (for one connected Γ_D). Too few leaves a null space (singular matrix). Too
many removes physical field content (wrong answer, no error).

---

## 3. What not to do, measured

### 3.1 A plain spanning tree, no collapsing, in any order — wrong fields

`n×A = 0` on the whole boundary of `cube_4` (125 nodes, 604 edges, 98 boundary
nodes, 288 boundary edges). Correct interior tree edges: `125 − 98 = 27`.

Test: a random field z that satisfies the boundary condition (z = 0 on boundary
edges), right-hand side `j = M z`. A correct gauge must recover a field with the
same curl as z — the same B.

| tree construction | tree edges | on boundary | interior | free A | error in B |
|---|---|---|---|---|---|
| (a) plain spanning tree, started inside | 124 | 54 | 70 | 246 | **36.5 %** |
| (b) plain spanning tree, started on the boundary | 124 | 60 | 64 | 252 | **26.0 %** |
| (c) **boundary surface first, then inward** | 124 | 97 | **27** | **289** | **1.2 × 10⁻¹⁵** |

All three are spanning trees over all nodes; only (c) is correct.

**Why (a) and (b) fail.** A plain tree reaches the boundary from the inside at
many separate places. Any two of those places are also joined *along* the
boundary, where A is already zero. The result is a closed loop on which every
edge has A = 0 — which forces **zero magnetic flux through that loop**. That is a
physical constraint, not a gauge choice. The matrix stays non-singular, so
nothing fails: the answer is simply wrong. In (a) the tree places 43 more
interior constraints than a gauge needs (70 vs. 27), and B is off by 36 %.

**Why (c) is right.** The boundary's own tree is built first, from boundary
edges only, so the tree touches Γ_D as one connected piece. No loop of zero-A
edges can form through the interior, and exactly `N − V_b` interior edges are
added.

### 3.2 One root per Dirichlet surface — singular matrix

The current `build_tree_cotree` makes each PEC body its own root. Two opposite
faces of `cube_4` as separate Dirichlet surfaces; nullity of the reduced
matrix by dense elimination (dropped pivots ≈ 10⁻¹⁴, kept ≥ 0.05):

| case | tree edges | Dirichlet edges | nullity |
|---|---|---|---|
| no Dirichlet surface | 124 | 0 | 0 |
| 1 surface, its edges *not* eliminated | 100 | 0 | **24** |
| 1 surface, `n×A = 0` imposed | 100 | 56 | 0 |
| 2 surfaces, **each its own root** (current code) | 75 | 112 | **1** |
| 2 surfaces, one tree entering each surface once | 76 | 112 | 0 |

With *k* separate roots, ψ equal to a different constant on each root's tree
passes every tree and Dirichlet constraint yet is a nonzero gradient: *k* − 1
null directions. The tree is exactly one edge short (75 vs. 76).

Together, §3.1 and §3.2 bracket the correct construction: the tree must
touch each Dirichlet surface **as one piece** (else over-constrained), and the
surfaces must be joined into **one** tree (else under-constrained).

---

## 4. The construction

### 4.1 Algorithm

```
input:  mesh, dirichlet_edge[e]    (true iff edge e lies ON a Dirichlet face)

1. Find the Dirichlet surfaces: connected components of the graph that uses
   only Dirichlet edges. Record, for each node on a surface, its component.

2. For each component: breadth-first search over Dirichlet edges only, from
   any of its nodes. Every discovering edge is a surface tree edge.
   (These edges are zero anyway, by n x A = 0; recording them keeps the tree
   a genuine spanning tree over all nodes.)

3. Interior growth, from ONE root component:
     mark the root component's nodes visited; queue them all
     while the queue is not empty:
       u = pop
       for each edge e = (u, v):
         if v is visited: continue
         mark e a tree edge
         if v lies on another, not-yet-entered component:
             mark EVERY node of that component visited; queue them all
             (the component is entered through e and only e)
         else:
             mark v visited; queue v

4. If nodes remain unvisited, the mesh has another connected piece: repeat
   step 3 from a component (or any node) in that piece.
```

Step 3's "mark every node of the component" is what keeps the tree from
entering a surface twice. No node is merged: the surface's nodes keep their
own identities and their own surface-tree edges from step 2.

### 4.2 Properties to assert

- Total tree edges = `N − (number of connected mesh pieces)` — a genuine
  spanning tree (forest only across genuinely disconnected pieces).
- Interior (non-Dirichlet) tree edges = `(N − V_b) + (k − 1)` per connected
  piece, for *k* Dirichlet components: one edge discovering each node off the
  surfaces, plus one entering each non-root surface. That is `N − V_b` for a
  single boundary, and 76 for the two-face case of §3.2
  (`125 − 50 + 2 − 1`).
- Every Dirichlet component is entered by exactly one interior tree edge,
  except the root.

---

## 5. Three rules that go with it

1. **"On a Dirichlet surface" means the edge lies on a Dirichlet face** — an
   edge of a tagged face, or of a boundary face (one adjacent tet). It does
   *not* mean "both endpoints are on the surface". An interior edge can join
   two surface nodes — across a corner tet, or through a thin layer — and must
   stay a free unknown. The current code unions tagged nodes whenever *any*
   edge joins them; in a thin substrate, one tet edge from a trace to the
   ground plane below would merge trace and ground into one body: a short
   circuit, in exactly the geometry this solver targets.
2. **Conductivity never enters the tree.** Finite-conductivity regions (the
   via, traces) are ordinary nodes and edges. Their edges become tree or
   cotree edges like any others; no conductor volume is ever collapsed. σ
   affects only matrix coefficients, and Φ's support at DC.
3. **The tree gauge removes no physical content.** Where A is zero on a tree
   edge, Φ — solved at every conductor node — carries that part of the field.
   Current distribution, skin effect and losses are fully represented; the two
   gauges in `tests/test_gauge_variants.cpp` already give different A with
   identical curl.

---

## 6. The Φ side: one unknown per equipotential surface

Independent of the tree. Many nodes, at most one Φ unknown:

| surface | Φ unknowns | how |
|---|---|---|
| ground (PEC, Φ = 0) | 0 | every node fixed at 0 |
| port terminal | 1 | every node maps to the port's `V_k` |
| floating PEC | 1 | every node maps to one unknown potential, net current 0 |
| port sheet, insulating wall | one per node | Φ natural, each node its own unknown |

In the local → global DOF map (`port_and_dof_plan.md` §5.3) a terminal node is
the single entry `(V_k, 1)`, a ground node has no entry, and a free node is
`(φ_i, 1)`. This is where "equipotential" lives — not in the tree.

---

## 7. Changes to the existing code

| where | change |
|---|---|
| `build_tree_cotree` | new signature taking `dirichlet_edge_mask`; implements §4.1; drops the `is_pec` node grouping and its edge-based union |
| `TreeCotreeResult` | `node_group` / `num_reference_groups` no longer describe the tree; record instead each node's Dirichlet component (or −1) and which interior edge enters each component |
| Method D (`compute_essential_incidence_matrix`) | its fundamental-cycle walk uses the group tree; re-express it on the new tree (Dirichlet components play the role of the old reference groups) |
| `tests/test_tree_cotree.cpp` | the "disconnected PEC bodies → 2 reference groups" test asserts the old, singular behaviour and is replaced |
| new helper | `dirichlet_edge_mask(mesh, surfaces)` — edges of the chosen faces; for the whole outer boundary, edges of faces with one adjacent tet |

---

## 8. Tests

1. **B recovery** (§3.1 as a permanent test): `n×A = 0` on the whole boundary,
   random admissible z, recovered curl matches to round-off; the two plain-tree
   constructions are kept as documented negative cases.
2. **Nullity** (§3.2 as a permanent test): two and three separate Dirichlet
   surfaces, nullity 0.
3. **Counts:** interior tree edges equal `N − V_b` on `cube_4` (27) and on the
   coax.
4. **Coax free A unknowns:** exactly **27,732**
   (34,539 edges − 2,565 boundary edges − 4,242 interior tree edges; the
   boundary is a closed genus-0 surface of 1,710 triangles, so it has 2,565
   edges and 857 nodes, and `5,099 − 857 = 4,242`).
5. **No conductor collapsing:** on the coax, no interior node of the via lies
   on a surface tree, and interior tree edges pass through the via.
6. **Thin layer:** two tagged surfaces joined by a single interior edge stay
   two separate components.
7. **Existing gauge checks** (curl recovery for both gauge variants) still pass.

---

## 9. Scope and order

- **First step, sufficient for the coax:** `n×A = 0` on the whole outer
  boundary. That is one connected Dirichlet surface, so only step 2 and step 3
  of §4.1 with a single component are exercised — but built the §4.1 way from
  the start, so the multi-surface case needs no second rewrite.
- **Then:** internal PEC bodies and terminals not touching the outer boundary,
  which exercise the "enter each component exactly once" rule.
- **Deliberately out of scope:** the ungauged formulation (all edges and all
  nodes kept). Its matrix is singular, though consistent; it can work with
  Krylov iterative solvers and is a candidate for the Phase 05 iterative
  track, not for the MUMPS direct path.
