# The gauge-free scalar potential

How `Φ` is recovered from the field instead of read out of the gauged system,
why that makes it independent of the spanning tree, and — at least as
important — what it therefore cannot do.

Implementation: `src/gauge_free_potential.cpp`,
`include/aphi_solver/gauge_free_potential.hpp`.
Measurements: `GAUGE_CHOICE.md` §15.10, `tools/gauge_spike/RESULTS_SPIKEC.md`.

---

## 0. Summary

**Off by default.** `[output] gauge_free_potential = yes` turns it on.

| | |
|---|---|
| **What it is, primarily** | A **diagnostic**: comparing the recovered `Φ` against the gauged one measures how much of a terminal's potential is a spanning-tree artefact — in **one run**, where a permutation sweep needs three |
| **What it also gives** | `R` exact and tree-free, in the quasi-static regime. A by-product, not the headline |
| **What it does not fix** | `E`, `B`, `H`, `J` — unchanged inputs — and the **reactance**, which it cannot see at all |
| **Does tree–cotree still run?** | **Yes, unchanged.** This is not a constraint on the gauged system; it is a separate system solved afterwards |
| **Is it symmetric?** | **Yes**, complex-symmetric, upper triangle only, same `LDLᵀ` as the main solve |
| **Cost** | 319 ms against a 219 s factorization — 0.15 % of the run |

### The diagnostic reading

The recovered `Φ` carries no reactance (§2.4), so only the **real** parts are
comparable. Their gap is the discretization difference *plus* whatever of the
gauged terminal potential is a tree artefact:

| | gap, recovered vs gauged |
|---|---|
| case 03, both ports on the outer boundary | **0.0014 %** — discretization alone |
| case 08, one terminal inside the domain | **0.118 %** — 85× larger |

An **85× signal, not a 10¹⁰ one**, and confounded by the discretization gap: a
screening indicator that says *look harder here*, not a proof.
`07_GaugeInvariance`'s permutation harness remains the definitive test.

### Why it is off by default

Three objections, all real, none resolved by this implementation:

1. it is a **second system solved after the first**, where one united system
   with a PDE gauge — Jochum's, `JOCHUM_PORT_DERIVATION.md` — is the right
   architecture;
2. the equation is the **Coulomb/MQS** form, so the `k²` term it drops is order
   unity by 1 GHz on a 40 mm domain and dominates at THz (§7). **Quasi-static
   only**;
3. its `Φ` is **P1** where the solve it corrects uses **P2** — a lower-order
   answer replacing a higher-order one, computed from fields built with the
   higher-order one.

It will be superseded entirely by a united formulation with a PDE gauge. Until
then it is a cheap check, deliberately opt-in rather than a default.

---

## 1. The problem

`Φ` from the gauged A‑Φ system is a **gauge representative**. Tree–cotree picks
one by setting `a = 0` on a spanning tree, which fixes `ψ` at every node
relative to the root through

```
    ψ(n) − ψ(root)  =  − ∫_path A_physical · dl
```

— a **path integral**, hence tree-dependent.

Where every terminal sits on the outer boundary this does no harm: `n × A = 0`
forces `A_tangential = 0` there, which makes the boundary an equipotential for
`ψ`, and `07_GaugeInvariance` measures the terminal quantities invariant to
1e‑12 across three spanning trees.

Where a terminal is **inside** the domain, nothing pins `ψ`.
`08_MixedPort_Interior` measures the consequence: the terminal inductance moves
**3.1 %** and the nearby current density **77 %** with the choice of tree.
`GAUGE_CHOICE.md` §15.9 then derives that **no tree can fix it** — every route
to pinning `ψ` on an interior face either leaves it path-dependent, and so
unprotected, or imposes a spurious flux constraint, and so wrong.

So `Φ` has to come from somewhere other than the gauged system.

---

## 2. The formulation

### 2.1 What is actually being computed

`E` is **gauge-invariant**. It is the one thing in the solution that does not
care which tree was chosen. So: take `E` as given and ask what scalar potential
it implies.

Split `E` into an irrotational part and a `β`-solenoidal remainder — the
`β`-weighted **Helmholtz decomposition**:

```
    E  =  −∇Φ  +  E_sol ,        with    ∇·(β E_sol) = 0        (1)
```

`Φ` is exactly the irrotational potential of that decomposition. Taking
`∇·(β · )` of (1) kills `E_sol` and leaves a scalar equation for `Φ`:

```
    ∇·(β ∇Φ)  =  − ∇·(β E)                                      (2)
```

equivalently `∇·[β(E + ∇Φ)] = 0` — **find the `Φ` that makes `E + ∇Φ`
divergence-free.**

`β = σ + jωε` is the material weight, the same one the main assembly uses for
the `Φ` block. The spike used `ε_eff = ε − jσ/ω`; the two differ by a factor
`jω` that appears on both sides of (2) and cancels. `β` was chosen so the two
potentials are written in the same units and are directly comparable.

### 2.2 Weak form

Multiply (2) by a test function `λ` and integrate by parts:

```
    ⟨β ∇Φ, ∇λ⟩  =  − ⟨β E, ∇λ⟩          ∀λ                      (3)
```

which is what the code assembles. This is Stysch §3.1/§6.3
(`GAUGE_CHOICE.md` §11.8) — and the **third independent route** in that file to
a unique `Φ`, after Chew's `χ` and Ostrowski & Hiptmair's electro-quasistatic
gauge.

### 2.3 Why it is gauge-free

Not approximately, not by good conditioning — **by construction**:

> The only input to (3) is `E`. `E` does not depend on the gauge. The boundary
> conditions do not depend on the gauge. Therefore `Φ` does not depend on the
> gauge.

**No tree appears anywhere in (3).** Contrast every other remedy tried: the
tree contraction of §15.9 changed *which* tree, and bought 10× where nine
orders were needed. This removes the tree from the question.

Relation to the gauged potential: with `E = −jωA − ∇Φ_gauged`, (1) gives

```
    ∇Φ_recovered  =  ∇Φ_gauged  +  (jω A)_irrotational
```

so the two differ by exactly the irrotational part of `jωA` — which is the
gauge-dependent piece, removed.

### 2.4 Why there is no reactance, and why that is the same fact

(1) discards `E_sol`. The inductive contribution `−jωA` is essentially
solenoidal, so the decomposition assigns it to `E_sol` and it never reaches
`Φ`.

Measured: `Im(Φ) ≈ 1e‑16` in **every** run — including case 03, where the
gauged answer is entirely correct and carries `L = 22.6 nH`.

This is not a defect to be fixed in a later version. **It is the same property
that makes the method gauge-free.** The tree-dependence lives in the solenoidal
part; `div` annihilates the solenoidal part; so the method is blind to the
tree-dependence *and* to the inductance, together, for one reason.

A consequence worth stating plainly, because it bounds what any post-process
could achieve: the three tree-permuted solutions differ by a field with
`div(β ΔE) = 0`. **Any** method whose only view of the solution is through a
divergence is blind to that difference. Two independent attempts confirmed it —
the §15.9 tree contraction (10×) and an energy-based inductance (**42 %**
spread across three trees). The post-processing avenue is closed; correct
*fields* at a mixed port require the gauge itself to change.

---

### 2.5 Provenance — the equation is published, the usage is ours

Worth separating, because the two carry different amounts of warranty.

**The equation is Stysch's (6.25)**, verbatim up to two terms we do not have:

> *"Find `φ_c ∈ H¹(Ω)` such that*
> `⟨ε_r grad φ_c, grad ψ⟩ + (2/r)⟨φ_c, ψ⟩_Γa = −⟨ε_r E, grad ψ⟩ + jμ₀ck⟨g, ψ⟩`*"*

| his (6.25) | ours, (3) |
|---|---|
| `⟨ε_r ∇φ_c, ∇ψ⟩` | `⟨β ∇Φ, ∇λ⟩` |
| `−⟨ε_r E, ∇ψ⟩` | `−⟨β E, ∇λ⟩` |
| `(2/r)⟨φ_c, ψ⟩_Γa` | absent — we have no absorbing boundary |
| `jμ₀ck⟨g, ψ⟩` | absent — his explicit source current; ours is already inside `E` |

His (6.26) is the Lorenz-gauged variant, carrying an extra `−k²⟨φ_c, ψ⟩` mass
term; (3) is his Coulomb form, which is the MQS limit of it. He also states the
well-posedness this relies on: *"The weak-form BVP (6.25) does not possess a
general LF breakdown. If the boundary is not purely magnetic… the resulting
system matrix always has full rank."*

**Three independent sources converge on this same scalar problem** — Stysch
(6.25), Ostrowski & Hiptmair's electro-quasistatic gauge, and Chew's `α → 0`
limit, which `RESULTS_SPIKEC.md` milestone 3 measured degenerating to exactly
it. That convergence is the main reason to trust the equation.

**The usage is not published, and that is the part to be careful with.**
Stysch never solves A‑Φ at all: his §3.1 explicitly chooses the "E approach"
over the "A‑Φ approach", so (6.25) is his formulation's **only** `Φ` — primary,
not recovered. We instead solve a tree-gauged A‑Φ system that *already produces
a `Φ`*, then recover a **second** one from the resulting `E`, and tell the
reader to take `R` from the second and `L` from the first. That is his equation
applied as a **repair to another formulation's output**, and it has not been
found reported anywhere.

Two consequences follow, and neither is resolved:

- **No joint consistency guarantee.** His `Φ` and his `E` come out of one
  formulation and satisfy one set of equations by construction. Ours come from
  two. Nothing guarantees that `Φ_recovered` and the gauged `E` together
  satisfy any single system. The 0.0014 % agreement on case 03 (§5.1) is
  reassuring, but it is one number on one geometry.
- **It invites the obvious objection.** If a separate scalar solve is needed
  anyway, why tree-gauge the first one? That *is* Stysch's argument, and it
  points at restructuring the formulation rather than patching it. The absence
  of this trick from the literature may mean it is an unremarkable shortcut, or
  may mean the field went the other way for a reason. **We do not know which.**

So: the formulation is published and well grounded; the application is ours and
is tested on two cases.

---

## 3. Discretization

### 3.1 Spaces

`Φ` is **P1 on mesh vertices**, deliberately, although the main solve uses P2.
The source `div(βE)` is piecewise constant on a tet, so P2 would double the
system and buy nothing. The price is visible and small: on case 03 the
recovered and gauged real parts differ by **0.0014 %**, which is that
discretization gap and not a gauge effect.

`E` enters as `fields.e_tet[t]` — one complex vector per tet, the value at the
centroid, already computed by `compute_fields`.

### 3.2 Element matrices

Per tet, with `β` from the tet's body and `V` its volume:

```
    stiffness     K_ij  =  β · V · (∇λ_i · ∇λ_j)
    load          b_i   =  − β · V · (E · ∇λ_i)
```

Both exact: `∇λ` is constant on a tet and `E` is taken constant on it, so no
quadrature rule is involved.

A tet with `β == 0` is skipped. At DC that is every insulator, which is correct
— `Φ` has no support there.

### 3.3 The system is complex-symmetric

`K_ij = β V (∇λ_i · ∇λ_j)` is manifestly symmetric in `i, j`: `β` is a scalar
per tet and the dot product commutes. So:

- only the **upper triangle** is stored, in `SparseSymmetricZ`. The assembly
  skips `cj < ri`, and `from_pattern` **throws** on a below-diagonal entry, so
  an asymmetric assembly could not go unnoticed;
- it is solved with `solve_symmetric`, the project's existing complex-symmetric
  `LDLᵀ` — the same factorization the main system uses, not a second solver.

Measured residual: 1.5e‑13 on case 08, 2.3e‑13 on case 03.

### 3.4 Terminals

`build_vertex_map` mirrors the main solve's rules exactly, so that the two
potentials are comparable rather than two conventions:

| vertex | treatment |
|---|---|
| on a **current** port | shares that port's **single** unknown — this is what makes a terminal a terminal |
| on a **voltage** port | prescribed at that port's value; leaves the system |
| on a **cut** | prescribed zero. `Φ` *jumps* across a cut, and single-valued P1 cannot represent that, so the vertex is excluded rather than quietly averaged across the jump |
| absent / fixed | eliminated, as in the main solve |
| otherwise | free |

An eliminated column's known value is lifted to the right-hand side in the
usual way (`rhs -= k * val`).

### 3.5 The reference — what fixes the additive constant

Only `∇Φ` appears in (3), so **`Φ` is determined only up to a constant.** What
pins it is the **voltage port**: its vertices are prescribed, and the Dirichlet
lift carries that value into the system.

In case 08 that is `P2`, `boundary_voltage` with `voltage = 0.0` on
`wire_bottom`, and the output reads back `port P2 Phi = 0 + 0j V` — exactly as
prescribed, which is a free check that the elimination is right.

A reference always exists: `input_file.cpp:809` **refuses** a file whose ports
are all boundary current sources, precisely because `Φ` would then be fixed
only up to a constant. So (3) is never a pure Neumann problem in practice.

Crucially it is the **same** reference the main solve uses, which is what makes
`Φ(P1) − Φ(P2)` directly comparable with the gauged terminal voltage.

---

## 4. Where it sits relative to the main solve

**Tree–cotree still runs, unchanged.** Nothing in the gauged system was
modified, no constraint was added to it, and its matrix, its gauge and its
answers for `A`, `E`, `B`, `H` and `J` are exactly what they were.

```
    ┌─ solve 1 ─────────────────────────────────┐
    │  the gauged A-Φ system, tree-cotree        │
    │  → A, Φ_gauged → E, B, H, J                │   unchanged
    └──────────────┬─────────────────────────────┘
                   │  E  (gauge-invariant)
    ┌──────────────▼─────────────────────────────┐
    │  solve 2: ⟨β∇Φ, ∇λ⟩ = −⟨βE, ∇λ⟩            │   new
    │  → Φ_gaugefree  →  R                       │
    └────────────────────────────────────────────┘
```

Two independent linear systems. The second reads one output of the first and
feeds nothing back. It can be removed without affecting anything else, and its
failure cannot corrupt the main solve.

---

## 5. Measured

### 5.1 Validation first, on the case whose answer is already right

Case 03, both ports on the outer boundary, where the gauged answer is correct
and a disagreement would mean the implementation is wrong:

| | gauge-free `Φ` | gauged `Φ` |
|---|---|---|
| Re | 9.79686e‑05 | 9.796995e‑05 — **agree to 0.0014 %** |
| Im | **−3.0e‑16 ≈ 0** | +7.0986e‑06 (`L` = 22.6 nH) |

### 5.2 On the mixed interior port

Case 08, three spanning trees, against the analytic
`R = L/(σ·A_poly) = 6.857804222298e‑05 Ω`:

| | gauged `Φ` | gauge-free `Φ` |
|---|---|---|
| base | 68.659272 µΩ | **68.57804222150 µΩ** |
| permA | 68.660045 | 68.57804222146 |
| permB | 68.606254 | 68.57804222188 |
| **error vs exact** | **0.12 % high** | **2.2e‑11** |
| **tree spread** | 7.7e‑04 on `R`, **3.1e‑02 on `L`** | **5.5e‑12** |
| Im | +6.4e‑06 | −1.0e‑16 |

About **eight orders better** on both accuracy and tree-stability.

*On the eleven figures:* a uniform bar has a `Φ` linear in `z`, which P1
represents **exactly**, so this geometry flatters the accuracy. Expect ordinary
discretization error on a non-uniform conductor. **The tree-invariance is the
geometry-independent part, and that is the result.**

---

## 6. How to use it, and how not to

| quantity at an interior mixed port | where to read it |
|---|---|
| **`R`**, resistive terminal drop | **`potential_gaugefree.out`** — exact, gauge-free |
| `L`, reactance | the gauged `potential.out` — **and it is wrong by 3.1 %** |
| `B`, `H`, local `J` | the field files — **and they are wrong by 50–91 %** |

This is **not** a replacement for the gauged potential and must not be read as
"mixed ports are fixed". It is a complementary extraction, reliable precisely
where the gauged one is not. Both file headers carry the same rule: **read `R`
from the gauge-free file, read `L` from the gauged one, and know what the
latter is worth.**

Its primary use is the one line it prints per port:

```
    port 'P1'  recovered Re(Phi) = 6.8578e-05 V,  gauged 6.86593e-05 V,  gap 0.118 %
```

A gap at the discretization level means the terminals are where the gauge can
protect them. A gap orders above it means at least one terminal is not, and the
permutation harness should be run.

---

## 7. Limits and known rough edges

- **No reactance, ever.** §2.4. Structural, not a gap to be filled later.
- **P1 against the solver's P2**, so the recovered `Φ` is one order lower than
  the gauged one. Visible as the 0.0014 % gap on case 03.
- **`E` at tet centroids** is first-order. On case 08 it still gave eleven
  figures, but that is the uniform-bar effect of §5.2 and should not be
  expected generally. A refinement study has **not** been done.
- **A cut's `Φ` is excluded, not recovered.** `Φ` jumps across an internal cut
  and single-valued P1 cannot represent the jump, so those vertices are
  prescribed zero. A case driven only by an internal cut therefore gets nothing
  useful from this file. Extending it would mean a two-sided recovery.
- **`PhiDof::Absent` is routed through `prescribed[] = 0`**, which is *not* the
  same as "has no equation". It is currently safe only because the assembly
  skips `β == 0` tets, and absent nodes occur only where `β == 0` in both
  regimes we run. **That is correct by accident rather than by design** and
  should be made explicit — either an assertion, or a branch that genuinely
  excludes them.
- **Not asserted in any regression test yet.** The numbers above are measured
  and recorded, not pinned.

---

## 8. Bibliography

The sources this rests on, all in `APhi_Papers/`, and what each one supplies.
`GAUGE_CHOICE.md` carries the full reading of each.

| | |
|---|---|
| **Stysch (2022)**, *Stable Broadband Finite Element Parasitic Extraction and Sensitivity Analysis*, PhD thesis, TU Darmstadt, open access (handle `tuda/8889`) — `Stysch2022_PhDThesis_StableBroadbandParasiticExtraction.pdf` | **The equation.** §6.3's (6.25) is (3) up to an absorbing-boundary term and an explicit source current; (6.26) is its Lorenz-gauged variant. §6.3 also states the full-rank property (3) relies on. §3.1 is his case for solving `E` first and `Φ` as a post-process — and, read carefully, an argument against the way we use it (§2.5). |
| **Stysch, Klaedtke & De Gersem (2022)**, *Broadband finite-element impedance computation for parasitic extraction*, Electrical Engineering **104**(2) 855–867, [doi:10.1007/s00202-021-01348-9](https://doi.org/10.1007/s00202-021-01348-9), preprint [arXiv:2009.08232](https://arxiv.org/abs/2009.08232) | The article form. Its §2 is the clearest published statement of why a terminal potential is gauge-dependent at all: `V = V_c + jωI·L_r` with `V_c = Φ(T_b) − Φ(T_a)`, and *"as it incorporates a partial inductance, `V_c` formally depends on the gauge condition."* That is the sentence this whole document is downstream of. |
| **Ostrowski & Hiptmair (2021)**, *Frequency-Stable Full Maxwell in Electro-Quasistatic Gauge*, SIAM J. Sci. Comput. **43**(4) B1008–B1028 (ETH SAM report 2020‑43) — `ETH2020_FreqStableFullMaxwell_EQSGauge.pdf` | The same scalar problem reached independently, as a **gauge condition** rather than a post-process: the electro-quasistatic problem used to fix the potentials in a two-step scheme. |
| **Chew (2014)**, *Vector Potential Electromagnetics with Generalized Gauge for Inhomogeneous Media: Formulation*, PIER **149** 69–84 | The third convergence. His decoupled `Φ` equation (2) with `χ = αμε²` degenerates to (3) as `α → 0`, measured in `RESULTS_SPIKEC.md` milestone 3 — where it also shows why `α ≠ 0` is unusable for a lossy conductor. |
| **Rapetti, Alonso Rodríguez & De Los Santos (2022)**, *J* **5**(1) 52–63 — `Rapetti2022_OnTheTreeGauge_J5.pdf` | Why the gauged `Φ` has no continuum limit to begin with: tree gauges *"are not a discretization of the Coulomb gauge and enforce no orthogonality"*, and the tree degrees of freedom are *"set arbitrarily"*. The reason §1's problem is structural rather than a matter of mesh quality. |
| **Ansari, Farquharson & MacLachlan (2017)**, GJI **210**(1) 105–129 — `Ansari2017_GaugedFE_CoulombLagrange_GJI210.pdf` | The galvanic/inductive decomposition of `E`, which is the same split (1) performs — `Φ` carries the galvanic part, `A` the inductive. He obtains unique potentials by an explicit Lagrange-multiplier Coulomb gauge rather than by recovery, so the route differs; the decomposition does not. |
| **COMSOL**, *Computing and Using Partial Inductance with COMSOL®* (vendor documentation) | Independent public confirmation of the limit in §6: total inductance comes from stored energy, `L_tot = 2W_m/I²`, while partial inductance is a free choice of subdivision and *"we can never measure any of these partial inductances, as only the total inductance of a closed loop is measurable."* |

**Not found reported anywhere:** recovering `Φ` this way as a *repair* applied
to a tree-gauged A‑Φ solve, which is what this implementation does. §2.5.
