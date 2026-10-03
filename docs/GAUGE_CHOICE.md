# Choosing a gauge

A decision report for the gauge that replaces tree–cotree, written against the
primary sources rather than summaries of them.

**Second edition.** The first draft recommended Zhao & Fu's dummy-variable
Coulomb gauge on cost-and-fit grounds. That recommendation is **withdrawn**. It
was made without the potential-BEM endgame in view and before the primary
sources had been obtained; both change the answer.

---

## 0. Decision

| | |
|---|---|
| **Recommended gauge** | **Chew's generalized gauge** — `∇·(εA) = −χ ∂Φ/∂t`, `χ = αε²μ` |
| **Why, in one line** | It is the only candidate under which the FEM half, the interface conditions, and the BEM half all come from one formulation. |
| **Is anything broken today?** | **Not for any geometry in the suite — and every one of them puts its terminals on the outer boundary.** `E`, `B`, `H`, `J`, terminal `V`, `I`, `R`, `L`, `Z` are correct and gauge-invariant there, measured to 1e‑12 by `07_GaugeInvariance`. That result does **not** extend to a port inside the domain: see §1.2. |
| **What is wrong** | Φ in the interior is not physical; tree–cotree postpones rather than removes low-frequency breakdown; and **at a port face inside the domain — especially one crossing conductor and dielectric — it has no principled rule at all, at any frequency**. |
| **What it costs us** | Less than the first edition assumed: **our Whitney‑1 `A` space is already the one Chew's gauge requires.** Φ moves P2 → P1 and the gauge term is added; the `A` discretization is untouched. |
| **Still open** | Whether `A` stays a Whitney 1‑form **at the boundary-integral surface**. This is now the real research question — see §6. |
| **Who decides** | The mixed-port benchmark the roadmap requires. This report narrows the field to one candidate; it does not retire the benchmark. |

---

## 1. Where we stand

Tree–cotree sets `A = 0` on a spanning tree's edges. That makes the matrix
non-singular and nothing more. `07_GaugeInvariance` measures the consequence on
one problem solved on three different trees:

| | change across trees |
|---|---|
| `Φ`, `A` | **0.14 – 0.22 relative** |
| `E`, `B`, `H`, `J` | 4.7e‑13 … 3.2e‑12 |
| terminal `V`, `R`, `L` | 2.7e‑13 |

`docs/FIELD_POSTPROCESSING.md` §9 quantifies how much of Φ is gauge: about
`ωL/|Z| = sin(arg Z)`. Case 06 sits at `ωL/R = 4.48`, case 04 at 12.8 — both
gauge-dominated, which is why neither asserts Φ.

**A consequence worth stating before anything else.** Our terminal quantities are
gauge-independent *only* because `n × A = 0` on the outer boundary forces the
gauge function ψ to be constant there — measured, `|ψ|` is 3.6e‑16 on the
boundary against 4.6e‑04 inside. `n × A = 0` is a condition on **A**, so it is
not gauge-invariant; that is exactly why it pins ψ.

**Replace it with a radiation, transmission or boundary-integral condition — all
of which are written in the gauge-invariant fields — and ψ stops being pinned at
all.** The moment the boundary treatment becomes correct for THz, tree–cotree
leaves not only interior Φ but the terminal quantities floating. The gauge must
become a genuine PDE condition *before* the boundary condition changes, not
after.

### 1.1 Why tree–cotree's Φ has no continuum limit — Rapetti et al. settle it

F. Rapetti, A. Alonso Rodríguez and E. De Los Santos, *On the Tree Gauge in
Magnetostatics*, **J 5(1), 52–63, 2022**, DOI 10.3390/j5010004.

The project's proposal attributes to this paper the claim that the tree gauge is
not a discretization of any orthogonality condition. **The claim is correct, and
the paper states it outright.** Their §4, immediately after defining the tree
gauge `a_t = 0`:

> It is worth noting that **these tree gauges are not a discretization of the
> Coulomb gauge** stated in (5) or (6). **They are not enforcing in any sense the
> orthogonality to the gradient** in conditions (5) and (6).

And their §5 opens by saying the tree degrees of freedom

> **is set arbitrarily**, eventually equal to zero, without affecting the
> corresponding field `B_h = curl A_h`.

That is case 07's result stated in the literature — for `A`, in magnetostatics,
from the algebraic side rather than by measurement.

**This is the deepest available answer to "why is tree–cotree's Φ not a physical
quantity."** The discrete Coulomb gauge *is* an orthogonality condition — their
eq. (11) says `A_h` is orthogonal to `grad(W⁰)`, and their §5 constructs it as a
projection onto `(ker S)^⊥ = Im(S)`. The tree gauge is not the discretization of
anything in the continuum; it is a choice of which algebraic degrees of freedom
to delete. A quantity fixed by such a choice **cannot converge to a definite
function under mesh refinement, because there is no continuum condition for it to
converge to.**

That is the structural reason behind every measurement in this report, and it is
why §9's first test — pointwise mesh convergence of Φ — is the one that matters.

### 1.2 The limit of what case 07 proves — interior and mixed ports

Case 07's result is **conditional on where the ports are**, and the condition is
easy to lose sight of because every case in the suite satisfies it.

The terminal quantities are gauge-independent because `n × A = 0` on the outer
boundary forces ψ constant there and a port pins it to zero. **Both terminals of
every case in this suite lie in the outer boundary** — the conductors span the
full height of their air boxes, by construction. Measured: `|ψ|` is 3.6e‑16 on
the boundary against 4.6e‑04 in the interior, twelve orders apart.

**A port face inside the domain has none of that protection.** ψ is unpinned
there, so Φ on such a face is gauge-dependent — which matters not only for
reading a voltage out but for *posing* a current-driven port, whose terminal
potential is a free unknown.

A face crossing **conductor and dielectric** is worse than merely interior. The
tree gauge grounds tree paths to conductors; on a face that is partly conductor
and partly dielectric there is no principled rule for which edges to take, and —
by Rapetti et al. (§1.1) — **no continuum condition to appeal to**, because the
tree gauge is not a discretization of one. There is nothing to be right or wrong
against.

> **This is not a high-frequency problem.** Ansari et al. (§3.1) demonstrate
> non-unique potentials in the **quasi-static** regime: a conductive prism in a
> conductive background, driven by a grounded wire — a galvanic low-frequency
> geophysics problem with material contrasts inside the domain. Their fields were
> fine and their potentials were not. **The failure mode is about material
> contrast and port placement, not about frequency.**

**What this project has and has not measured.** We have measured, thoroughly,
that tree–cotree is sound for boundary-mounted ports. **We have never built a
case with an interior or mixed-material port, so we have not measured the failure
either** — it is predicted by the mechanism above and documented in the
literature, not yet observed here. §9 proposes closing that gap, and it should be
closed *before* the gauge change, so that there is a failing baseline to compare
against.

### 1.3 Measured: an interior cut's LOCAL field is tree-dependent at 8 %

The strongest evidence in this report is not from the literature. It is a
measurement on `05_Loop_1A_50Hz`, a case that **passes all its assertions**, at
**50 Hz**, through a cut that lies **entirely inside a conductor** — the most
favourable interior-port configuration there is, with one side of the cut even
pinned at 0 V.

Case 05's ring is driven by an internal cut and asserts `L_nH = 19.7854 ± 0.05`,
obtained from `Z = V / I_drive`, so the assertion rides on a potential. Running
`07_GaugeInvariance`'s permutation harness on it, three spanning trees:

| | base → permA | base → permB |
|---|---|---|
| `L` (asserted 19.7854 ± 0.05 nH) | **+0.0035 nH** | +0.0009 nH |
| `\|V\|` | ~3e‑6 relative | ~3e‑6 relative |
| **`J` within 2 mm of the cut** | **8.3 %** | — |
| **`B` within 2 mm of the cut** | **7.4 %** | — |

and the field change decays sharply away from the cut — `J` goes
**8.3 % → 2.3e‑4 → 7.2e‑6** across the 0–2, 2–5 and 5–10 mm bands. It is
localized at the port, not global.

**Three checks say it is real rather than numerical noise:**

- all three solves are exact — backward errors 9.1e‑23, 9.0e‑23, 9.2e‑23;
- `ΔE = ΔJ/σ` to a ratio of **1.0000**, so the `E` and `J` changes are one
  physical change reported against different global maxima;
- it decays smoothly over ~11,000 nodes, so it is not the 205 cut nodes being
  mislabelled in the output.

**The integrated quantities survive and the local ones do not.** `L` moves by
0.018 %, comfortably inside a tolerance it was never in danger of failing — which
is exactly why case 05 passes and why nothing flagged this. **Case 05 is not
wrong.** What is wrong is that the current distribution within a couple of
millimetres of the feed depends on an arbitrary spanning tree, and current
distribution near a feed is precisely what a phased-array analysis is for.

> **Practitioner experience, recorded as corroboration.** From the project
> owner's own use of a commercial A‑Φ solver: a mixed-material port ran to
> completion without crashing and returned wrong solutions, and a Coulomb-type
> gauge improved matters without fully fixing them. The working assumption there
> was that **internal ports lying inside a conductor were fine under tree–cotree,
> and only mixed-material interior ports were problematic.** (Observed
> user-visible behaviour only; nothing here is a claim about that tool's
> internals.)
>
> **This measurement refines that assumption rather than contradicting it.** For
> the quantities such a port is normally asked for — terminal voltage, R, L, Z —
> the assumption holds, and holds well: 0.018 % on `L`. It is the **local field**
> that moves, by 8 %, and no integrated check would ever reveal it. So the
> conductor-only interior port is not "fine"; it is **fine for circuit
> quantities and unreliable for field quantities**, and the distinction has
> probably gone unnoticed because the things people assert about such ports are
> all integrals.

**Why this is the report's most useful result.** It is ours, it is measured, it
is at low frequency, and it is on a case that passes — so it cannot be dismissed
as a pathological configuration or a high-frequency effect. It also sets the bar
for the replacement: a new gauge has to make the local field near an interior
port tree-independent, not merely keep `L` inside tolerance.

---

## 2. Chew 2014 is the keystone

W. C. Chew, *Vector Potential Electromagnetics with Generalized Gauge for
Inhomogeneous Media: Formulation* (Invited Paper), **Progress In Electromagnetics
Research, Vol. 149, pp. 69–84, 2014**, DOI 10.2528/PIER14060904. Open access;
now in `APhi_Papers/Chew2014_GeneralizedGauge_PIER149.pdf`.

### 2.1 The gauge

His eq. (15), with (18):

```
    ∇·(εA) = −χ ∂Φ/∂t        χ = α ε² μ
```

where **α may be a function of position**. `α = 1` recovers the generalized
Lorenz gauge `ε⁻¹∇·(εA) = −με ∂Φ/∂t` (his eq. 14); a homogeneous medium with
`α = 1` recovers the ordinary Lorenz gauge.

> **`χ = αε²μ` is confirmed from the source.** The first edition of this report
> flagged the roadmap's rendering as unverified because it had been read out of
> garbled PDF text. It is correct.

The resulting equations (16) and (17) decouple, and (16) is derivable from (17)
by taking the divergence and using charge continuity. Chew states plainly that,
unlike the vector wave equations for the fields, **these do not exhibit the
low-frequency catastrophe**, so ordinary FEM/FD solvers work across scales and
bandwidth.

He also rules out the alternative explicitly: under the **Coulomb gauge** the
scalar potential has infinite velocity and *the vector potential equation is not
completely decoupled from the scalar potential equation*. That single sentence
disqualifies every Coulomb-type candidate for a potential-BEM hybrid, for the
reason developed in §4.

### 2.2 It dissolves the port problem — which is why you went looking

§3, on the interface between two media:

> the boundary conditions at the interface of two homogeneous media are **also
> embedded in these equations**. Therefore, when one solves the PDEs directly in
> an inhomogeneous medium, **one need not stipulate the boundary conditions**.
> The solutions naturally obey the boundary conditions if they are arrived at
> correctly via a numerical method.

The conditions he then derives are the classical ones:

```
    n̂ × A₁ = n̂ × A₂                              (23)   tangential A continuous
    n̂ × (1/μ₁)∇×A₁ = n̂ × (1/μ₂)∇×A₂   ≡  n̂ × H₁ = n̂ × H₂     (24),(25)
    n̂ · (ε₁A₁) = n̂ · (ε₂A₂)                      (26)   normal εA continuous
    n̂ · (ε₁∇Φ₁) = n̂ · (ε₂∇Φ₂)                    (27)
```

together giving `n̂·(ε₁E₁) = n̂·(ε₂E₂)`, the usual normal-`D` condition.

**Equation (26) is the one that decides our element choice**, and §6 is about why:
it is `εA`, not `A`, whose normal component is continuous, so **`n̂·A` jumps
wherever ε jumps.**

**So a port face that crosses conductor and dielectric needs no special
treatment.** The gauge carries the interface physics. `ROADMAP.md` Phase 03
identifies exactly this as the reason tree–cotree is inadequate; this is the
paper that answers it.

He also treats the PEC limit: as `ε₂ → ∞`, `n̂ × ∇Φ₂ = 0`, so `Φ₂` is constant in
the PEC and `n̂ × E₁ = 0` **even at ω = 0**.

### 2.3 It already contains the boundary-integral half

§6 derives the **generalized Green's theorem, the extinction theorem and the
surface equivalence principle** for the A‑Φ formulation; §7 derives a **surface
integral equation for a PEC scatterer** with its matrix representation. Per the
abstract, those integral equations exhibit **neither the low-frequency
catastrophe nor the frequency imbalance** of the classical `E`–`H` ones.

**This is the decisive architectural fact: the FEM half, the interface
conditions, and the BEM half are one formulation under one gauge.** No other
candidate offers that.

Two implementation details from §6–7 worth recording now:

- For the SIE derivation Chew assumes the **simple Lorenz gauge**, because SIEs
  apply to piecewise-*homogeneous* regions where ε is constant. The generalized
  gauge is what the inhomogeneous FEM volume needs; the simple one is what the
  homogeneous exterior needs. They agree where it matters.
- The equivalence principle needs **six scalar quantities** on the surface: two
  tangential components of `H₁`, two tangential components of `A₁`, `Φ`, and
  `n̂·A₁`. He notes this parallels the augmented EPA, where six components are
  likewise needed to keep low-frequency stability.
- The surface current basis must be **divergence-conforming** (RWG-type) so that
  the `A` it produces is divergence-conforming too.

---

## 3. A published statement of the problem case 07 measures

S. Sharma and P. Triverio, *Electromagnetic Modeling of Lossy Interconnects From
DC to High Frequencies With a Potential-Based Boundary Element Formulation*,
IEEE Trans. Electromagn. Compat. **64**(4), 2022 — [arXiv:2112.07360](https://arxiv.org/abs/2112.07360).

They use the **usual Lorenz gauge outside** and a **modified Lorenz gauge inside**
a lossy conductor:

```
    ∇·Ã = −jωε₀μ₀ Φ                 in free space        (17)
    ∇·Ã = −(jωε + σ) μ Φ_r          in the conductor     (19)
    Φ_r = Φ − Φ_a ,   Φ_a = average of Φ over the object's surface   (20),(21)
```

Their stated reason is the one this project measured independently:

> in order to determine Φ uniquely within V, one must specify a boundary
> condition on S and a reference point; if the conventional Lorenz gauge is
> used, **a change in reference for Φ will lead to a change in `∇·Ã|_V`** — in
> other words `∇·Ã|_V` cannot be determined uniquely until a reference is set
> for Φ.

Subtracting the object's own average surface potential makes the interior gauge
reference-independent, and both regions then yield clean Helmholtz equations.

**This is a close published analogue to case 07's result** — the potential
*reference* contaminating the vector potential, where ours is the *tree* doing
it. The disease is the same: a bookkeeping choice leaking into a solved
quantity. It is worth citing in the case 07 README.

### 3.1 Ansari et al. is the closest analogue of all

M. Ansari, C. G. Farquharson and S. P. MacLachlan, *A gauged finite-element
potential formulation for accurate inductive and galvanic modelling of 3-D
electromagnetic problems*, **Geophys. J. Int. 210(1), 105–129, 2017**,
DOI 10.1093/gji/ggx149. Paywalled; abstract and technical content read from the
publisher's article page, PDF not yet obtained.

They demonstrate **"non-unique, incorrect potentials"** from the usual
incompletely-gauged system — and the demonstration is the striking part. They
show it by

> showing **inconsistent results obtained from iterative and direct linear
> equation solvers**

— the same physical fields, but **different potentials depending on which linear
solver was used**. After introducing the Coulomb gauge explicitly as an extra
equation and augmenting the Helmholtz equation with the gradient of a Lagrange
multiplier, "both the iterative and direct solvers produce the same responses for
the potentials, demonstrating the uniqueness of the numerical solution."

Two further points make this the most directly transferable result in the
collection:

- **They use edge elements for `A` and nodal elements for Φ — our exact
  discretization.**
- They identify the root cause as the **normal component of `A` being
  discontinuous across material interfaces**, which edge elements permit and
  which an incompletely-gauged system therefore leaves unconstrained.

**Three independent triggers, one disease.** Collecting what the literature and
this project have now each observed:

| observed by | the potentials change with… | fields unchanged? |
|---|---|---|
| Ansari et al. (2017) | the **linear solver** (iterative vs direct) | yes |
| Li, Sun, Dai & Chew (2015) | the **initial guess** of an iterative solve | yes |
| Sharma & Triverio (2022) | the **potential reference** inside a conductor | yes |
| **`07_GaugeInvariance`** | the **spanning tree** | yes, to 1e‑12 |

Each trigger is different and each is incidental; the common cause is that
nothing in the formulation pins the gauge.

Note also that they *exploit* gauge freedom rather than merely tolerating it:
they leverage the gauge invariance of the potentials to devise simpler boundary
conditions for `A` than existing potential-integral-equation methods use.

---

## 4. Why a potential-BEM forces a Lorenz-type gauge

This is physics, independent of any citation.

In the **Lorenz** gauge the exterior potentials satisfy

```
    (∇² + k₀²) A = 0        (∇² + k₀²) Φ = 0
```

so **every Cartesian component of `A`, and `Φ`, satisfies the same scalar
Helmholtz equation** — which is exactly what admits the scalar Green's function
`G₀` and a boundary representation built from standard single- and double-layer
operators.

In the **Coulomb** gauge this collapses: Φ satisfies Poisson's equation
(instantaneous, no `k₀`) and `A`'s equation retains a `∇(∂Φ/∂t)` coupling term —
which is precisely Chew's objection in §2.1. The clean scalar-Green's-function
BEM is simply unavailable.

Corroborated by every potential-BEM paper obtained:

| | gauge used |
|---|---|
| Chew 2014 §6–7 | simple Lorenz for the homogeneous regions |
| Sharma & Triverio 2022 | Lorenz outside, modified Lorenz inside |
| Li, Fu & Shanker 2017 | Lorenz; decoupling is what gives low-frequency stability |
| Abdrabou & Gomez 2026 | Lorenz, stated in the abstract |

---

## 5. Do hybrids still need gauging? Abdrabou & Gomez answer it

A. Abdrabou and L. J. Gomez, *A Hybrid DEC-SIE Framework for Potential-Based
Electromagnetic Analysis of Heterogeneous Media*, J. Comput. Phys. 2026 —
[arXiv:2507.02099](https://arxiv.org/abs/2507.02099).

**They use no tree–cotree and no discrete gauging at all.** The gauge is imposed
in continuous form before discretization, and they represent **`A` as three
discrete 0‑forms (nodal scalars)** rather than the conventional 1‑form. Their own
framing: this sidesteps gauge-freedom issues at the discrete level by treating
the components as independent scalars that satisfy the Lorenz constraint
identically.

**The null space you gauge away is a property of the 1‑form / edge-element
representation.** Drop it and there is nothing to gauge.

That is the honest answer to "are gauges still needed in a hybrid": **not if you
abandon Whitney 1‑forms for `A`. If you keep them — and we do — then yes.**

Their scalar reformulation reduces the surface integral operators **from fourteen
to two** (single- and double-layer). The reduction is performed by expressing the
problem in Cartesian components and their normal derivatives; **the Lorenz gauge
is the necessary enabler**, because it is what makes each component satisfy the
same decoupled scalar equation in the first place. Gauge and reformulation are
two steps, not one.

> **A constraint their method carries — and how they dispose of it.** The
> coupling surface Γ must lie entirely in **free space** for the component-wise
> continuity of `A` and Φ to hold naturally: "the surface Γ has to be only in
> free space." Their fix is a **buffer region**: the DEC domain is the
> inhomogeneous body *plus* a homogeneous free-space collar, `Ω₂ = Ω_s ∪ Ω_a`,
> so that `Γ = ∂Ω₂` is embedded entirely in free space, which in their words
> greatly simplifies enforcing the boundary conditions.
>
> **This is directly transferable.** A THz array on a substrate puts the
> boundary-integral surface outside the substrate with a free-space collar,
> rather than cutting through it. It costs volume elements in the buffer and
> removes the material-crossing problem at the BI surface entirely.
>
> Their reported condition number stays on the order of **10³** across the tests,
> including as `k₀a → 0`.

### 5.1 Their method does not cover our problem

**The paper contains no conductivity.** Permeability is `μ₀` throughout, the only
material variation is `ε(r)`, and every validation case is dielectric: a sphere
(`ε = 2.25`), a multilayered sphere, and dielectric elliptical cylinders on a
dielectric slab. A text search of the full paper returns **zero** occurrences of
conductivity, conductor, lossy or σ.

A THz phased-array feed is the opposite case — lossy metal with a skin depth,
sharp conductor edges, and conductor/dielectric interfaces where their
permittivity-gradient coupling term becomes a surface delta. **The 0‑form route
is not something we could adopt off the shelf**, and §6 explains why it is not
merely a matter of porting it.

---

## 6. The gauge and the element space are not independent choices

This is the central finding of the second edition, and it reframes the question.

### 6.1 The two representations

| | Whitney 1‑form | Whitney 0‑form |
|---|---|---|
| FEM name | **lowest-order Nédélec edge element** | P1 nodal (Lagrange) |
| basis | `w_ij = λ_i∇λ_j − λ_j∇λ_i` | hat function `λ_i` |
| degree of freedom | `∫_edge A·dl`, one per edge | the value at a node |
| continuity enforced | **tangential only** | **full** |
| what we use for `A` | ✅ this one | — |

In DEC the correspondence is exact: a discrete 1‑form's DOF is the integral along
a primal edge — the Nédélec edge DOF — and a 0‑form's DOF is the nodal value.

### 6.2 Abdrabou & Gomez did not use Nédélec elements — they avoided them

In their own words: *instead of modelling `A` on the boundary as a 1‑form, they
treat its Cartesian components `A_x, A_y, A_z` as independent 0‑forms*, which
they acknowledge departs from the standard DEC formulation. Their algorithm's
output is literally `discrete 0-forms a_x, a_y, a_z, φ_s`.

In finite-element terms that is **node-based vector elements** — three copies of
P1 Lagrange — the construction electromagnetic FEM normally forbids, for three
reasons: spurious modes in curl–curl eigenproblems, inability to represent the
singular field at a reentrant corner, and **being wrong at material interfaces**,
where full continuity is imposed on a quantity whose normal component physically
jumps.

### 6.3 Why they get away with it — the gauge decides

The third objection is the decisive one, and **which gauge you choose determines
whether it applies.**

```
Chew generalized gauge      ∇·(εA) = −χ ∂Φ/∂t
    ⟹  ∇·(εA) finite at an interface
    ⟹  n̂·(ε₁A₁) = n̂·(ε₂A₂)          Chew eq. (26)
    ⟹  n̂·A JUMPS wherever ε jumps
    ⟹  A must live in H(curl): Whitney 1-form / Nédélec

Abdrabou gauge              ∇·A = i k₀² ε(r) Φ        (divergence of A, not of εA)
    ⟹  ∇·A finite at an interface
    ⟹  n̂·A is CONTINUOUS
    ⟹  A may live in H¹: nodal 0-forms are legitimate

Coulomb gauge               ∇·A = 0
    ⟹  A solenoidal, n̂·A CONTINUOUS
    ⟹  nodal-natural; Ansari cites Biro & Preis for exactly this
```

**So the gauge and the basis cannot be chosen independently.** That is the real
content of "1‑form vs 0‑form", and it is not stated anywhere in the project's
existing documents.

| gauge | `n̂·A` at a material interface | element space it suits |
|---|---|---|
| Coulomb `∇·A = 0` | continuous | nodal |
| Abdrabou Lorenz `∇·A = ik₀²εΦ` | continuous | nodal |
| **Chew generalized `∇·(εA) = −χ∂Φ/∂t`** | **jumps ∝ 1/ε** | **edge / Whitney‑1** |

`A` is gauge-dependent, so its *continuity class* is gauge-dependent too. There
is no interface behaviour of `A` that is true independently of the gauge — only
a behaviour consistent with each.

### 6.3.1 This explains why Ansari needed a Lagrange multiplier — and why we would not

Ansari's diagnosis (§3.1) now reads as a **mismatch between gauge and element
space**, not as a defect of either:

- He wanted the **Coulomb** gauge, under which `n̂·A` must be **continuous**.
- He discretized `A` with **edge elements**, which permit `n̂·A` to **jump** — in
  his words, the normal component of an edge-element vector potential "is not
  necessarily continuous across the boundary interfaces", and linear edge basis
  functions "do not necessarily form a solenoidal vector potential."
- So nothing in the discretization enforced what the gauge required, and the
  potentials came out non-unique. The Lagrange multiplier is the machinery that
  repairs the mismatch.

**Chew's generalized gauge has no such mismatch.** It *wants* `n̂·A` to jump at an
ε discontinuity, which is precisely what edge elements naturally permit. The
gauge and our existing element space agree by construction rather than by
enforcement.

That is an argument for Chew over every Coulomb-family candidate that is
independent of the BEM argument in §4 — and it applies even if the hybrid never
gets built.

Their resulting system makes the consequence visible: each component gets a plain
scalar Helmholtz operator `∇²Ã_ν + k₀²ε Ã_ν`, and the **only** coupling between
components and Φ is through `∂_ν ε`, the gradient of permittivity. In a
homogeneous region the three components decouple into independent scalar
problems — which is exactly the regime nodal elements are built for.

### 6.4 What this means for us

**Our existing `A` space is already the one Chew's gauge requires.** The first
edition under-weighted this. Adopting Chew's generalized gauge means:

| | |
|---|---|
| `A` — Whitney 1‑form edge elements | **unchanged** |
| Φ — P2 nodal | → P1, to sit in Whitney‑0 and restore the de Rham structure |
| the gauge term in the `A` equation | new: needs the Hodge/Whitney-mapping machinery |
| tree–cotree | removed |

Switching instead to the 0‑form route would mean **replacing the `A`
discretization entirely**, adopting a formulation with **no demonstrated
conductor support** (§5.1), and accepting nodal elements at the sharp metal edges
of a THz feed — where objection (2) above still stands, gauge or no gauge.

### 6.5 The open question, stated precisely

Not *which gauge* — that is settled — but **whether `A` stays a Whitney 1‑form at
the boundary-integral surface.**

The fourteen operators arise *specifically* from the 1‑form trace at the
boundary. Abdrabou & Gomez's answer was to stop using it. Coupling Whitney 1‑form
`A` to a potential BEM means taking on exactly the problem they chose to avoid:
**that is the project's novelty claim and its principal technical risk, and the
proposal should say so** rather than treating the fourteen-operator trap as
something the gauge disposes of.

Chew's own §6 suggests a middle path worth evaluating before committing: his
equivalence principle needs **six scalar surface quantities**, not fourteen
operators — two tangential components of `H₁`, two of `A₁`, Φ and `n̂·A₁` — and
his surface currents are **divergence-conforming**, which is far closer to
Whitney-compatible than a Cartesian 0‑form decomposition.

---

## 7. Caveats, including two from Chew himself

**Chew prefers `E`–`H` in genuinely wave-dominated regimes.** §5:

> An advantage of (41) and (42) is that they are derivable from each other,
> whereas for the A‑Φ formulation, **only the Φ equation is derivable from the A
> equation.** Hence, when the `E` and the `H` fields are equally strong, and
> strongly coupled to each other, the `E`–`H` formulation is preferred, as it
> describes wave physics better.

This does not sink a THz solver — a phased-array feed is a multiscale problem
where local low-frequency breakdown is the binding constraint — but **the A‑Φ
advantage is multiscale robustness, not superior wave physics**, and the proposal
should claim the former rather than the latter.

**The low-frequency *inaccuracy* problem is distinct from breakdown.** At `ω = 0`
the electrostatic solution lives in (16) and the magnetostatic solution in (17),
so **both must be solved in tandem** to retrieve both accurately. Chew repeats
this for the integral equations: the vector-potential SIE captures
magnetoquasistatics, and a *scalar-potential* SIE is needed alongside it for
electroquasistatics. Budget for solving both, not one.

**Topology.** At `ω = 0` a PEC object of genus > 0 admits non-trivial
magnetostatic solutions with no excitation — superconducting loops. Our case 05
is a ring, so this is live for us, not hypothetical.

**What is still not verified.** Equations here were read from text extracted out
of PDF content streams; the prose is reliable, the typeset mathematics is not.
Every equation above should be checked against the typeset paper before it is
implemented.

**All cited papers have now been obtained and read.** Ansari and Rapetti were
supplied directly and are in `APhi_Papers/`; everything else was downloaded from
open-access sources. Nothing in this report now rests on a summary of a paper
rather than the paper.

---

## 8. Candidates, final ranking

| | verdict |
|---|---|
| **Chew generalized gauge** | **Recommended.** Only candidate where FEM, interface conditions and BEM are one formulation — and **our `A` space is already the one it requires** (§6.4). Costs: Hodge/Whitney-mapping machinery, a sparse approximate inverse, and moving Φ from P2 to P1. |
| Zhao & Fu dummy-variable Coulomb | **Withdrawn.** Symmetric, direct-solvable and a good architectural fit — but Coulomb-gauged, so it cannot couple to a potential BEM. Implementing it means implementing a gauge twice. |
| Ansari explicit Lagrange-multiplier Coulomb | Same Coulomb objection, plus a specific one: it exists to repair a **gauge/element mismatch we would not have** (§6.3.1). Retain as the roadmap's benchmark comparison — it is the documented failure mode a new gauge must not reproduce. |
| Improve tree–cotree (root choice, MOR) | Cheap, treats the symptom only. Lee & Jin themselves report the condition number still grows at low frequency. |
| Do nothing | Viable *only* while every port sits on an `n × A = 0` outer boundary. That excludes interior and mixed conductor/dielectric ports **at any frequency** (§1.2), and ends entirely when the boundary treatment changes (§1). |

**A design decision the roadmap does not currently mention:** Chew's construction
places Φ in the Whitney‑0 space, i.e. **P1 nodal**. **Ours is P2.** The
compatibility argument rests on the Whitney/de Rham structure, so adopting the
method means either moving Φ to P1 or redoing the analysis for P2.

That is the *whole* discretization change. `A` stays exactly as it is — which,
given that `n̂·A` must be free to jump under this gauge (§6.3), is not a
coincidence but a consequence of having picked edge elements correctly in the
first place.

---

## 9. What to measure before committing

1. **Is Φ mesh-independent?** Refine twice; Φ must converge pointwise.
   Tree–cotree fails this by construction. **Case 07's permutation harness is
   directly reusable** — under a real gauge, its *potential-moves* assertions
   should invert while every field assertion holds unchanged.
2. **Condition number against frequency**, to the static limit, against the
   tree–cotree baseline.
3. **The mixed conductor/dielectric port**, with Chew's embedded interface
   conditions — does the port behave with no special-casing?
4. **Both quasi-static limits in tandem** (§7): electroquasistatic from the Φ
   equation, magnetoquasistatic from the `A` equation, both accurate at once.
5. **Factorization cost and memory** at case 05/06 scale.
6. **A genus > 0 case** — case 05's ring — for the null magnetostatic solution.
7. **Solver-independence**, following Ansari (§3.1): solve the same gauged system
   with a direct and an iterative solver and require the *potentials* to agree,
   not just the fields. It is a cheap check and it is the one that exposed the
   problem for him.

**And one measurement to take first, before any gauge work begins.** Build a
minimal case with a port face **inside** the domain crossing conductor and
dielectric, and run case 07's tree-permutation harness on it. Two outcomes, both
worth having:

- **Φ and the terminal quantities move with the tree** — the predicted failure is
  now *measured* in this project rather than inherited from the literature, and
  there is a quantitative baseline the new gauge must fix.
- **They do not move** — then the premise behind the whole gauge programme needs
  re-examining before anything is rewritten.

This is cheap: the harness exists, the solver exists, and the case is smaller
than case 06. It should run at low frequency, where §1.2 argues the failure
already appears.

---

## 10. Corrections to the project's own bibliography

Four entries in `Gemini/APhi_FEM_BEM_Comprehensive_Report_Updated.tex` are wrong
and have propagated into the proposal's novelty claims:

| claim | correction |
|---|---|
| "Bogaert, Cools & Andriulli, arXiv:1705.00265" | That arXiv ID is **Li, Fu & Shanker**, *Potential Integral Equations in Electromagnetics*. |
| Abdrabou & Gomez use "structured dual grids" | They use **unstructured meshes**. The real distinction is **0‑form Cartesian components vs Whitney 1‑forms** — the novelty claim must be re-grounded on that axis. |
| The generalized Lorenz gauge is what reduces the boundary system to `V` and `K` | The **Cartesian scalar reformulation** performs the reduction; the gauge enables it. Two steps. |
| Chew 2014: "Generalized Lorenz Gauge … Anisotropic Media", pp. 69–77 | Actual title *Vector Potential Electromagnetics with Generalized **Gauge** for **Inhomogeneous** Media: Formulation*, pp. 69–**84**. |

Also: `AphiFreqDomainTreeCotree_SeungJFlee2003.pdf` is mis-filed — it is Lee,
Lee & Lee, *Hierarchical Vector Finite Elements for Analyzing Waveguiding
Structures*, IEEE T‑MTT **51**(8) 2003, a 2‑D waveguide eigenanalysis paper that
says nothing about Φ's non-uniqueness.

---

## 11. The Darmstadt / ABB / Siemens line — the literature does settle the spike's finding 2

Added after `GENERALIZED_LORENZ_GAUGE.md` §11 reported that `ε_eff = ε − jσ/ω`
wrecks conditioning below ~1 GHz and concluded **"this is now the main obstacle
and nothing should be implemented before it is settled."** That conclusion was
reached without checking the roadmap's own reference list, and it is wrong in its
premise: the phenomenon is named, published, and has remedies. It is a *known*
obstacle, not a new one.

Four papers, one continuous research line (TU Darmstadt + ABB + Siemens Digital
Industries), none of which appeared in §8's ranking.

### 11.1 Balian et al. (2023) — our finding 2, by name, with a remedy

> *"When simulating resistive-capacitive circuits or electro-quasistatic problems
> **where conductors and insulators coexist**, one observes that large time steps
> or **low frequencies lead to numerical instabilities, which are related to the
> condition number of the system matrix**. Here, we propose several stable
> formulations **by scaling the equation systems**."*

That is the spike's finding 2 restated as a paper's abstract, and the mechanism
they name is the same one: the `σ/ωε` ratio between conducting and insulating
regions.

**The remedy is a scaling applied analytically, before assembly.** Multiply
equation block `k` by `a_k` and substitute the unknown `ψ_k = b_k⁻¹ φ_k`:

| | |
|---|---|
| (i) symmetric | `a₂ = b₂ = ω^(−1/2)`, `a₁ = b₁ = 1` |
| (ii) non-symmetric | `a₂ = ω⁻¹`, rest 1 |
| (iii) **symmetric with material** | `a₁ = b₁ = (σ₁ + jωε₁)^(−1/2)`, `a₂ = b₂ = (ε₂ jω)^(−1/2)` |
| (iv) non-symmetric with material | `a₁ = (σ₁ + jωε₁)⁻¹`, `a₂ = (ε₂ jω)⁻¹`, `b = 1` |

Three things in that table matter to us:

1. **`a_k = b_k` is a congruence** — it scales the equation *and* the unknown, so
   symmetry survives. The implementation spec's §11 framed our choice as
   "symmetric or well-conditioned, not both". Variants (i) and (iii) are both.
   §11.4 below is the qualification on how far that carries to our matrices.
2. **Variant (iii) scales by `(σ + jωε)^(−1/2)`** — exactly our `β`. The material
   combination the spike found fatal is the one they scale *by*.
3. **It must be done before assembly**, symbolically: *"the products of powers of
   ω must be determined before matrix assembly to avoid numerical errors or
   division by zero."* A post-assembly Jacobi preconditioner is variant (iv)
   applied too late (their §III, citing [17, §4.1]).

Measured: the unscaled RC system breaks down below 1e10 Hz; (i) and (iii) stay
flat to 0 Hz, with (iii) better than (i) by `1/(2RC)` ≈ 5e11. Two caveats the
authors state themselves — **(i) destabilizes above 1e10 Hz**, where the problem
becomes essentially capacitive, and under (i)/(iii) **`φ` cannot be recovered at
exactly ω = 0** because the unknown is the scaled one, *"a natural consequence of
the fact that it is not well-defined from the start."* Variants (ii)/(iv) keep
the original unknowns and lose symmetry.

This covers the **Φ equation only.**

### 11.2 Herles et al. (2025) — the A equation, and it keeps tree–cotree

[arXiv:2502.13588](https://arxiv.org/abs/2502.13588), IEEE Trans. Magn. The
companion that stabilizes the other half; they say so outright — *"Effective
modifications are proposed in [14] to improve the condition number. The focus of
our paper is the stabilization of (11)"*, (11) being the curl–curl system.

Their scheme is Ostrowski & Hiptmair's two-step (§11.3): solve an
electro-quasistatic problem for `φ`, then use it as the source for `A`, with
**the EQS problem itself serving as the gauge condition**. Taking `div` of the
`A` equation and dividing by `iω` gives the implicit constraint

```
    div(κ A) = 0,          κ = σ + iωε                                  (17)
```

a **generalized Coulomb gauge with complex conductivity**. It degenerates in the
insulator as `ω → 0`, so they split it by region:

```
    α div((σ + iωε) A) = 0      in Ω_C   (conductor)
    β div(ε A)         = 0      in Ω_A   (insulator)                    (18)
```

**Each piece is frequency-independent in its own region.** This is the structural
option the spike never considered: *a different gauge condition inside conductors
than outside*. Chew carries one `χ` everywhere, and the spec's §11 listed "a
different `χ` inside conductors" as a speculative line item — here it is, built
and measured. Their weights are `α = 1 + ω` and
`β = (1 + ω) max σ + σ_art max ε` with `σ_art = 1e−6`: **the scaling carries the
peak material magnitude**, the same idea as Balian's variant (iii).

The discretization is the part worth copying. Build the weighted divergence
matrix `(S_?)ᵢⱼ = ∫_? div(w_j) v_i dV`, take the tree–cotree split of the edge
DOFs, then **replace the redundant tree rows of the curl–curl system with the
gauge-constraint rows**:

```
    ⎡ W^(RR)   W^(RT) ⎤ ⎡ a^(R) ⎤   ⎡ j^(R)(u) ⎤
    ⎣ S^(R)    S^(T)  ⎦ ⎣ a^(T) ⎦ = ⎣    0     ⎦                        (26)
```

They prove `W^(RR)` has full rank even at `ω = 0`, and the second row of the
cotree system is automatically satisfied, so nothing is lost by the swap. The
symmetric alternative is the Lagrange-multiplier saddle point (21), which they
reject on size.

**Tree–cotree is not discarded — it is what tells them which rows are
redundant.** That is a very different and much cheaper change than replacing the
gauge: we already have a spanning tree.

Measured, on a mesh of **three conducting bars in a dielectric box** with `φ = 0`
and `φ = 1` on opposite faces:

| | cond at f = 0 |
|---|---|
| original | **singular** |
| stabilized | 4.87e5 |

and on a copper planar coil (`σ = 6e7`), 1.9e9 against singular. The discrete
gauge residual `‖S·a‖₂` stays below 1e−11 across the whole sweep and is 1.52e−12
at DC, against ~1e−2 unstabilized. Fig. 8 is the one to look at: the unstabilized
`‖B‖` and `‖E‖` are **plotted logarithmically because they reach 1e70 and 1e75** —
garbage of exactly the kind a mixed-material port produces.

Two further details bear on our problem directly:

- **Their DOF partition assigns interface nodes to the conductor.** Index sets
  (12)–(13): `I_v^(A)` is the nodes whose support does *not* intersect the
  conductor, and `I_v^(C)` is *everything else*, so a node on the
  conductor/dielectric interface is a conductor node. They do not try to make a
  per-element material factor behave — they partition, then scale blockwise.
- **Their §V-D is close to our case 05.** A thin dielectric slit inside a
  conducting loop (`ε_r = 7.2e15`, acting as a lumped capacitor), with the current
  closing through the gap via `iωD_m`. A dielectric gap interior to a conductor,
  carrying current — structurally our internal cut — handled by the generalized
  gauge.

### 11.3 Ostrowski & Hiptmair (2020/2021) — the origin

ETH SAM Research Report 2020-43, published as SIAM J. Sci. Comput. **43**(4)
B1008–B1028 (2021). Introduces EQS-as-gauge-condition and the two-step procedure,
with frequency-stable weak forms for both steps:

> *"the electro-quasistatic fields can be corrected for magnetic/inductive
> phenomena at any frequency in a second step. The combined field from both steps
> is a solution of the full Maxwell's equations… Electro-quasistatics serves as a
> gauge condition in this semi-decoupled procedure."*

Structurally this is Chew's move — an independent `Φ` equation, then an `A`
equation taking `Φ` as a source — reached from the low-frequency side with a
different gauge condition and no `χ`.

### 11.4 What transfers to Chew's formulation, and what does not

Being careful here, because these papers stabilize the **generalized Coulomb**
gauge, not Chew's generalized Lorenz, and §8's ranking must not be rewritten on a
false equivalence.

**Transfers directly:**

| | |
|---|---|
| the diagnosis | `σ/ωε` contrast as the mechanism is confirmed and published — not a bug in our spike. Herles §III-B states it in a line: *"The ratio σ/ωε ≫ 1 in Ω_C increases this issue even further."* |
| the technique | analytic pre-assembly scaling by material-dependent powers of `ω` is formulation-independent, and applies to the §7 block system as written. |
| `a_k = b_k` | a two-sided congruence keeps symmetry. The spike's dichotomy was too quick. |
| per-region gauge | a different gauge condition in conductor and insulator is a real, measured option, not speculation. |
| partition, don't diagonalize | assign interface DOFs to the conductor set and scale blockwise. |

**Does not transfer, and this is the honest limit:**

- Our symmetry problem is **not** a block scaling.
  `K_NE[m,n] = −∫(ε/χ)∇λ_m·ω_n` against `K_EN[m,n] = ∫ε ω_m·∇λ_n`: the factor
  `ε/χ` sits *inside the element integral*. If it is constant per region, then for
  a node strictly interior to one region every element in its support shares the
  factor, and that row of `K_NE` is exactly `−c` times the matching row of
  `K_ENᵀ`. **For a node on a material interface it is not**, under any diagonal
  scaling. So a congruence repairs symmetry everywhere except at interface
  nodes — which is precisely where a mixed port lives. Balian and Herles sidestep
  this by assigning those nodes to one block, not by exactness.
- **Neither paper addresses a port on a mixed-material face.** Herles's dielectric
  slit is an interior gap with no potential boundary condition on it. §11.2's
  evidence is suggestive for case 05, not an answer for the coax via.
- Their target is DC-to-MHz industrial devices. **Nothing here is validated at
  THz**, which is the stated end goal.

### 11.5 What this changes

1. **`GENERALIZED_LORENZ_GAUGE.md` §11 is corrected** — the `σ` question is a known
   problem with published remedies, not a blocker on implementation.
2. **A cheaper route exists and belongs in §8's ranking**: keep our spanning tree,
   add the weighted divergence matrix `S`, swap the tree rows per Herles (26).
   That is incremental on the solver we have, where Chew is a rewrite. It gives up
   symmetry as (26) is written, and it is a Coulomb-family gauge, so it inherits
   the `n̂·A` continuity problem of §6 — the reason we went looking in the first
   place. **It is not a substitute for the mixed-port fix; it is a substitute for
   the low-frequency-conditioning half of the argument**, and those two
   motivations should stop being bundled.
3. **Reference code exists**: Herles, *Low Frequency Stable Full Maxwell*,
   [doi:10.5281/zenodo.14810885](https://doi.org/10.5281/zenodo.14810885) —
   GeoPDEs / Octave. Worth reading before writing our own `S`.

### 11.6 Still not obtained

| | why it matters |
|---|---|
| **Eller, Reitzinger, Schöps & Zaglmayr (2017)**, SIAM J. Sci. Comput. **39**(4) B703–B731, [doi:10.1137/16M1077817](https://doi.org/10.1137/16M1077817) | *"monolithic, symmetric, low-frequency stable, broadband… no auxiliary variables… stable even if the frequency equals zero."* Symmetric **and** stable **and** monolithic answers all three of the spike's complaints at once. Paywalled, no preprint — but **not worth buying**, because §11.8 obtained its method from two open sources. |
| **Jochum, Farle & Dyczij-Edlinger (2015)**, IEEE Trans. Magn. **51**(3) 7402304, and the 2016 SCEE companion *A symmetric and low-frequency stable potential formulation* | Previously flagged, still unread. |
| **Demerdash & Wang (1990)** | Coulomb-gauge breakdown at permeability contrast. Still unread. |
| **Zhu & Jiao (2010)**, IEEE Trans. Adv. Packag. **33**(4) 1043–1050 | *"theoretically rigorous full-wave FEM solution of Maxwell's equations from dc to high frequencies"* — the reference both Darmstadt papers cite for the breakdown itself. |
| **Manges & Cendes (1995)**, IEEE Trans. Magn. **31**(3) 1342–1347 | *A generalized tree-cotree gauge for magnetic field computation* — origin of the term Herles et al. use. |
| **Hiptmair, Kramer & Ostrowski (2008)**, IEEE Trans. Magn. **44**(6) 682–685 | *A robust Maxwell formulation for all frequencies.* |
| **Chew's [28] White & Koning (2002)** and **[29] Dai, Chew & Jiang (2013)** | Sources for the `Ḡ⁻¹` elimination. Still unread. |

### 11.7 New in `APhi_Papers/`

| file | |
|---|---|
| `Balian2023_LowFreqStab_ConductorsInsulators.pdf` | Balian, Merkel, Ostrowski, De Gersem & Schöps, [arXiv:2302.00313](https://arxiv.org/abs/2302.00313), IEEE Trans. Dielectr. Electr. Insul. **30**(6), 2023. §11.1. **Downloaded under `DyczijEdlinger2023_…` and renamed** — the roadmap attributes it to Dyczij-Edlinger, who is not an author. |
| `TwoStep_GeneralizedTreeCotree_LowFreqStability_2025.pdf` | Herles, Mally, Ostrowski, Schöps & Merkel, [arXiv:2502.13588](https://arxiv.org/abs/2502.13588). §11.2. |
| `ETH2020_FreqStableFullMaxwell_EQSGauge.pdf` | Ostrowski & Hiptmair, ETH SAM 2020-43. §11.3. |
| `TwoStep_TimeDomain_Stabilized_2025.pdf` | [arXiv:2507.18235](https://arxiv.org/abs/2507.18235), same group, time-domain extension. Not relevant to a frequency-domain solver; filed for completeness. |
| `Clemens2022_DarwinTypeQuasistatic.pdf` | [arXiv:2204.06286](https://arxiv.org/abs/2204.06286), Darwin-type quasistatic formulations. Cited by both; unread. |
| `ShinFan2013_EigenvalueEngineering_OE21.pdf` | Source of Chew's `α > 0` eigenvalue claim. Extracted, unread. |

---

### 11.8 Eller's method, obtained without Eller — and a competing design for our exact application

Both sources below are open access and were obtained after §11.6 listed Eller et
al. as the paper we most wanted. **Eller no longer needs to be bought.**

| | |
|---|---|
| `Stysch2022_BroadbandFEM_Impedance_ParasiticExtraction.pdf` | Stysch, Klaedtke & De Gersem, *Electrical Engineering* **104**(2) 855–867, 2022, [doi:10.1007/s00202-021-01348-9](https://doi.org/10.1007/s00202-021-01348-9). Also [arXiv:2009.08232](https://arxiv.org/abs/2009.08232) — it was free all along. |
| `Stysch2022_PhDThesis_StableBroadbandParasiticExtraction.pdf` | Stysch, PhD thesis, TU Darmstadt, 2022, handle `tuda/8889`. **Chapter 6, 24 pages, is a full treatment of Eller's scheme**, including §6.6.1 tree–cotree splitting, §6.6.2 gradient-space splitting, and §6.3.2 *Stabilizing the Lorenz-gauged system matrix* — which Eller's own paper, being E-field only, does not cover. |

Robert Bosch GmbH + TU Darmstadt (De Gersem is also a co-author of Balian et
al., §11.1). The application is **parasitic extraction for EMC** — ours.

#### Eller's mechanism, as the article states it

Split the trial and test space three ways, `H(curl, Ω) = V ⊕ W ⊕ U`, with

```
    V :  ∫_Ω |curl v|² dV  ≠  0                                      (36a)
    W :  curl w = 0   and   ∫_Ωc |w|² dV  ≠  0                       (36b)
    U := { u ∈ H(curl,Ω) : curl u = 0  ∧  u = 0 in Ω_c }             (36c)
```

so `V` is the part with real curl, `W` the gradient fields that live in the
conductor, `U` the gradient fields that vanish there. Then scale the three parts
by **different powers of ω**:

```
    E  =  jω E_V  +  (jω)^(1/2) E_W  +  E_U                           (37)
```

and test the equation separately with each subspace, which scales each of the
three resulting equations independently. **That fractional `ω^(1/2)` is what
Balian et al. meant** by *"a scaling by fractional powers of ω was similarly
applied in [12]"* — their [12] is Eller et al. Same family as §11.1's variant
(i), applied to a three-way space split rather than a two-way material one.

**At lowest order the `V` basis is found by a tree–cotree split** (their [19]).
That is now the third independent paper in this section that keeps the spanning
tree and changes only what is done with it.

#### The part that bears on case 06 and 07 directly

§2 of the article is the cleanest published statement of what our gauge tests
measured. They start from the problem by name — *"To reconcile the
path-dependent voltage concept of electromagnetic field theory with the
path-independent voltage concept of electrical circuits poses a challenge"* —
and resolve it by subtracting the return path's partial inductance:

```
    V_c  :=  V − jω I L_r  =  ∫_c E·dl + jω ∫_c A·dl  =  Φ(T_b) − Φ(T_a)
```

their (5)–(7). So the terminal Φ-difference **is** a path-independent quantity,
and it is *not* the full terminal voltage: it is the voltage minus the return
path's partial inductance. And then, explicitly:

> *"The values of these partial inductances depend on the gauge condition chosen
> for the magnetic vector potential A and electric scalar potential Φ."*
> *"As it incorporates a partial inductance, V_c formally depends on the gauge
> condition of potentials Φ and A."*

**This is our ωL/R mechanism, published.** The gauge-sensitive part of a terminal
Φ-difference is exactly the split of the loop inductance into the conductor's
partial inductance and the return path's — `L_loop = L_c + L_r` is gauge-invariant
while the two terms separately are not (their (4)). `07_GaugeInvariance` measured
the same thing from the other end.

The consequence for how we present our results: a tree-dependent Φ is **not
simply an error to be driven to zero**. For parasitic extraction the partial
inductance is the quantity you want, and the gauge is what *defines* how the loop
inductance is attributed between conductor and return path. It is a modelling
convention that must be fixed, not a defect. What is a defect is leaving it
*implicitly* fixed by a spanning-tree traversal order, which is what we do today.

#### How they make Φ unique — the answer to the question that started this file

They impose the gauge as **an explicit PDE for Φ**, and solve it. Lorenz gauge
(their (12), `div(ν_r A) + jωΦ/c² = 0`) combined with `E = −grad Φ − jωA`
eliminates `A` and leaves a scalar boundary-value problem:

```
    −div(ν_r grad Φ) − (ω²/c²) Φ  =  div(ν_r E)                       (13)
```

with, after the inductive compensation term is added (their (27)):

```
    Φ_c = const.          on Γ_el
    n̂ · ν_r grad Φ_c = 0  on Γ_mag
```

**`Φ = const.` on the electrode, with the constant not prescribed, is a floating
potential condition** — the "single potential port, no ground defined" behaviour
you described from Maxwell, as a boundary condition rather than a special port
type.

So Φ is unique because it solves its own well-posed BVP. No tree, no gauge
ambiguity in Φ at all. The gauge condition has become the equation that
determines it. Structurally this is Chew's decoupled Φ-equation and Ostrowski &
Hiptmair's EQS-as-gauge (§11.3) reached a third way.

Their full procedure is three sequential solves — `g`, then `E`, then `Φ_c`
(their Fig. 2) — and the MQS limit is reached by dropping the `ω²` term, which
they note *"is equivalent to choosing the Coulomb gauge condition for the
calculation of Φ"*, necessary at high frequency for consistency between the two
PDEs.

#### Why this is a competing design and not just a reference

They considered our formulation and rejected it. Quoting their §3.1 on the two
options for computing Φ:

> *"It is more advantageous to use the latter 'E approach' for several reasons:
> The two fields E and Φ can be calculated in sequence… thereby avoiding a
> computationally more expensive coupled boundary value problem (BVP), which
> occurs in the 'A‑Φ approach'. Furthermore, the E approach allows for an easy
> treatment of conductors modeled as perfect electric conductors… Finally, the use
> of the E-field formulation allows for an elegant stabilization of the
> low-frequency instability."*

Three reasons, from the same institute, for the same application. This belongs in
§8's ranking as a genuine alternative to A‑Φ, not a footnote. Against it, for our
purposes: `E` in `H(curl)` gives up the potentials we need for a potential-based
BI coupling (§4), and their own §4.1.1 shows the volumetric discretization
saturating the resistance above ~10 MHz, which is a hard limit for a THz target.

#### Still unread here

Chapter 6 of the thesis. Neither our extractor nor TUprints' own text layer
recovers it — TUprints' stops mid-chapter 5, and the thesis PDF yields only 14 k
characters from 148 streams because of its font encoding. The PDF renders
normally and can be read page by page; chapter 6 spans pp. 63–86.

---

### 11.9 Where this leaves the decision

Three papers now keep tree–cotree and change what is done with it (Herles,
Eller/Stysch, and the Manges–Cendes line they both cite), and two independent
groups obtain a unique Φ by **making the gauge condition an explicit equation for
Φ** rather than by fixing degrees of freedom on a tree (Chew via `χ`, Stysch via
the Lorenz PDE, Ostrowski & Hiptmair via EQS). That convergence is the strongest
signal in this document, and it points at the same place from three directions:
**stop setting tree DOFs arbitrarily and give Φ its own equation.**

What that does *not* settle is the mixed-material port, which is still the
project's actual blocker. Stysch's `Φ_c = const.` on an electrode is the closest
published thing to the floating-potential port, and it is a boundary condition on
a scalar BVP — which is only available once Φ has its own equation. The two
problems may therefore have one solution, but that is a conjecture and is not
demonstrated by anything read so far.

---

### 11.10 Chapter 6 of the thesis, read in full — and it contains the floating port

§11.8 left chapter 6 unread because the thesis PDF defeated our extractor. It is
a CID-font document: text is hex-coded glyph ids, not characters, and `ex2.py`
only understands literal strings, so it recovered 14 k characters from a
155-page thesis. `tools/pdftext_cid.py` decodes through each font's `/ToUnicode`
CMap (306 of them here) and splits output by page: 267 k characters, and
chapter 6 is pp. 63–86, which are PDF pages 75–98.

Two sections of it are **directly implementable in our solver today, whichever
gauge we end up with.**

#### 11.10.1 What Eller's method actually is

Not a new formulation — *a synthesis of two existing ones* (thesis p. 64):

| from | what it contributes |
|---|---|
| **Hiptmair, Kramer & Ostrowski (2008)** | the splitting of the scalar space `H¹`. Their own "generating systems" approach enforces Gauss's law in the non-conducting domain only and the Coulomb gauge in the conductor — but yields matrices that are **singular and non-symmetric**. |
| **Jochum, Farle & Dyczij-Edlinger (2015/16)** | the view that **tree–cotree *is* a Helmholtz-type decomposition of `H(curl)`**. Their formulation is symmetric and regular even in the lossy case, by treating conducting and non-conducting regions with separate equations plus interface conditions on `∂Ω_c` — at the cost of extra DoFs in conductors and *"significantly less sparse matrices, increasing the computation times considerably."* |

Eller unites the two splittings. Note what this does to §11.6's wish-list: the
two Jochum papers we still have not read are *inside* this, and **in the PEC case
Eller's equations and Jochum's are equivalent** (p. 64).

The decomposition is two-stage. First a Helmholtz split off the first level of
the de Rham complex, then a split of the gradient part by whether it survives in
the conductor:

```
    H(curl,Ω)′  =  V ⊕ Y                                          (6.11)
        V :  ⟨curl v, curl v⟩ ≠ 0   (functions with real curl)
        Y :  curl w = 0             (gradient fields)

    Y  =  W ⊕ U                                                   (6.12)
        W :  ⟨w, w⟩_Ωc ≠ 0          (gradients that live in the conductor)
        U :  u = 0 in Ω_c           (gradients that vanish there)

    H(curl,Ω)′  =  V ⊕ W ⊕ U                                      (6.13)
```

Then scale the three parts by different powers of the wave number `k = ω/c` and
test each with its own subspace:

```
    E  =  jμ₀c ( k E_V  +  √k E_W  +  E_U )                        (6.14)
```

`j` is folded into the scaling deliberately, so that the discretized system
becomes **real** in the PEC case.

**The sentence that matters most to us** (p. 67):

> *"The magnetic flux density can be recovered with the simple expression
> `B = −μ₀ curl E_V`, illustrating that **`E_V` is essentially a MVP**, while
> `E_W` and `E_U` can be understood as **gradient fields produced by ESPs**."*

So Eller's three-way split of `E` *is an A‑Φ decomposition in disguise* — `E_V`
plays `A`, and `E_W`/`E_U` are `∇Φ` separated by whether the potential survives
inside the conductor. The apparent fork between "the E-field route" and "our A‑Φ
route" is much narrower than §11.8 implied: the space structure is the same one,
and the scaling (6.14) is a statement about the **relative weighting of A against
∇Φ**, which is a thing we can apply directly.

The resulting weak system (6.15a–c) is **symmetric**. It is LF-stable in the
sense that the matrix has full rank at every frequency — but the right-hand
sides of (6.15b,c) still diverge as `k → 0`, so a second stabilization solves
the static fields `F_U` (electrostatic) and `F_W` (static current) first, in the
static limit, and writes `E` as those plus non-static corrections (6.17). That
gives (6.22), with the same operator and a bounded RHS. A useful special case:
**if both terminals lie on the same conductor, `F_U = 0`** and the impedance can
be evaluated exactly at DC.

#### 11.10.2 Our `ε_eff` instability, named, in their equation (6.1)

The spike's finding 2 appears here as a one-line diagnosis. Enforcing Gauss's law
explicitly over all of `Ω` in a potential formulation gives

```
    div[ (ε + σ/jω)(grad φ + jωA) ]  =  (1/jω) div J_s             (6.1)
```

and the thesis says of it: *"This equation contains the unstable term
`(jω)⁻¹ div σ grad φ` on its left-hand side."* Plus, for a source current that is
not divergence-free, the RHS diverges too.

`ε + σ/jω` is our `ε_eff`. **The unstable object is specifically
`(jω)⁻¹ div(σ ∇φ)`** — not the material contrast in the abstract, but that one
term. That is a sharper statement than `GENERALIZED_LORENZ_GAUGE.md` §11 reached
by measurement, and it says where to look: the `σ`-weighted gradient block of the
`Φ` equation, at the conductor boundary.

#### 11.10.3 Measured

Conductor segment, 10 mm × 1 mm², 680 elements, so the condition number and the
*rank deficiency* are both computed exactly:

| | unstable (6.8) | stabilized (6.22) |
|---|---|---|
| matrix dimension | 4518 | 4519 — **one extra DoF** |
| condition number | rises sharply below 10 kHz, then plateaus on machine precision | approximately **constant at all frequencies** |
| rank deficiency at 10 kHz | 0 | 0 |
| rank deficiency at 1 kHz | **630** | 0 |

On a capacitor–coil model the unstabilized impedance is *"strongly scattered"*
below ~100 kHz; at 5 kHz the unstabilized `|E|` runs to 1e8–1e15 V/m with an
implausible distribution (their Fig. 6.6).

Two results worth more than the stabilization itself:

- **Fig. 6.7: at high frequency the stabilized and conventional systems give
  identical values**, so *"no switching of systems is necessary."* One
  formulation covers DC to resonance. That is the property our roadmap wants and
  that Balian's variant (i) explicitly does *not* have (§11.1: it destabilizes
  above 1e10 Hz).
- **Cost, their Table 6.1**, on a 520 k-DoF wire model:

| system | N (stab/unstab) | non-zeros | time | ratio |
|---|---|---|---|---|
| full wave | 520 k / 520 k | 27.3 M / 22.7 M | 86.5 s / 67.9 s | **1.27×** |
| MQS | 435 k / 520 k | 13.7 M / 13.9 M | 36.1 s / 43.0 s | **0.84×** |
| MS (PEC) | 336 k / 520 k | 8.71 M / 9.54 M | 38.2 s / 44.7 s | **0.86×** |

Stabilization costs 27 % in the full-wave system and **saves** 14–16 % in the
quasistatic ones. The +20 % non-zeros in the full-wave row are attributed
specifically to constructing `U_h` — §11.10.5 below, which is the part that
densifies the matrix.

#### 11.10.4 The tree-construction rules Eller omitted

§6.6.1, citing Klis 2015, and flagged in the thesis as *"crucial details
especially with regards to the tree creation that were omitted in [Ell+17;
Ell17]"* — i.e. **the thesis is strictly better than the paper for implementing
this**, which settles §11.6's question about buying Eller.

The rules, for a boundary split into magnetic, electric and absorbing parts:

| | |
|---|---|
| `Γ_mag` | needs no special treatment |
| `Γ_el` | edge DoFs on it are set to zero |
| `Γ_a` (absorbing) | its basis functions **must** belong to `V_h`, by (6.11b) |
| consequence | **edges of `Γ_el` and `Γ_a` must be assigned to the tree** |
| compatibility | for each *disconnected* part `Γ_i` of `Γ_el ∪ Γ_a`, **contract all of its mesh nodes to a single vertex `v_i` in the graph** before spanning it |

In the PEC case each conductor's surface is one such disconnected boundary. Then
the ordering rule, which is the kind of thing that only shows up when someone has
actually built it:

> *"numerical artifacts can arise if in the non-conducting domain `Ω₀` the tree
> would not be equivalent to a PEC case tree that contains only one vertex `v_i`
> per disconnected boundary. The problem is avoided by **first constructing a
> tree in `Ω₀` respecting the boundaries of (6.50) and subsequently adding the
> edges of trees constructed in the individual conductors**."*

So: contract each boundary component to one vertex, span the **insulator first**,
then span inside each conductor. Our tree is built by a traversal over the whole
mesh with no such structure, and `07_GaugeInvariance` measures the consequence.
**This is implementable now and is independent of the gauge decision.**

#### 11.10.5 The floating-potential port, as a congruence

§6.6.2. This is the answer to the port question that has been open since you
described Maxwell's double-potential floating ports, and it arrives from an
unexpected direction — it is not a port model at all, it is how they *construct
the space* `U_h` of gradient fields that vanish in the conductors.

A potential `ψ_u` whose gradient vanishes in `Ω_c` must be **constant on each
disconnected conductor** `Ω_c,i`, and on each disconnected part of `Γ_el` and
`Γ_a`. Call the number of such entities `n_γ`. Then:

> *"The first constant can be chosen zero, and all `n_γ − 1` others **must become
> DoFs in the linear system**."*

One scalar unknown per floating entity, with the first grounded. That *is* a
floating-potential port. The implementation (6.51):

1. discretize with the ordinary nodal gradient basis, giving `Â`, `b̂`;
2. zero the DoFs of the first entity, and zero all **higher-order** DoFs on the
   remaining `n_γ − 1` entities, since first order already suffices to represent
   a constant;
3. aggregate the first-order DoFs of each remaining entity onto one DoF with a
   0/1 matrix `P` (their example (6.52) is a row of ones over one conductor's
   nodes), and form

```
    A = P Â Pᵀ ,        b = P b̂ ,        recover  x̂ = Pᵀ x        (6.51)
```

**`A = P Â Pᵀ` is a congruence, so symmetry is preserved exactly** — and `P` is
*not diagonal*. That is the general form of the device §11.4 was looking for
when it concluded a diagonal congruence fails at material-interface nodes: the
right object is an **aggregation** matrix, not a scaling. It costs sparsity (the
+20 % non-zeros of §11.10.3) because each aggregated row is dense over that
conductor's nodes.

Three things follow:

- **A floating port needs no new formulation and no new gauge.** It is one
  aggregated nodal DoF per floating conductor, applied as a congruence. This can
  be built on the solver we have.
- It answers your question *"do we work on floating double potential ports first,
  then the gauge, or both together?"* — **they are separable**, and the port is
  the smaller, lower-risk piece.
- It does **not** by itself fix the mixed-material port. Aggregating a face's
  nodes to one DoF presumes that face is an equipotential, which is what a
  conductor electrode is and a mixed conductor/dielectric face is not. For the
  coax via the dielectric annulus must stay free (`COAX_PORT_ANALYSIS.md` line
  66 already says so). So this gives us the floating *conductor* port, and leaves
  the mixed port where it was.

#### 11.10.6 And the Φ equation is the easy half

§6.3 states plainly that the BVP for `Φ` *"does not possess a general LF
breakdown"*: if the boundary is not purely magnetic the matrix has full rank at
all frequencies, and if it is, only the spatially constant component of `Φ` is
undetermined — fixed by pinning `Φ` to zero at one arbitrary point (their (6.29)).

The one exception is the **Lorenz-gauged** `Φ` equation with a purely magnetic
boundary, where the `−k²⟨Φ,ψ⟩` mass term fixes the constant at finite frequency
but vanishes as `k → 0`, leaving a rank deficiency of exactly 1 — harmless for an
iterative solver, fatal for a direct one. Their fix (§6.3.2) is the same move
again, one dimension smaller:

```
    H¹(Ω)  =  H¹_f(Ω) ⊕ ℂ                                         (6.30)
    Φ_c    =  jμ₀c ( φ_v  +  φ₀ / k )                             (6.31)
```

— split off the constant mode as **its own unknown, scaled by `1/k`**, and test
with all of `H¹_f` plus the single constant function `1`. The result (6.32) is
*"still symmetric but less sparse"*. Same pattern as §11.10.5 and as Balian's
variant (i): isolate the mode that degenerates, give it its own row, and scale it
by the power of `ω` that keeps it finite.

#### 11.10.7 New references from this chapter

| | |
|---|---|
| **Albanese & Rubinacci (1988)** | *"a seminal paper"* — the origin of tree–cotree, earlier than the Manges & Cendes 1995 generalization in §11.6. |
| **Klis (2015)**, *Schnelle Finite-Elemente-Verfahren für verschiedene Klassen magneto-quasistatischer Probleme*, dissertation, Saarland Univ., [doi:10.22028/D291-23134](https://doi.org/10.22028/D291-23134) | **The source of §11.10.4's tree rules.** Open access, German, now in `APhi_Papers/` as `Klis2015_PhDThesis_FastFEM_MagnetoQuasistatic_DE.pdf`. Uses the **A‑V‑A formulation** — i.e. an A‑Φ method, not an E-field one, so its tree treatment should transfer to us more directly than Stysch's does. |
| **Dyczij-Edlinger & Bíró (1996)**, *A joint vector and scalar potential formulation for driven high frequency problems using hybrid edge and nodal finite elements*, IEEE T‑MTT **44**(1) | cited `[DB96]` for **gauged *and* ungauged non-lossy potential formulations**. "Joint vector and scalar potential" with "hybrid edge and nodal finite elements" is *exactly our discretization*, which makes this the most directly comparable prior formulation found so far. Not obtained. |
| **Hiptmair, Kramer & Ostrowski (2008)** | already in §11.6; now known to be one of Eller's two parents, and to give singular non-symmetric matrices. |
| **Jochum et al. (2015)** and **Joc+16** | already in §11.6; now known to be the other parent, PEC-equivalent to Eller, and expensive through matrix fill. |
| **Eller (2017)**, *A Low-Frequency Stable Maxwell Formulation in Frequency Domain and Industrial Applications*, dissertation, TU Darmstadt, [tubiblio/92810](http://tubiblio.ulb.tu-darmstadt.de/92810/) | the long form of the SIAM paper. **Record only — no full text on TUprints**, and per the thesis it omits the tree-creation details too, so it is not worth chasing either. |
| **Kruskal (1956)** | the spanning-tree algorithm they use. |

---

## Bibliography

Papers obtained and read for this edition, all now in `APhi_Papers/`:

| | what it contributes |
|---|---|
| **Chew (2014)**, PIER **149**, 69–84 | The generalized gauge `∇·(εA) = −χ∂Φ/∂t`, `χ = αε²μ`; interface conditions embedded in the PDEs; generalized Green's theorem, extinction theorem, equivalence principle and SIEs for A‑Φ. **The keystone.** |
| **Sharma & Triverio (2022)**, [arXiv:2112.07360](https://arxiv.org/abs/2112.07360) | Modified Lorenz gauge inside lossy conductors; explicit statement that the Φ reference contaminates `∇·A`. DC to high frequency, surface-only conductor modelling. |
| **Sharma & Triverio (2021)**, [arXiv:2108.02764](https://arxiv.org/abs/2108.02764) | Companion: potential-based BEM for lossy materials. |
| **Abdrabou & Gomez (2026)**, [arXiv:2507.02099](https://arxiv.org/abs/2507.02099) | DEC + SIE under Lorenz gauge; no discrete gauging because `A` is three 0‑forms; 14→2 operator reduction; free-space coupling-surface constraint. |
| **Li, Fu & Shanker (2017)**, [arXiv:1705.00265](https://arxiv.org/abs/1705.00265) | Potential integral equations, Lorenz gauge, decoupling gives low-frequency stability; second-kind, well-conditioned operators. |

| **Ansari, Farquharson & MacLachlan (2017)**, GJI **210**(1), 105–129 | Explicit Lagrange-multiplier Coulomb gauge. Demonstrates non-unique potentials **by showing direct and iterative solvers disagree on Φ and A while agreeing on the fields**. Uses edge elements for `A` and nodal for Φ — our discretization. Identifies the cause as `n̂·A` being discontinuous across interfaces while the Coulomb gauge requires it continuous (§6.3.1). |
| **Rapetti, Alonso Rodríguez & De Los Santos (2022)**, *J* **5**(1), 52–63 | High-order tree–cotree in general domains. States outright that **tree gauges are not a discretization of the Coulomb gauge and enforce no orthogonality**, and that the tree degrees of freedom are **set arbitrarily** without affecting `B`. The structural reason tree-gauged Φ has no continuum limit. |

Previously read, retained from the first edition:

| | |
|---|---|
| **Munteanu**, *Tree-cotree condensation properties* | Variants A–E and their conditioning. §IV: non-uniqueness of **A** is "in principle irrelevant"; gauging exists for regularity. §V.E: tree-independence possible, costs sparsity. **A-only, magnetostatic — contains no Φ and cannot speak to its uniqueness.** |
| **Li, Sun, Dai & Chew (2015)**, IEEE Trans. Magn. | Generalized Coulomb gauge (static limit of the Lorenz one). Ungauged systems let iterative solvers converge to different answers; constant gradients live in the edge space. |
| **Li, Sun, Dai & Chew (2016)**, IEEE TAP **64**(10) | FEM implementation of the generalized-Lorenz A‑Φ formulation. Two null-space-free Helmholtz equations; conditioning ~1e6 against 1e17. Needs Hodge operators and a sparse approximate inverse. |
| **Zhao & Fu (2017)**, IEEE Trans. Magn. **53**(6) | Coulomb gauge via a dummy scalar proven zero; symmetric, direct-solvable, 1.38 M unknowns. Good fit, wrong gauge family. |
| **Lee & Jin (2008)**, MOTL **50**(6) | Tree–cotree conditioning; splitting is not unique; root choice matters; conditioning still degrades at low frequency. |
| **Su Yan (2021)**, PIER M **106** | Auxiliary scalar potential to eliminate tree–cotree graph searches while keeping all-frequency stability. |

**No paper in this collection states that tree–cotree yields a non-unique Φ.**
The nearest are Munteanu's tree-dependence remark (about **A**) and Sharma &
Triverio's reference-dependence of `∇·A`. Our 15 % measurement in
`07_GaugeInvariance` remains our own result and should keep being presented as
such.
