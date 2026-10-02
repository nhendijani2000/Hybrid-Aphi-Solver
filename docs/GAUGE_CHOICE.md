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
| **Is anything broken today?** | No. `E`, `B`, `H`, `J`, terminal `V`, `I`, `R`, `L`, `Z` are correct and gauge-invariant — `07_GaugeInvariance` measures 1e‑12. |
| **What is wrong** | Φ in the interior is not physical; tree–cotree postpones rather than removes low-frequency breakdown; and it offers nothing at a conductor/dielectric port. |
| **Still open** | Whether `A` stays a Whitney 1‑form at the boundary-integral surface. This is now the real research question — see §6. |
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
    n̂ × A₁ = n̂ × A₂                              (23)
    n̂ × (1/μ₁)∇×A₁ = n̂ × (1/μ₂)∇×A₂   ≡  n̂ × H₁ = n̂ × H₂     (24),(25)
```

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

**This is the closest published analogue to case 07's result** — the potential
*reference* contaminating the vector potential, where ours is the *tree* doing
it. The disease is the same: a bookkeeping choice leaking into a solved
quantity. It is worth citing in the case 07 README.

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

---

## 6. The open question is no longer which gauge

It is **1‑form vs 0‑form `A` at the boundary-integral surface.**

The fourteen operators arise *specifically* from treating `A` as a 1‑form at the
boundary. Abdrabou & Gomez's answer was to stop doing that. Our plan is to couple
**Whitney 1‑form `A`** directly to a potential BEM — i.e. to take on exactly the
problem they chose to avoid.

That is simultaneously the project's novelty claim and its principal technical
risk, and the proposal should say so rather than treating the fourteen-operator
trap as something the gauge disposes of.

Chew's §6 offers a middle path worth evaluating: his equivalence principle needs
six scalar surface quantities, not fourteen operators, and his surface currents
are divergence-conforming — which is closer to Whitney-compatible than a
Cartesian 0‑form decomposition.

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
implemented. `Ansari, Farquharson & MacLachlan (2017)` and `Rapetti, Alonso
Rodríguez & De los Santos (2022)` are still absent from `APhi_Papers/`; the first
is the roadmap's required benchmark comparison.

---

## 8. Candidates, final ranking

| | verdict |
|---|---|
| **Chew generalized gauge** | **Recommended.** Only candidate where FEM, interface conditions and BEM are one formulation. Costs: Hodge/Whitney-mapping machinery, sparse approximate inverse, and the Φ-space question below. |
| Zhao & Fu dummy-variable Coulomb | **Withdrawn.** Symmetric, direct-solvable and a good architectural fit — but Coulomb-gauged, so it cannot couple to a potential BEM. Implementing it means implementing a gauge twice. |
| Ansari explicit Lagrange-multiplier Coulomb | Same Coulomb objection. Retain as the roadmap's benchmark comparison; the paper is still missing. |
| Improve tree–cotree (root choice, MOR) | Cheap, treats the symptom only. Lee & Jin themselves report the condition number still grows at low frequency. |
| Do nothing | Viable *only* while the outer boundary stays `n × A = 0`. §1 explains why that ends. |

**A design decision the roadmap does not currently mention:** Chew's construction
places Φ in the Whitney‑0 space, i.e. **P1 nodal**. **Ours is P2.** The
compatibility argument rests on the Whitney/de Rham structure, so adopting the
method means either moving Φ to P1 or redoing the analysis for P2.

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

## Bibliography

Papers obtained and read for this edition, all now in `APhi_Papers/`:

| | what it contributes |
|---|---|
| **Chew (2014)**, PIER **149**, 69–84 | The generalized gauge `∇·(εA) = −χ∂Φ/∂t`, `χ = αε²μ`; interface conditions embedded in the PDEs; generalized Green's theorem, extinction theorem, equivalence principle and SIEs for A‑Φ. **The keystone.** |
| **Sharma & Triverio (2022)**, [arXiv:2112.07360](https://arxiv.org/abs/2112.07360) | Modified Lorenz gauge inside lossy conductors; explicit statement that the Φ reference contaminates `∇·A`. DC to high frequency, surface-only conductor modelling. |
| **Sharma & Triverio (2021)**, [arXiv:2108.02764](https://arxiv.org/abs/2108.02764) | Companion: potential-based BEM for lossy materials. |
| **Abdrabou & Gomez (2026)**, [arXiv:2507.02099](https://arxiv.org/abs/2507.02099) | DEC + SIE under Lorenz gauge; no discrete gauging because `A` is three 0‑forms; 14→2 operator reduction; free-space coupling-surface constraint. |
| **Li, Fu & Shanker (2017)**, [arXiv:1705.00265](https://arxiv.org/abs/1705.00265) | Potential integral equations, Lorenz gauge, decoupling gives low-frequency stability; second-kind, well-conditioned operators. |

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
