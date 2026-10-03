# The generalized-Lorenz gauged A‑Φ formulation

An implementation specification, transcribed from

> Y.-L. Li, S. Sun, Q. I. Dai and W. C. Chew, *Finite Element Implementation of
> the Generalized-Lorenz Gauged A‑Φ Formulation for Low-Frequency Circuit
> Modeling*, **IEEE Trans. Antennas Propag. 64(10), 4355–4364, 2016**,

with the continuum theory from

> W. C. Chew, *Vector Potential Electromagnetics with Generalized Gauge for
> Inhomogeneous Media: Formulation*, **PIER 149, 69–84, 2014**.

Both are in `APhi_Papers/`. `docs/GAUGE_CHOICE.md` is the decision record for
*why* this gauge; this document is *what to build*.

**Status: specification, nothing implemented.** Equations (1)–(37) below are
transcribed from the typeset paper. §7 is **my own derivation and is not in
either paper** — it is flagged as such and must be checked before it is trusted.

> **§11 reports a spike that measured this, and it changes §7 and §8.** The
> formulation conditions beautifully for dielectrics — flat across ten decades —
> and **badly for lossy conductors at low frequency**, where `ε_eff = ε − jσ/ω`
> creates a material contrast of 1e17. `χ` must also be position-dependent, which
> breaks §7's symmetrization. **Read §11 before implementing anything.**

---

## 0. The convention trap — read this first

**Chew uses `e^{−iωt}`. We use `e^{+jωt}`.** His eq. (1) is `E = iωA − ∇Φ`; ours
is `E = −jωA − ∇Φ`.

```
                every  i  in the paper   ->   −j  here
```

Our convention is fixed by `src/element_matrix.cpp` (`α = jωσ − ω²ε`) and
`src/problem.cpp` (`A(cos φ + j sin φ)`), and is stamped into every `.out`
header. **Nothing below is a transcription error waiting to happen except this**,
so every equation in this document is already written in *our* convention, and
anyone cross-checking against the paper must apply the mapping above rather than
reading the signs off the page.

The purely real matrices in §5 are unaffected — they contain no `i`.

---

## 1. The continuum formulation, in our convention

```
    E  =  −jωA − ∇Φ                                              (1)

    ∇·(ε∇Φ) + χω²Φ  =  −ρ                                        (2)

    ∇×(ν ∇×A) − ω²εA − ε∇[ (1/χ) ∇·(εA) ]  =  J                  (3)

    χ  =  α μ ε²                                                 (4)

    (1/χ) ∇·(εA)  =  −jωΦ          the generalized-Lorenz gauge   (5)
```

with `ν = 1/μ` as in our existing code.

**(2) and (3) are completely decoupled** — (2) contains no `A`, (3) contains no
`Φ`. That is the whole point of this gauge: two independent systems instead of
one coupled block system, with `E` recovered afterwards from (1). (2) is a scalar
Helmholtz–Poisson equation; (3) is a vector Helmholtz equation carrying a gauge
term.

**On α.** Chew's §II.A allows `α` to be a function of position and notes that all
eigenvalues are positive for `α > 0`, so (3) yields a positive-definite system.
**He recommends `χ = με²`, i.e. `α = 1`**, for the smallest condition number.

> **An inconsistency to resolve before coding.** §II.A defines `χ = αμε²`, which
> is position-dependent through `μ` and `ε`. §II.C then says *"Assuming
> `χ = μ₀ε₀²`"* and carries that constant into (21), (24) and (27). So the
> **discretization as published uses a constant, free-space `χ`**, not the
> position-dependent one. Both are defensible — the gauge is a free choice — but
> they are different gauges and we must pick one deliberately. §8 lists this as
> an open decision.

### Excitation

```
    ∇·J  =  −jωρ                                                 (13)

    ∇·(ε∇Φ) + χω²Φ  =  ∇·J / (jω)                                (14)
```

(14) is (2) with `ρ` eliminated through continuity, and is how the Φ-equation is
driven in a current-driven problem. Chew's remark is worth keeping: **(3)
resembles the current–voltage relation of an inductor and (14) that of a
capacitor** — the magneto-quasi-static and electro-quasi-static halves, solved
separately and to the same accuracy.

### Boundary conditions

```
    Φ = Φs                      r ∈ Γ_D       Dirichlet            (9)
    n̂·ε∇Φ = −ρs                 r ∈ Γ_N       Neumann             (10)
    n̂×A = As                    r ∈ Γ_D       Dirichlet           (11)
    n̂×(ν ∇×A) = Js              r ∈ Γ_N       Neumann             (12)
```

Chew notes that **with Whitney‑0 for Φ and Whitney‑1 for A, the conditions on the
material interfaces are satisfied automatically** — which is the property that
makes a mixed conductor/dielectric port tractable, and the reason this gauge was
chosen. See `GAUGE_CHOICE.md` §2.2.

### Conductors

**Chew 2016 contains no `σ`.** Its conductors are PEC, imposed as boundary
conditions; the volumes are dielectric. Finite conductivity enters the usual way,
by folding it into a complex permittivity — in our convention

```
    ε_eff  =  ε − j σ / ω
```

so that `jω ε_eff E = (σ + jωε) E = J_total`. Everything below then holds with
`ε → ε_eff`. **This is not in the paper** and has one consequence worth flagging
early: `χ = αμε²` becomes complex, and whether the positive-definiteness argument
of §II.A survives that is not established. §8.

---

## 2. Why the gauge term cannot be discretized directly

The third term of (3) applies a divergence to `A`. Discretized naively on
Whitney‑1 elements it **vanishes identically**, because

```
    ∇·ω_mn  =  ∇·( λ_m ∇λ_n − λ_n ∇λ_m )  =  0                    (15)
```

so a Ritz or Galerkin discretization of that term gives zero matrix entries and
(3) silently degenerates to the ungauged double-curl equation.

But Whitney‑1 elements are **not** solenoidal, because a constant gradient lies
in their span:

```
    ∇λ_m  =  Σ_n ω_κm                                            (16)
```

(the same fact that makes the tree–cotree null space exist at all, and the one
Chew 2015 flags as the reason non-uniqueness survives curl-conforming elements).

**The resolution is an intermediate quantity.** Define

```
    d  :=  (1/χ) ∇·(εA)
```

and expand *it* in Whitney‑0 elements. Chew's Fig. 2 is the justification, and it
is a differential-forms argument:

```
    A  --*ε-->  εA  --∇·-->  ∇·εA  --*μ,*ε-->  (1/με²)∇·εA
  1-form      2-form        3-form              0-form
```

so `d` lands in the same space as `Φ`, exactly as (5) requires. Table I of the
paper:

| quantity | space | geometric constituent |
|---|---|---|
| `Φ` | 0-form, potential space | primal nodes |
| `A` | 1-form, field space | primal edges |
| `E` | 1-form, field space | primal edges |
| `(1/με²)∇·εA` | 0-form, potential space | primal nodes |

**No dual mesh is involved**, which is what keeps this compatible with an
ordinary primal-mesh FEM code.

---

## 3. Expansions

```
    Φ  =  Σ_{n=1}^{Nn}  φ_n λ_n(r)          Whitney-0 (P1 nodal)   (17)

    A  =  Σ_{n=1}^{Ne}  a_n ω_n(r)          Whitney-1 (edge)       (18)

    d  =  Σ_{n=1}^{Nn}  d_n λ_n(r)          Whitney-0 (P1 nodal)   (19)
```

> **Φ is P1 here. Ours is P2** — `src/element_matrix.cpp` builds a 10×10 Φ‑Φ
> block over 4 vertices plus 6 edge midpoints. The de Rham/Whitney argument above
> is what licenses the whole construction, and it is a statement about Whitney‑0.
> Moving Φ to P1 is the main discretization change this gauge asks for. §8.

Chew notes (20) that one could instead expand `(1/ε)∇·εA`, but then the
divergence theorem is needed and a surface integral term appears; he keeps (19)
because it makes the third term of (3) integrable over the domain and yields a
simple matrix representation.

---

## 4. Eliminating `d`

Testing (19) with `λ_m` and using `χ = μ₀ε₀²`:

```
    Σ_n ⟨ λ_m , (1/χ) ∇·(ε ω_n) ⟩ a_n  =  Σ_n ⟨ λ_m , λ_n ⟩ d_n   (21)

    d  =  Ḡ⁻¹ K_NE a                                             (22)
```

---

## 5. The matrices

All real; none of them carries the time convention.

```
    [Ḡ]_mn     = ∫_Ω λ_m λ_n dΩ                                              (23)

    [K_NE]_mn  = ∫_Γ (ε/χ) λ_m ω_n·n̂ dΓ  −  ∫_Ω (ε/χ) ∇λ_m·ω_n dΩ           (24)

    [S1]_mn    = ∫_Ω ∇λ_m·ε∇λ_n dΩ  −  ∫_Γ λ_m n̂·ε∇λ_n dΓ                   (26)

    [S2]_mn    = ∫_Ω χ λ_m λ_n dΩ                                            (27)

    [b_s]_m    = ∫_Ω λ_m ρ dΩ                                                (28)

    [K1]_mn    = ∫_Ω ν (∇×ω_m)·(∇×ω_n) dΩ                                    (30)

    [K2]_mn    = ∫_Ω ε ω_m·ω_n dΩ                                            (31)

    [K_EN]_mn  = ∫_Ω ε ω_m·∇λ_n dΩ                                           (32)

    [b_k]_m    = ∫_Ω ω_m·J dΩ                                                (33)
```

### The two systems

```
    Φ-equation:   ( S1 − ω² S2 ) φ  =  b_s                                   (25)

    A-equation:   ( K1 − ω² K2 − K_EN Ḡ⁻¹ K_NE ) a  =  b_k                   (29)
```

### Dimensions

| | |
|---|---|
| `S1, S2, Ḡ` | Nn × Nn |
| `K1, K2` | Ne × Ne |
| `K_EN` | Ne × Nn |
| `K_NE` | Nn × Ne |
| `φ, b_s` | Nn × 1 |
| `a, b_k` | Ne × 1 |

For a large 3‑D mesh `Nn/Ne ≈ 1/7.3`, so the A-system costs about **13.7 %** more
unknowns than an E formulation. Chew measures 15.04 % and 15.35 % on his two
examples.

### Static bounds, for preconditioning and for sanity

```
    K  =  K1 − ω²K2 − K_EN Ḡ⁻¹ K_NE          K' =  K1 − K_EN Ḡ⁻¹ K_NE    (35),(37)
    S  =  S1 − ω²S2                          S' =  S1                    (34),(36)
```

`K'` and `S'` are the `ω = 0` limits and are **frequency-independent**. Chew's
reported condition numbers: `K` stays bounded by `K'` at `5.94e5` and `1.21e6`
across his two cases and the whole frequency range, against the E formulation's
`1.04e17` at 10 kHz.

---

## 6. Solving — direct is explicitly supported

Worth recording plainly, because the first reading of this paper suggested
otherwise. §II.D: *sparse direct solvers such as UMFPACK and CHOLMOD are suitable
candidates* for problems of tens of thousands of unknowns, and §III states that
**UMFPACK multifrontal LU is applied first** to make the solution as accurate as
possible. The ILU/GMRES material is an optimization for wideband sweeps — the
frequency-independent `K'` can be factored once and reused at every frequency —
**not a requirement of the method.**

So this composes with our MUMPS direct path rather than fighting it.

---

## 7. Avoiding the sparse approximate inverse — our variant

> **This section is my derivation. It is not in either paper and has not been
> tested.** It must be checked symbolically and numerically before it is relied
> on. Everything above is transcribed; everything here is not.

`Ḡ⁻¹` is the problem. `Ḡ` is the nodal mass matrix, so `Ḡ⁻¹` is dense, and
`K_EN Ḡ⁻¹ K_NE` formed explicitly would destroy sparsity. Chew's answer is a
**sparse approximate inverse** — on his examples `Ḡ⁻¹` is replaced by an SAI with
9057 and 7147 nonzeros, 98.63 % and 99.98 % sparse. That is an approximation
adopted to keep the system `Ne × Ne` for iterative solving.

**We factor directly, so we can avoid it entirely** by keeping `d` as unknowns.
From (22), `Ḡ d = K_NE a`, and (29) becomes `(K1 − ω²K2) a − K_EN d = b_k`:

```
    ⎡ K1 − ω²K2    −K_EN ⎤ ⎡ a ⎤   ⎡ b_k ⎤
    ⎢                    ⎥ ⎢   ⎥ = ⎢     ⎥
    ⎣ K_NE          −Ḡ   ⎦ ⎣ d ⎦   ⎣  0  ⎦
```

Everything in it is sparse and exact. The cost is `+Nn` unknowns — about 14 % on
top of the A-system, i.e. roughly the same again as the formulation already costs
over an E formulation.

**It can be made symmetric.** Ignoring the surface term, (24) is

```
    [K_NE]_mn  =  −(1/χ) ∫_Ω ε ∇λ_m·ω_n dΩ  =  −(1/χ) [K_EN]_nm
```

so `K_NE = −(1/χ) K_ENᵀ`. Scaling the second block row by `χ`:

```
    ⎡ K1 − ω²K2    −K_EN ⎤ ⎡ a ⎤   ⎡ b_k ⎤
    ⎢                    ⎥ ⎢   ⎥ = ⎢     ⎥        symmetric
    ⎣ −K_ENᵀ       −χ Ḡ  ⎦ ⎣ d ⎦   ⎣  0  ⎦
```

which is what MUMPS's complex-symmetric path wants.

**Three things to check before trusting this:**

1. **The surface term of (24)** — `∫_Γ (ε/χ) λ_m ω_n·n̂ dΓ` — is dropped above.
   It must either vanish under our boundary conditions or be carried, and if
   carried it breaks the exact symmetry.
2. **`χ` inside the scaling** must be constant for the row scaling to be a
   congruence. If we adopt a position-dependent `χ` (§1), this derivation needs
   redoing.
3. **Conditioning.** The block system is larger and has a zero-ish structure in
   the `(2,2)` corner; whether its condition number behaves as well as Chew's
   reported `K` is an empirical question, and is the first thing the spike should
   measure.

If any of these fails, Chew's SAI route remains available and is the fallback.

---

## 8. Open decisions

| | |
|---|---|
| **`χ` constant or position-dependent?** | §II.A defines `χ = αμε²`; §II.C discretizes with `χ = μ₀ε₀²` constant. Constant makes §7's symmetrization a clean congruence. Position-dependent is the "generalized" gauge as stated. **Pick deliberately and record why.** |
| **Φ: P2 → P1?** | The Whitney/de Rham argument is about Whitney‑0. Our Φ is P2. Either move Φ to P1 or redo the compatibility analysis for P2 — the latter is research, not implementation. |
| **Lossy conductors** | Fold `σ` into `ε_eff = ε − jσ/ω`. Not in the paper. Makes `χ` complex; the positive-definiteness argument of §II.A is then unproven. |
| **`Ḡ⁻¹`: block system or SAI?** | §7 versus Chew's SAI. Decide on measured conditioning, not on preference. |
| **Per-region boundary conditions on a port face** | Needed for the mixed port, which is the acceptance test (§10). The formulation supplies it through the `Γ_D`/`Γ_N` partition of (9)–(12); our parser rejects per-surface BCs today. Parser and binding work, independent of the gauge, and can proceed in parallel. |

---

## 9. What maps onto existing code

From `include/aphi_solver/element_matrix.hpp`, we already build per-tet:

| ours | theirs | status |
|---|---|---|
| A‑A curl part, `ν ∫(∇×W_i)·(∇×W_j)` | **`K1`** (30) | **identical** |
| A‑A mass part, `α ∫W_i·W_j`, `α = jωσ − ω²ε` | **`−ω²K2`** (31) | same shape; our `α` already carries `σ` via the complex-`ε` route |
| A‑Φ coupling, `β ∫W_i·∇S_b`, `β = σ + jωε` | **`K_EN`** (32) | same shape, different coefficient and Φ space |
| Φ‑Φ, `β ∫∇S_a·∇S_b` | **`S1`** (26) | same shape, coefficient `ε` not `β` |

**New work:**

- **`Ḡ`** (23) — nodal mass. Trivial; the quadrature exists.
- **`S2`** (27) — nodal mass weighted by `χ`. Trivial given `Ḡ`.
- **`K_NE`** (24) — **the genuinely new one.** A volume term that is `K_ENᵀ` up to
  the `ε/χ` weighting, plus a **surface** term over `Γ`. Our assembly has no
  surface-integral path for this block today.
- **The gauge term's coupling** — either the block system of §7 or an SAI.
- **Decoupling the solve** — (25) and (29) are independent, so the driver changes
  from one factorization to two.
- **Removing tree–cotree** from the A path.

**What goes away:** `tree_cotree.cpp`, `gauge_variants.cpp` and the gauging in
`dof_map.cpp` cease to be load-bearing. Keep `07_GaugeInvariance` — it becomes the
test that Φ *no longer* moves.

---

## 10. Validation plan

### The acceptance test is the mixed-material port

**This is why the gauge is being replaced at all.** `ROADMAP.md` Phase 03 names a
port face with conductor on one side and dielectric on the other, inside the
domain, as the reason tree–cotree is inadequate; it is corroborated by the
project owner's experience of a commercial solver running such a port to
completion and returning wrong answers (`GAUGE_CHOICE.md` §1.3). Everything else
in this list is a necessary condition. **This one is the goal**, and the gauge
has not succeeded until it passes.

**The formulation already provides what a mixed port needs**, which the first
draft of this section under-stated. Eqs (9)–(12) partition the boundary into
`Γ_D` and `Γ_N` **separately for Φ and for A**. A mixed port face is exactly that
partition:

| part of the face | Φ | A |
|---|---|---|
| conductor | `Φ = Φs` — Dirichlet (9) | `n̂×A = As` — Dirichlet (11) |
| dielectric | `n̂·ε∇Φ = −ρs` — Neumann (10) | `n̂×(ν∇×A) = Js` — Neumann (12) |

and Chew's remark that with Whitney‑0 Φ and Whitney‑1 A the **interface
conditions are satisfied automatically** is what makes the conductor/dielectric
boundary inside the face unproblematic.

> **The gap is on our side and it is bounded.** `src/input_file.cpp` currently
> fails with *"per-surface boundary conditions are not supported yet"*, and every
> port type we have imposes a single potential over the whole tagged surface. So
> a mixed port needs **per-region boundary conditions on a port face** — Dirichlet
> on the conductor part, natural on the dielectric part. That is parser and
> binding work, not a reformulation, and it is independent enough to proceed in
> parallel with the gauge.

### The baseline cannot be built yet, and that sets the order

**A mixed port is not expressible by this solver today.** `[boundary]` "takes no
name and applies to the whole outer boundary" — per-surface conditions are
explicitly unimplemented — and every port type imposes one potential over the
whole tagged surface. So building the coax now would produce one of two useless
results:

- tag the whole cross-section as one port → the dielectric is shorted to the
  conductor's potential → the ill-posed garbage the first case 08 attempt
  produced;
- tag only the conductors → it is not a mixed port, and case 05 has already
  measured what an interior conductor cut does.

**Either way it would measure the port model's inability to pose the problem, not
the gauge's failure to solve it** — the same confound that wrecked the first
attempt, on a larger mesh. "Measure the failure before fixing it" does not rescue
it: a failure cannot be measured on a problem that cannot be stated.

So the dependency is strict rather than parallel:

```
    per-region port BCs   ->   coax baseline   ->   mixed-port acceptance test
```

> **The first case 08 attempt is worth remembering here.** A straight rod with an
> internal cut has no return path, so the system was ill-posed; both the mixed
> face and the conductor-only control came back wrecked, and it was the *control*
> that showed the geometry was at fault rather than the gauge. Any rebuild needs
> a closed circuit — coax, up the centre and back through the shield.

### In order

Steps 1 and 2 are genuinely independent and can run in parallel. Nothing after
them can start early.

1. **Spike the gauge term** — assemble `K_NE` and `Ḡ` on `07_GaugeInvariance`'s
   8,365-tet mesh. Is the §7 block system complex-symmetric, does MUMPS factor
   it, and how does its condition number compare with the tree-gauged system?
   **No new geometry, no port work, commits to nothing, and it retires the largest
   unknown in the plan.**
2. **Per-region boundary conditions on a port face** — Dirichlet on the conductor
   part, natural on the dielectric part, shaped by the `Γ_D`/`Γ_N` partition of
   (9)–(12). Parser and binding work, not reformulation. **On the critical path
   regardless of the gauge**, because without it there is no acceptance test.
3. **Case 08 built and failing** — the coax mixed port, demonstrated broken under
   tree–cotree, with numbers. Nearly free once step 2 exists.
4. **Reproduce Chew's Fig. 5** qualitatively — condition number flat against
   frequency where the E formulation blows up. His `5.94e5` and `1.21e6` are the
   target order.
5. **`07_GaugeInvariance` inverted** — Φ must now be mesh-independent. The
   permutation harness stays; the *potential-moves* assertions flip.
6. **Case 05's local field** — the 8 % tree-dependence of `J` within 2 mm of the
   cut (`GAUGE_CHOICE.md` §1.3) must go away. A necessary condition, and the only
   one measurable on an existing case today.
7. **Cases 01–06 unchanged** on every field quantity. This is a reformulation, not
   a change of physics.
8. **The mixed port passes** — step 3's baseline, now correct. **Tree-independence
   is not the test here**: once the tree is gone there is no tree to depend on, so
   the standard is that the local field on and near the port face **converges
   under mesh refinement**.

---

## 11. Spike results — measured, and they change §7 and §8

`tools/gauge_spike/` assembles §5's matrices on 518 tets and computes exact
condition numbers by dense SVD. It ran before any C++ was written, which is the
only reason the findings below cost a day rather than a month.

### The formulation does what Chew claims — for dielectrics

His own materials: dielectric everywhere, no `σ`.

| freq (Hz) | cond gauged | cond ungauged |
|---|---|---|
| 1e0 | **2.492e+02** | 1.884e+18 |
| 1e4 | **2.492e+02** | 2.584e+14 |
| 1e10 | 3.907e+02 | 5.729e+02 |

**Flat across ten decades** while the ungauged system degrades by sixteen orders
as `ω → 0`. That is Fig. 5 of the paper reproduced in character, and it is what
validates the assembly — every other number here depends on it.

### It does not survive `σ` folded into `ε`

| freq (Hz) | cond gauged | cond ungauged | ε contrast |
|---|---|---|---|
| 1e0 | **2.046e+22** | 1.655e+06 | 2.317e+17 |
| 1e4 | 9.272e+17 | 5.745e+02 | 2.317e+13 |
| 1e10 | 1.527e+10 | 1.118e+09 | 2.317e+07 |

The condition number **tracks the ε contrast** at about `contrast × 1e4`, and
below roughly 1 GHz **the gauged system is worse than the ungauged one**.

Not a bug: `ε_eff = ε − jσ/ω` gives `σ/ωε = 2.3e17` between copper and dielectric
at 1 Hz. **Chew has no finite conductivity** — PEC boundary conditions, dielectric
volumes — so his `ε` varies by 4.5. `ε_eff` is §1's addition, and it is what
breaks the conditioning.

> **This qualifies the roadmap's low-frequency rationale, sharply.** Phase 03
> wants this gauge partly to remove the low-frequency floor. For a *dielectric*
> problem it does, completely. For **lossy metal at low frequency it makes
> conditioning worse**, monotonically in `σ/ωε`.
>
> It does not invalidate the choice — we want this gauge for the **embedded
> interface conditions** at a mixed port (§10), which the above does not touch.
> But "fixes low-frequency breakdown" must not be claimed for conductors without
> qualification.
>
> **At THz the problem recedes.** The contrast falls as `1/ω`: 2.3e7 at 10 GHz and
> about 2.3e5 at 1 THz, where conditioning should be unremarkable. The difficulty
> is specific to the low-frequency end with good conductors.

### `χ` must be position-dependent — and that breaks §7's symmetry

| | cond at 1 Hz |
|---|---|
| `χ = μ₀ε₀²` constant — his §II.C | 7.026e+37 |
| `χ = με²` per element — his §II.A | 2.046e+22 |

Fifteen orders, because `ε²/χ = 1/μ = ν` self-normalizes to the curl–curl
coefficient **only when `χ` carries the local `ε`**. §8's open question is
therefore answered: **position-dependent**.

But §7's symmetrization is a global row scaling by `χ`, a congruence only for
constant `χ`. Measured per element, `K_NE` is no longer a scalar multiple of
`K_ENᵀ` — residual **O(1)**, not 1e‑16.

**As derived, the block system is symmetric or well-conditioned, not both.** The
options are a per-element scaling that is still a congruence, an unsymmetric
system, or Chew's SAI route.

> **Qualified by `GAUGE_CHOICE.md` §11.4, added after this section.**
> Balian et al. (2023) give the congruence explicitly — scale the
> equation *and* the unknown by the same `(σ + jωε)^(−1/2)`, so symmetry
> survives. It does **not** fully rescue us: `ε/χ` sits inside the
> element integral, so a diagonal congruence reproduces
> `K_NE = −c·K_ENᵀ` for every node interior to one material and **fails
> at material-interface nodes** — exactly where a mixed port lives. The
> literature's answer is to assign interface DOFs to the conductor block
> and scale blockwise, rather than to seek exactness.

### What §7 got right, and what it got wrong

**Right:** with constant `χ` the block system is exactly complex-symmetric —
measured 1e‑21 to 8e‑17. The derivation `K_NE = −(1/χ)K_ENᵀ` holds.

**Wrong:** §7 ignored scaling. `χ ≈ 1e‑28` puts the `(2,2)` block 1e40 below the
`(1,1)` block, and the first spike run reported 1e46 because of it. The scaling
`c` is free — it cancels in the Schur complement — so it can be chosen for
conditioning, and `c² = ‖K_uu‖ / (χ‖G‖)` recovered 26 orders. **Any block form
must be scaled deliberately**; the naive one is unusable.

### Why the controls are the whole point

The first two spike runs concluded the formulation was unusable, at 1e38. Both
were wrong — the first was my block scaling, the second was finite conductivity,
which is not in the paper. Only running **Chew's materials unchanged** showed the
assembly had been right the whole time.

Without that control this document would now be recommending against a sound
method.

### Revised open questions

| | |
|---|---|
| **How to handle `σ`** | `ε_eff` is measured to wreck conditioning below ~1 GHz. **This is a known, published phenomenon with published remedies** — see `GAUGE_CHOICE.md` §11, added after this section was written. Balian et al. (2023) name the same `σ/ωε` mechanism and fix it by scaling equations *and* unknowns analytically before assembly; Herles et al. (2025) use **a different gauge condition inside conductors than outside**, which is the "different `χ` inside conductors" line below, built and measured. Not a blocker. |
| **Symmetry with per-element `χ`** | As above. Symmetric or well-conditioned, not both, as currently derived. |
| **The surface term of (24)** | Still omitted in the spike. Must vanish under our boundary conditions or be carried; carrying it may break symmetry independently of `χ`. |
