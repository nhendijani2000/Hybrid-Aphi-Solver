# A terminal on Γ in Jochum's formulation — derivation report

**Status: §6 has been tested. It gets the physics right and the gauge wrong.**
§§1–5 are worked and stand. §6's construction reproduces the correct terminal
impedance at every frequency — and is **as tree-dependent as the tree–cotree
solve it was meant to replace**. §10 has the measurement and the diagnosis.
Nothing is implemented in C++.

Why this document exists: `tools/gauge_spike/RESULTS_SPIKEC.md` milestone 2
found that Jochum, Farle & Dyczij-Edlinger's formulation — the only one in the
literature carrying a materially-weighted gauge interface condition at a
conductor/insulator boundary, which is exactly where our failing port sits — has
**nowhere to put a port**. This works out why, and what would have to change.

---

## 1. The formulation, in our convention

Jochum writes `e^{+ik₀·}` with a scaled potential; translated to ours
(`e^{+jωt}`, `E = −jωA − ∇Φ`), with `κ = σ + jωε` and `ν = 1/μ`:

```
    A  =  A_c + ∇ψ ,       B = ∇×A_c ,                             (15),(16)
    E  =  −∇V − jω(A_c + ∇ψ)                                       (17)
```

Three unknowns: `A_c` on cotree edges, `ψ` nodal, `V` nodal. `A_c` lives in the
reduced space of the **inexact Helmholtz splitting** (14), realized discretely
by tree–cotree — so the tree is still there, but it decomposes the space rather
than fixing the gauge.

In the lossy sub-domain `Ω_C`:

```
    ∇×(ν∇×A_c) + κ[jω(A_c+∇ψ) + ∇V]  =  J_i          Ampère        (18a)
    ∇·[κ(A_c+∇ψ)]                     =  0            GAUGE         (18b)
```

In the lossless sub-domain `Ω_N`:

```
    ∇×(ν∇×A_c) + jωε[jω(A_c+∇ψ) + ∇V] =  J_i          Ampère        (19a)
    ∇·{ε[jω(A_c+∇ψ) + ∇V]}            =  ρ            Gauss         (19b)
```

Since `κ|_Ω_N = jωε`, (19a) is (18a) with the local `κ` — one expression covers
both. The asymmetry is entirely in the **second** equation: a gauge in the
conductor, physics in the insulator.

---

## 2. Which equation lives in which row — the fact everything turns on

Jochum tests each equation against a different space (his §5.2):

| weak form | equation tested | test function | produces the row for |
|---|---|---|---|
| (21) | Ampère (18a) | `w_c ∈ Q̃H(curl)` | **`A_c`** |
| (22) | Ampère (18a) | **`∇ψ̄`**, `ψ̄ ∈ H¹₀` | **`ψ`** |
| (23) | **gauge (18b)** | `V̄ ∈ H¹₀` | **`V`** |

and states that (22) is a weak form of the **continuity equation**.

So, and this is the crux:

> **`ψ`'s row is continuity. `V`'s row is the gauge.**

That cross-placement is what makes the off-diagonal blocks transposes and the
matrix complex-symmetric — it is a feature, not an accident. But it **inverts
the usual A‑Φ correspondence**, where a current-driven port's row is the
continuity equation and its unknown is the potential.

---

## 3. Deriving (32) without a terminal

Test (18b) by `V̄` over `Ω_C` and integrate by parts:

```
    −∫_Ω_C κ(A_c+∇ψ)·∇V̄  +  ∮_∂Ω_C V̄ n̂·[κ(A_c+∇ψ)]  =  0
```

Test (19b) by `V̄` over `Ω_N`, noting `ε[jω(A_c+∇ψ)+∇V] = −εE`:

```
    −∫_Ω_N (−εE)·∇V̄  +  ∮_∂Ω_N V̄ n̂·(−εE)  =  ∫ ρ V̄
```

Let `n̂` point from `Ω_C` into `Ω_N`, so `n̂_C = n̂` and `n̂_N = −n̂` on Γ. The two
surface integrals contribute, on Γ,

```
    ∮_Γ V̄ { n̂·[κ(A_c+∇ψ)]_C  +  n̂·(εE)_N } dS                      (★)
```

With `V` single-valued and **free** on Γ, `V̄` is arbitrary there, so (★)
vanishes only if its integrand does pointwise:

```
    n̂·[κ(A_c+∇ψ)]_C  +  n̂·(εE)_N  =  0        on Γ                 (32)
```

**This is Jochum's (32), and it is not an extra assumption — it is the condition
under which the two sub-domain weak forms assemble into one system with no
leftover surface term.** It is materially weighted on both sides: `κ = σ + jωε`
from the conductor, `ε` from the insulator. That weighting is why §12.3 ranks it
above the ε-blind Coulomb-family conditions.

---

## 4. What a terminal is, formally

A current-driven terminal on `Γ_T ⊆ Γ` is two statements:

1. **`V` is single-valued on `Γ_T`** — one unknown `V_T`, not one per node. In
   the discrete setting this is the aggregation `A = P Â Pᵀ` (Stysch §6.6.2,
   `GAUGE_CHOICE.md` §11.10.5), which is a congruence and so preserves symmetry.
2. **The total current through `Γ_T` is prescribed**:

```
    I  =  ∫_Γ_T n̂·(κ E)_C dS                                        (33)
```

Statement 1 collapses the test space: on `Γ_T`, `V̄` is a single constant rather
than an arbitrary function.

---

## 5. The obstacle, located exactly

With `V̄` constant on `Γ_T`, (★) no longer forces its integrand to vanish there.
It forces only the **integral**:

```
    ∫_Γ_T { n̂·[κ(A_c+∇ψ)]_C + n̂·(εE)_N } dS  =  0                   (34)
```

That is one scalar equation, which is right — one equation for one unknown
`V_T`. But now look at what it contains against what we need to impose:

| | normal flux appearing |
|---|---|
| **(34)**, the `V_T` row | `n̂·(κ **A**)_C` — the **gauge** flux |
| **(33)**, the current | `n̂·(κ **E**)_C` — the **current** |

**These are different physical objects.** `A` is not `E`. The gauge equation's
natural boundary datum is the normal component of `κA`, which is a statement
about the gauge, and no choice of right-hand side turns it into a current.

And the current flux `n̂·(κE)` does appear in the weak form — but in the
**`ψ` row**, since (22) is continuity and its surface term is `∮ n̂·(κE) ψ̄`.

So the terminal's two requirements land in two different rows:

```
    the terminal POTENTIAL we want to read   is   V_T
    the row V_T owns                         is   the GAUGE        (34)
    the row carrying the CURRENT             is   ψ's              (continuity)
```

**That is the whole obstacle**, and it is structural rather than a matter of
assembly detail.

### 5.1 It explains the measurement exactly

Spike C milestone 2 put the current drive on `V`'s row — the natural thing to
do, and wrong here. Measured consequence: a terminal potential of
`−1.149425e‑05` V identical across three trees to 2.4e‑14, with **zero
reactance**, **unchanged from 50 Hz to 50 GHz**, and a V-coupling of
`G_epsN = 4.2e‑02` against gauge terms `G_kap = 7.5e+07` — a factor 5.7e‑10.

Every one of those follows from §5: the drive perturbed the gauge condition
instead of injecting current, leaving a decoupled DC conduction problem whose
solution is linear in `z` and therefore exact in P1. The spike did not fail; it
measured this.

### 5.2 Why his paper never hit it

His §6.1 cavity sets `φ = 0` and `φ = 1` on faces, and §6.2 is a voltage-driven
RLC field model. Both are **Dirichlet conditions on `V`**, which are trivial
here — prescribe `V` and the row disappears. **A voltage port in Jochum's
formulation is free; only a current-driven floating terminal is a problem.**
That is a useful narrowing: the gap is one port type, not the whole concept.

---

## 6. Proposed resolution — aggregate `ψ` as well as `V`

**Not verified. The dimension count balances and the physics is plausible; that
is all that can be claimed so far.**

If `Γ_T` needs two conditions — the gauge (34) and the current (33) — then give
it **two unknowns**: aggregate `ψ` on `Γ_T` as well as `V`.

```
    V_T  =  V        on Γ_T      one unknown    row:  aggregated gauge   (34)
    ψ_T  =  ψ        on Γ_T      one unknown    row:  aggregated
                                                      continuity  =  I   (33)
```

**Count.** With `N_T` nodes on `Γ_T`, before: `N_T` `V`-unknowns and `N_T`
`ψ`-unknowns against `N_T` gauge rows and `N_T` continuity rows. After: 2
unknowns against 2 rows. Balanced.

**Legitimacy of aggregating `ψ`.** `ψ` is the gauge function. Forcing it
constant on a surface is a *gauge choice*, not a physical constraint — it is
precisely what contracting a surface in the tree does (§15.9), and what
`n̂ × A = 0` achieves on the outer boundary. So it costs nothing physical.

**Symmetry.** Both aggregations are congruences `P(·)Pᵀ`, applied to the `V` and
`ψ` blocks. A congruence preserves complex symmetry, so the system should remain
symmetric — the property §13.2 identifies as Jochum's main advantage over Chew.

**What `V_T` then means.** `V_T` is the terminal potential, read as the answer.
`ψ_T` is a gauge constant on the terminal, carrying no physical meaning — the
analogue of `ψ = 0` on the outer boundary.

### 6.1 What is not yet established

- **Uniqueness.** Does (33)+(34) determine `(V_T, ψ_T)` uniquely, or is the pair
  singular in some limit? The count balances; that is necessary, not sufficient.
- **Whether (32) survives off `Γ_T`.** The derivation in §3 assumed `V̄` free on
  all of Γ. On `Γ \ Γ_T` it still is, so (32) should hold there — but the two
  regions share an edge, and the rim of `Γ_T` needs checking.
- **The low-frequency limit.** `κ → σ` in the conductor and `κ → 0` in the
  insulator as `ω → 0`. §11.10 measured that an insulator-interior `ψ` row
  vanishes at DC, which is why Stysch's §5.5 restricts `ψ` to the interface.
  `Γ_T` **is** on the interface, so `ψ_T` should survive — but that is an
  expectation, not a measurement.
- **Whether `I` belongs on `ψ_T`'s row with a factor.** (33) is the current;
  (22) is continuity tested by `∇ψ̄`. The aggregated continuity row equals the
  net current through `Γ_T` only up to a sign and possibly a scaling from the
  non-dimensionalization. That needs doing carefully, and `RESULTS_COMPARE.md`
  records what dropping such a factor costs.

---

## 7. How to test it, cheaply

`tools/gauge_spike/spikeC.py` already has everything but this: a validated
miniature (milestone 1: `R` exact to 0.06 %, terminal `Φ` moving 3.8 % with the
tree), a Jochum assembly (milestone 2), and a terminal aggregation.

The change is to `solve_jochum`: add a second aggregated unknown on the terminal
face for `ψ`, move the drive from `V`'s row to `ψ`'s, and re-run. The acceptance
test is unchanged and already written:

| | must be |
|---|---|
| `R` | 11.494253 µΩ — the analytic DC value |
| reactance | present and frequency-dependent, **not** constant across decades |
| terminal `Φ` across three trees | invariant |
| symmetry residual | ~1e‑16 |

Milestone 2's failure signature is specific enough to be diagnostic: if the
reactance is again zero and frequency-independent, the drive is still in the
gauge row and §6 is wrong.

---

## 8. If §6 fails

Two fallbacks, both more expensive:

- **Re-derive (32) with the terminal present from the start**, rather than
  aggregating afterwards. (32) came from requiring the surface terms to cancel;
  with a terminal, the correct statement may be a *modified* interface condition
  on `Γ_T` rather than the same one aggregated. This is the honest version of
  "the derivation", and §6 is the shortcut worth trying first.
- **Take the condition from outside.** A boundary-integral exterior supplies the
  relation between `A` and `Φ` on a coupling surface from the exterior solution,
  which is what `GAUGE_CHOICE.md` §0 wants Chew for. That removes the need to
  invent an interface condition at all, at the cost of the whole BI programme.

---

## 9. Sources

| | |
|---|---|
| **Jochum, Farle & Dyczij-Edlinger (2016)**, *A Symmetric and Low-Frequency Stable Potential Formulation…*, in *Scientific Computing in Electrical Engineering* (SCEE 2014), Mathematics in Industry **23**, 63–71, [doi:10.1007/978-3-319-30399-4_7](https://doi.org/10.1007/978-3-319-30399-4_7) — `APhi_Papers/Jochum_SymmetricLowFreqStablePotentialFormulation.pdf` | Equations (14)–(19), (21)–(23), (32); §5.2's test-function assignment, which §2 above turns on; §6's examples, which are voltage-driven and so never meet this problem. |
| **Stysch (2022)**, PhD thesis, TU Darmstadt | §6.6.2's aggregation `A = P Â Pᵀ`, which §4 and §6 use for the terminal; §5.5's restriction of `ψ` to the interface, which §6.1 relies on. |
| `tools/gauge_spike/RESULTS_SPIKEC.md` | Milestone 2's measurement, which §5.1 explains. |
| `GAUGE_CHOICE.md` §13 | The reading of Jochum this builds on. |

---

## 10. Result: §6 tested

`tools/gauge_spike/spikeC.py`, milestone 5. Two findings, and the second
undoes the first.

### 10.1 The readout in §6 was wrong, and fixing it gives correct physics

§6 said `V_T` is the terminal potential. **It is not.** From (17),

```
    E  =  −∇V − jk₀(A_c + ∇ψ)  =  −∇(V + jk₀ψ) − jk₀A_c
```

so the scalar content is **split between `V` and `ψ`**, and the effective
potential — the one that plays the usual `Φ` against `A_c` as the vector
potential — is

```
    Φ_eff  =  V + jk₀ψ                                               (35)
```

Read that way, the construction is **exactly right**:

| freq | Jochum `Φ_eff` | tree–cotree A‑Φ, same quantity |
|---|---|---|
| 50 Hz | 1.150074e‑05 + 1.043363e‑06j | 1.150074461e‑05 + 1.043363113e‑06j |
| 50 kHz | 3.011561e‑04 + 3.567000e‑04j | 3.011561e‑04 + 3.567000e‑04j |
| 50 MHz | 1.212789e‑03 + 3.528559e‑06j | 1.212789e‑03 + 3.528559e‑06j |

Identical to every digit printed, at every frequency — `R` right, **reactance
present and frequency-dependent**, symmetry residual 3.6e‑24. Milestone 2's
failure signature is gone: the drive now injects current instead of perturbing
the gauge, exactly as §5 predicted it would once moved to `ψ`'s row.

So §6's *structure* — two aggregated unknowns, current on the continuity row —
is correct. Reading `V_T` alone was the error, and it explains why milestone 5
first reproduced milestone 2's symptom.

### 10.2 And it is exactly as tree-dependent as what it replaces

| | `Φ_eff`, base→permA | base→permB |
|---|---|---|
| 50 Hz | **3.843e‑02** | 3.533e‑02 |
| 50 kHz | 5.318e‑01 | 5.771e‑01 |
| 50 MHz | 2.823e+02 | 3.246e+02 |

The tree–cotree A‑Φ solve moves the same quantity by **3.84e‑02 and 3.53e‑02**.
Identical. And it degrades with frequency, to a factor of 282 at 50 MHz.

**This fails the acceptance test.** The whole point was a formulation whose
terminal potential does not depend on the tree.

### 10.3 Why — and one hypothesis already falsified

First suspicion: Stysch's §5.5 restriction of `ψ` to the interface and
conductor, which I had carried over. `A = A_c + ∇ψ` spans the full edge space
only if `ψ` spans all nodes — cotree is `#edges − #nodes + 1`, gradients are
`#nodes − 1` — so restricting `ψ` punches a hole in the span, and the missing
gradients would be tree-dependent.

**Falsified.** With `ψ` on every interior node the numbers are identical to
every digit: 3.843e‑02 and 3.533e‑02. Only the conditioning changed, 1.8e+10 →
5.6e+26, which reproduces `RESULTS_COMPARE.md`'s finding that an
insulator-interior `ψ` wrecks the conditioning.

The remaining explanation, and it is a caution about the whole approach:

> **`A_c` lives in the cotree space, and that space is defined by the tree.**
> Jochum's (14) is an *"inexact Helmholtz splitting"* which he realizes
> discretely **by tree–cotree**. So the tree is baked into the trial space
> before the gauge condition ever acts. The gauge then determines `ψ` *given* a
> tree-dependent `A_c`, and the dependence survives into `Φ_eff`.

That would make the gauge condition a statement about how `A` is *split*
between `A_c` and `∇ψ`, not about which `A` the discrete problem selects — which
is precisely the distinction Rapetti draws when he says tree gauges *"are not a
discretization of the Coulomb gauge and enforce no orthogonality."*

**Stated as a limit of this work, not of his paper.** My construction may be
unfaithful: his examples are voltage-driven with no interior terminal (§5.2), so
this configuration is outside what he demonstrates, and the aggregation is mine.
What is measured is that **the construction in §6 does not deliver a
tree-independent terminal**, and that the cause is not the `ψ` restriction.

### 10.4 What this leaves

| | |
|---|---|
| §6's structure | **right** — correct `R`, correct reactance, symmetric, drive in the right row |
| §6 as a fix for the gauge | **wrong** — tree-dependence identical to tree–cotree's, worsening with frequency |
| the `ψ` restriction as the cause | **falsified** |
| likely cause | `A_c`'s space is tree-defined, so the tree precedes the gauge |

If that diagnosis is right, no arrangement of terminals rescues it: the fix
would have to stop realizing the Helmholtz splitting with a tree. That points at
§8's second fallback — taking the condition from an exterior boundary-integral
solution — rather than at any further rearrangement of Jochum's rows.

### 10.5 The control that would settle it was attempted, and failed

The diagnosis in §10.3 needs one check: is **Jochum without a terminal**
tree-independent on this mesh? If his published configuration is also
tree-dependent here, my construction is faithful and the conclusion holds. If it
is not, the terminal is what broke it, and §6 needs revisiting rather than
abandoning.

I built it — his configuration, Dirichlet `V` on faces, no terminal, measuring
the dissipated power, which is gauge-invariant — and **it does not work.**
Two attempts:

| attempt | symptom | cause | outcome |
|---|---|---|---|
| 1, `ψ` on all interior nodes | `max\|V\| = 431 V` for a 1 V drive, `cond = 7.7e+26` | the system is numerically singular at double precision — the same ψ-everywhere conditioning collapse `RESULTS_COMPARE.md` measured at 5.6e+26 | **not the lift** |
| 2, `ψ` restricted per Stysch §5.5 | `max\|V\| = 1.000` ✓, `cond = 1.8e+10` ✓, but `P = 2.4e‑19 W` against an exact `4.35e+04 W` | `E ≈ 0` inside the conductor, where 1 V across 6 mm should give ~167 V/m | **still wrong, different bug** |

Attempt 2 establishes that **the Dirichlet lift was correct all along** — `V`
comes out in the right range once the conditioning is fixed. Something else in
the voltage-driven assembly or in the `E` reconstruction is wrong, and I have
not found it.

The 98–99 % tree-dependence both attempts report is therefore **meaningless**:
it is 98 % of a quantity that is numerically zero.

Worth stating plainly: the **current-driven** assembly in §10.1 is right — it
reproduces the A‑Φ solve digit-for-digit at three frequencies. The bug is
specific to the voltage-driven variant built for this control, which shares the
operator but not the right-hand side or the drive.

**So §10.3 remains a hypothesis**, and the two possibilities it was meant to
separate are both still open:

- my construction is faithful, and Jochum's formulation inherits tree-dependence
  through `A_c`'s tree-defined space; **or**
- my construction is unfaithful in some way the terminal exposes, and his
  published configuration is fine.

What is **not** in doubt, because it was measured three ways and cross-checked
digit-for-digit against the A‑Φ solve: §6 gives the right impedance at every
frequency, and it does not remove the tree-dependence. Whichever explanation is
right, **§6 is not the fix.**

Fixing that control is the cheapest next step, and it is a prerequisite for
trusting §10.3 either way.
