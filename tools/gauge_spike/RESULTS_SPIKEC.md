# Spike C — can any gauge protect an interior terminal?

```bash
gmsh stub.geo -3 -o stub.msh
"C:\Program Files\ParaView 6.1.1\bin\pvpython.exe" spikeC.py
```

A copper stub standing on the domain floor: bottom face on the outer boundary as
the Φ = 0 reference, top face interior with conductor below and air above — the
mixed terminal. 382 nodes, 1333 tets, dense algebra, three trees from shuffling
the edge scan.

`08_MixedPort_Interior` measures the real solver losing this port (L moves 3.1 %,
J 77 %) and `GAUGE_CHOICE.md` §15.9 derives that **no tree** can fix it. So the
gauge has to change. This asks which one.

## The acceptance test

A candidate must do all three:

1. reproduce the exact DC resistance, `R = L/(σA) = 11.494253 µΩ`;
2. produce a **sensible reactance** (the tree–cotree solve gives L ≈ 3.3 nH);
3. hold the terminal Φ **invariant across trees**.

**Nothing has met all three.** What follows is why, which turned out to be worth
more than a pass would have been.

## Milestone 1 — the test bed, validated

| | |
|---|---|
| R | **11.5007 µΩ** against exact 11.4943 — **0.06 %** |
| terminal Φ across trees | **3.84 %**, 3.53 % |
| symmetry residual | 1.6e‑19 |

The real solver moves L by 3.1 %; the miniature moves the terminal Φ by 3.8 %.
Same order, and the resistance is right to 0.06 %, so the bed is **validated**
rather than merely self-consistent.

**Two bugs the magnitude check caught**, each of which would have been a wasted
week in C++:

- the continuity row is `div[β(jωA + ∇Φ)] = 0`, so it carries a `jω` that the
  first version dropped. Result: a terminal potential of `−4.3e8 − 1.5e9j` V for
  1 A into a 6 µΩ stub — fourteen orders out, **while still showing a plausible
  66 % tree dependence**. Dividing the row by `jω` also restores symmetry, the
  same device Zhao & Fu apply as `−j/ω` (§12.1);
- the first geometry was `embedded.geo`, whose conductor is fully enclosed, so
  injected current could only leave through `σ = 0` dielectric. No return path,
  ill-posed — the same flaw that killed the first attempt at
  `08_MixedPort_Interior`.

## Milestone 2 — Jochum has nowhere to put a port

Bolting a current terminal onto his third row gives a terminal potential of
`−1.149425e‑05` V that is identical across three trees to **2.4e‑14**, has
**zero reactance**, and **does not change from 50 Hz to 50 GHz**. A terminal
impedance constant over nine decades is not physics.

The cause is structural. **Inside the conductor Jochum's third equation is the
gauge**, `div[κ(A_c + ∇ψ)] = 0`, which contains no V at all; V enters the third
equation only in the insulator, as Gauss's law. Measured on the terminal row:

| | |
|---|---|
| gauge terms, `G_kap` | 7.5e+07 |
| the only place V enters, `G_epsN` | 4.2e‑02 |
| ratio | **5.7e‑10** |

So the drive perturbs the gauge instead of injecting current, and exact R with no
reactance is the signature of a decoupled DC conduction problem — whose solution
is linear in z and therefore exact in P1.

His paper has no ports: §13.5's examples are a cavity and an RLC field model,
both driven by boundary conditions. And his gauge interface condition (32) is
*derived* by requiring the boundary integrals of (23) and (27) to cancel — a
terminal on Γ changes exactly those integrals, so **(32) would have to be
re-derived with a port present.** That is a derivation, not a spike.

## Milestone 3 — Chew has a place for a port, but α decides the answer

His (2), `div(ε∇Φ) + χω²Φ = −ρ`, contains **no A**, so a terminal is an ordinary
boundary condition on that equation rather than a perturbation of a gauge row —
which is exactly what Jochum lacks. And no tree enters it, so:

| | |
|---|---|
| tree-dependence | **0.000e+00 — exactly zero**, symmetry residual exactly 0 |
| R at α = 1 | 1.094475e‑05 Ω against exact 1.149425e‑05 — **4.8 % low** |

Invariance by construction, not by conditioning. But the answer moves with
**α, a free gauge parameter**:

| α | \|Φ\|/I (Ω) | vs exact |
|---|---|---|
| **1** — Chew's recommendation | 1.094475e‑05 | **4.8 % low** |
| 1e‑3 | **1.149425e‑05** | exact |
| 1e3 | 2.595849e‑07 | 97.6 % wrong |

A gauge parameter must not change a physical potential. Cause: `χ = αμ ε_eff²`,
and with `ε_eff ≈ −jσ/ω` that gives `χω² = −αμσ² ≈ −4.2e9` at α = 1 — a huge
spurious screening term in what should be a conduction equation. **Chew's paper
has no finite conductivity** (PEC conductors, `ε` varying by 4.5), so `χ` was
never meant to carry `σ`.

At α → 0 the term vanishes and the equation becomes `div(ε_eff ∇Φ) = 0`, which is
also Ostrowski & Hiptmair's EQS-as-gauge (§11.3) and Stysch's MQS Φ BVP (§11.8) —
three sources converging. Measured there: **R = 11.494253 µΩ, all eight digits,
at every frequency from 50 Hz to 50 GHz, tree-dependence exactly 0 — and X ≈ 0.**
No magnetic coupling in the equation, so no inductance.

## Milestone 4 — Φ as a post-process of E: R exact and invariant, no reactance

Stysch §3.1/§6.3 drives Φ's own BVP with the field:
`⟨ε∇Φ, ∇λ⟩ = −⟨εE, ∇λ⟩`. Since **E is gauge-invariant**, Φ should inherit that.

| | Φ from its own BVP | Φ read from the gauged A‑Φ system |
|---|---|---|
| **R** | **11.494253 µΩ — exact to 8 digits** | 11.5007 µΩ |
| **tree-dependence** | **1.1e‑14** | **3.8e‑02** |
| **X** | **−6.5e‑20 ≈ 0** | 1.043e‑06 → L ≈ 3.3 nH |

**This is a real result: the resistance can be made exactly correct and exactly
tree-independent, and it is cheap.** Step 1 is the ordinary tree–cotree solve;
step 2 is one extra scalar Poisson solve.

But the reactance is structurally absent. `div(εE)` keeps only the
**irrotational** part of E; the inductive `−jωA` is essentially solenoidal and
the divergence discards it. Φ recovers the **galvanic** part exactly and loses
the **inductive** part — Ansari's galvanic/inductive split, measured.

## What the four milestones together say

Stysch §2, already recorded in §11.8 and not previously connected to this:
`V = V_c + jωI·L_r` with `V_c = Φ(T_b) − Φ(T_a)`, and *"as it incorporates a
partial inductance, V_c formally depends on the gauge condition."*

**Partial inductance is a gauge-dependent modelling convention.** `L_loop =
L_c + L_r` is invariant; the split between conductor and return path is not. So:

| | |
|---|---|
| **R** | gauge-independent, and computable exactly — milestone 4 |
| **L at a terminal** | needs a convention for attributing loop inductance. **The gauge *is* that convention.** |

A boundary port gets the convention free: `n×A = 0` pins ψ, which is why case
03's L is invariant to 1e‑12. An interior port has nothing supplying it. That is
the 3.1 %, and it is why §15.9's derivation holds.

So the remedy must supply a **physical condition at the port face** that fixes
the split. Jochum's (32) is the only thing in the literature with the right
shape — materially weighted, at the conductor/insulator interface — and the BI
formulation would supply it from the exterior solution. Neither is available as
published.

**And the common root of all four failures: neither formulation was written for
a lossy conductor carrying a terminal.** Chew has PEC conductors and no `σ`;
Jochum has `σ` but no ports.

## Limits of this spike

- 382 nodes, P1 Φ (the solver uses P2), one geometry, `ν_r = 1` so no
  permeability contrast.
- The trees come from shuffling the edge scan, not from permuting the mesh file
  as `07_GaugeInvariance` does. Same lever, not the same code path.
- Milestone 4's step 2 reconstructs `E` at tet centroids from Whitney‑1 `A` and
  P1 `∇Φ`, which is first-order accurate; the 8-digit agreement on R is
  therefore better than the reconstruction deserves and probably reflects the
  exactness of a linear-in-z solution in P1 rather than general accuracy.
- Nothing here is in C++, and nothing should be until a candidate clears all
  three acceptance criteria.
