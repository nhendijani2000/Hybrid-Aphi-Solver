# Chew against Jochum, same mesh — `compare.py`

```bash
gmsh embedded.geo -3 -o embedded.msh
"C:\Program Files\ParaView 6.1.1\bin\pvpython.exe" compare.py
```

374 nodes, 1330 tets, 1973 edges — a copper cube fully embedded in a dielectric
box, so the conductor/dielectric interface is closed and there are three distinct
node populations: inside the conductor, **on the interface**, and inside the
insulator. Small enough that every condition number is an exact dense SVD.

`tiny.geo` could not be used: its conductor slab reaches the outer boundary and
the mesh is coarse enough that *every* interior node lands on the interface,
which makes the one measurement this script exists for vacuous.

## What it was for

`spike.py` measured Chew alone and reported three problems. `GAUGE_CHOICE.md`
§13 then found Jochum, Farle & Dyczij-Edlinger (SCEE 2014, pp. 63–71), which
claims to solve all three at once in our own discretization. This runs both on
the same mesh and the same element matrices.

It also tests a **falsifiable prediction** derived from their §5.5 one-liner
("set all `ψ` coefficients in the interior of `Ω_N` to zero… (26) still
contributes to unknowns on `Γ`"). Reading the blocks says *why*: for a `ψ` DOF
strictly inside the insulator, as `k₀ → 0`,

```
    (2,1) = jk₀C_κᵀ → jk₀·jk₀C_ε → 0
    (2,2) = jk₀G_κ  → −k₀²G_ε     → 0
    (2,3) = G_κ     → jk₀G_ε      → 0
```

— the **entire row vanishes**. A `ψ` DOF *on* the interface touches conductor
elements where `κ → σ`, which is finite, so its row survives. So the recipe is a
necessity, not an optimization. If the measurement had disagreed, §13's reading
was wrong.

## Result 1 — §5.5's restriction is necessary. Nine orders.

Copper against `ε_r` 4.5:

| freq (Hz) | cond, `ψ` per §5.5 | cond, `ψ` everywhere |
|---|---|---|
| 1e0 | **1.85e+17** | 4.78e+26 |
| 1e4 | **8.73e+17** | 9.59e+24 |
| 1e10 | **1.50e+18** | 1.84e+27 |

**Nine orders at 1 Hz.** And in the `σ = 1` control the predicted mechanism is
visible directly — `ψ`-everywhere degrades monotonically as the frequency falls,
1.19e18 at 10 GHz → 1.31e23 at 1 Hz, while the restricted version does not move.
That is the vanishing row, measured.

## Result 2 — complex-symmetric with losses, to round-off

Symmetry residual **5.7e‑17 … 1.9e‑16** across both sweeps and all ten decades,
*after* the equilibrating congruence.

This is the combination `GENERALIZED_LORENZ_GAUGE.md` §11 concluded was
unavailable — it found the §7 block system was "symmetric or well-conditioned,
not both". Jochum gets both, and not by a scaling applied afterwards: symmetry
comes from **where the test functions go** (Ampère tested by both `w_c` and
`∇ψ̄`, the gauge by `V̄`), so the off-diagonal blocks are transposes by
construction.

## Result 3 — the decisive one: it does not track the material contrast

| | `κ` contrast at 1 Hz | cond at 1 Hz |
|---|---|---|
| Jochum, `σ = 1` (his §6.1 cavity) | 1.80e+10 | 2.67e+17 |
| Jochum, `σ = 5.8e7` (copper) | **2.32e+17** | **1.85e+17** |
| Chew, same mesh, `ε_eff` | 2.32e+17 | **3.32e+22** |

Seven orders more conductivity, and **Jochum's condition number does not
change.** Chew's, on the identical mesh and element matrices, tracks the contrast
across the whole sweep — 3.3e22 at contrast 2.3e17 down to 1.1e9 at 2.3e7.

**That is the spike's finding 2 resolved, structurally.** Not by scaling `ε_eff`
better, but by never forming it: the gauge lives only in the conductor with
weight `σ + jk₀ε`, and the insulator carries physical Gauss with weight `ε_r`.
`σ/ωε` never appears as a ratio inside one operator.

At 1 Hz, Chew 3.32e22 against Jochum 1.85e17 — **five orders**, and Jochum is
symmetric where Chew's §7 form is not with per-element `χ`.

## Two errors in this run, both the same error

Recorded because it is now the **third** time in this project:

**Run 1 reported the condition number *rising* with frequency**, 8.8e6 at 1 Hz to
3.2e23 at 10 GHz — backwards. Cause: I transcribed Jochum's block structure into
our SI variables and dropped his non-dimensionalization. His `A` is a *scaled*
potential and his coefficients are dimensionless, and that scaling is part of the
stabilization. Substituting `A = Ã/c₀` and multiplying by `μ₀c₀ = η₀` reproduces
his (18a) exactly, with `σ₀ = η₀σ`, `ε_r`, `ν_r` and `k₀ = ω/c₀`. Over this sweep
`k₀` spans 2e‑8…2e2 where `ω` spans 6…6e10.

**Run 2's control came out worse than the copper case** — 1e19–1e21 against
1e17 — which is the tell, since the control has seven orders *less* contrast.
Cause: raw diagonal blocks differ by eleven orders. With `h ≈ 1.5e‑3 m`, the
curl–curl block scales as `1/h ≈ 670` while `(2,2) = jk₀G_κ ≈ 1e‑8` at 1 Hz.
Fixed by `equilibrate()`: a **congruence** `D M D` with `D` block-constant,
`d_k = 1/√‖M_kk‖` — Balian et al.'s prescription (`GAUGE_CHOICE.md` §11.1),
scaling equation *and* unknown by the same factor so symmetry is preserved. The
symmetry residual is reported *after* the scaling precisely so that is checked
rather than assumed.

`spike.py` run 1 was the same error (`χ ≈ 1e‑28` put the (2,2) block 1e40 below
the (1,1) and it reported 1e46).

> **Standing rule for this project.** A block system assembled from mixed
> differential orders and mixed material weights is **never** meaningfully
> conditioned as assembled. Equilibrate by congruence before reporting any
> condition number, and run the paper's own materials as a control. Two of the
> five spike runs so far have produced confident, wrong conclusions without
> these.

## What this does not show

- **The absolute level, ~1e17, is unexplained and higher than expected.** Flat,
  symmetric and contrast-insensitive are all solid; the absolute number is not a
  clean win and is not claimed as one. Candidates: the `(3,3)` block `G_epsN` is
  rank-deficient by construction (zero rows for nodes supported entirely in the
  conductor); per-block uniform equilibration cannot fix spread *within* a block;
  or the mesh. Jochum's Fig. 1b is described as "almost constant" without an
  absolute value in the text, so there is nothing to compare against yet.
- **No port, no excitation.** This measures operators, not solutions. The
  mixed-material port — the project's actual blocker — is untouched, and
  `GAUGE_CHOICE.md` §12.3's gap stays open.
- **Chew is measured on the tree–cotree-reduced edge space here**, not the full
  space `spike.py` used, so these numbers are comparable to Jochum's above but
  not to `GENERALIZED_LORENZ_GAUGE.md` §11's.
- `ν_r = 1` throughout. Permeability contrast is untested, which is where
  Demerdash & Wang is still wanted.

## Corroboration from practice

Reported independently by the user from Ansys Maxwell: its A‑Φ matrices were
"very ill conditioned", and **"with coulomb gauge was even more ill
conditioned"**. Three sources now agree that naive gauging *degrades*
conditioning:

| | |
|---|---|
| this run | Chew gauged 3.32e22 against ungauged 2.55e4 at 1 Hz — gauging costs **eighteen orders** |
| Ansari's own Table 1 | gauged 9.0 GB / 1102 s against ungauged 234–379 s, and with conventional preconditioning the iterative solver *"can fail to converge to the correct solution"* |
| Maxwell, in practice | as above |

The lesson is not "do not gauge". It is that **how** the gauge is placed and
scaled decides everything: split by region, weight by the local material, and
apply the balancing as a congruence.
