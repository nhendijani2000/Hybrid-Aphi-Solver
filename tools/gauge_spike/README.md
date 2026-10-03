# Gauge spike

A throwaway investigation, kept because of what it found.

```bash
gmsh tiny.geo -3 -o tiny.msh
"C:\Program Files\ParaView 6.1.1\bin\pvpython.exe" spike.py
```

Seconds to run. 164 nodes, 518 tets — small enough that condition numbers come
from a dense SVD and are **exact**, not estimated.

## What it was for

Before committing to implement `docs/GENERALIZED_LORENZ_GAUGE.md`, three
questions were worth a day rather than a month:

1. Is the §7 block system exactly complex-symmetric, so MUMPS's `LDLᵀ` applies?
2. Does the gauge condition as well as Chew reports?
3. Does any of that survive finite conductivity?

**Yes, yes, and no.**

## What it found

### The formulation works, and spectacularly, for dielectrics

Chew's own material setup — dielectric everywhere, no `σ`:

```
freq      | cond gauged  | cond ungauged
1e+00     | 2.492e+02    | 1.884e+18
1e+04     | 2.492e+02    | 2.584e+14
1e+10     | 3.907e+02    | 5.729e+02
```

**Flat across ten decades** while the ungauged system degrades by sixteen orders
as `ω → 0`. That is Chew's Fig. 5 reproduced in character, and it is what
validates the assembly.

### It does not survive `σ` folded into `ε`

```
freq      | cond gauged  | cond ungauged | eps contrast
1e+00     | 2.046e+22    | 1.655e+06    | 2.317e+17
1e+04     | 9.272e+17    | 5.745e+02    | 2.317e+13
1e+10     | 1.527e+10    | 1.118e+09    | 2.317e+07
```

The condition number **tracks the ε contrast**, at roughly `contrast × 1e4`, and
below about 1 GHz **the gauged system is worse than the ungauged one.**

The cause is not a bug. `ε_eff = ε − jσ/ω` makes `σ/ωε = 2.3e17` at 1 Hz between
copper and dielectric. **Chew's paper contains no finite conductivity** — his
conductors are PEC boundary conditions — so his `ε` varies by 4.5 and never by
1e17. `ε_eff` is our addition, and it is what breaks the conditioning.

### `χ` must be position-dependent, and that costs the symmetry

| | |
|---|---|
| `χ = μ₀ε₀²` constant, his §II.C | 7.0e37 at 1 Hz |
| `χ = με²` per element, his §II.A | 2.0e22 at 1 Hz |

About **fifteen orders**, because `ε²/χ = 1/μ = ν` self-normalizes to the
curl–curl coefficient only when `χ` carries the local `ε`.

But §7's symmetrization is a *global row scaling by χ*, which is a congruence
only when `χ` is constant. Measured with `χ` per element, `K_NE` is no longer a
scalar multiple of `K_ENᵀ` — the residual is **O(1)**, not 1e‑16. So as derived,
the block system is **symmetric or well-conditioned, not both**.

## Why the controls are the point

The first two runs concluded the formulation was unusable, at 1e38. Both were
wrong:

- run 1 — my own block scaling, not the formulation: `χ = μ₀ε₀² ≈ 1e‑28` put the
  (2,2) block 1e40 below the (1,1) block. The scaling `c` is free (it cancels in
  the Schur complement), so it can be chosen for conditioning;
- run 2 — finite conductivity, which is not in the paper at all.

Only running **Chew's materials unchanged** showed the assembly had been correct
the whole time. Without that control the spike would have reported a sound
method as broken.

## Limits

- The **surface term of (24)** is omitted. It must vanish under our boundary
  conditions or be carried; carrying it may break symmetry independently of `χ`.
- 518 tets. The *trends* across ten decades are the result; absolute condition
  numbers will differ on a real mesh. Chew reports 5.94e5 on his.
- `n × A = 0` is imposed by dropping boundary edges. No ports, no excitation —
  this measures operators, not solutions.
- Dense throughout. It cannot be pointed at a real mesh without being rewritten
  sparse.
