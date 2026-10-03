# 08_MixedPort_Interior — the mixed-material interior port, measured

```bash
gmsh stub_rod.geo -3 -o stub_rod.msh
solve_mesh.exe mixed_port.aphi              # ~5 min, 188554 unknowns
pvpython measure_gauge.py                   # after the two permuted solves
```

Case 02's copper cylinder with its top cap pulled **inside** the box, carrying a
**single-potential current port** on a face with conductor on one side and air on
the other. The configuration `GAUGE_CHOICE.md` §1.2 argues tree–cotree has no
principled rule for, and the one the project went looking for another gauge to
fix.

Case 03 is the control: same wire, same material, same 1 A, both caps on the
outer boundary.

## The port works, and the physics is right

Before any gauge claim, the case has to be sound — the previous attempt at this
measurement was abandoned because its geometry had no return path and the
conductor-only control was equally wrecked.

| check | result |
|---|---|
| Φ spread across the port face, 373 nodes | **exactly 0.0** — one shared DOF, as `dof_map.cpp` specifies for `PhiDof::Port` |
| port direction `d` | `(0, 0, −1)` — into the conductor |
| backward error, three solves | 2.95e‑21, 2.94e‑21, 2.87e‑21 |
| **R** | **68.659 µΩ** against **68.58 µΩ** expected (case 03's 40 mm rod gives 97.97 µΩ; this one is 28 mm) — **0.12 %** |

Resistance scales with length as it must, so the interior mixed port is
physically correct and not merely well-posed.

## The measurement

Three spanning trees, from permuting the node order of the same mesh (97.9 % and
97.7 % of nodes moved).

### Terminal quantities

| | R (µΩ) | L (nH) | ΔR | ΔL |
|---|---|---|---|---|
| base | 68.659272 | 20.382779 | — | — |
| permA | 68.660045 | **19.747798** | 1.13e‑05 | **3.1e‑02** |
| permB | 68.606254 | 20.435042 | 7.72e‑04 | 2.6e‑03 |

**L moves by 3.1 % with the choice of spanning tree.** `R` barely moves at all.

That split is the `ωL/R` mechanism of `06_TwoWires_Quadrature` §7 behaving
exactly as derived: here `ωL/R = 0.093`, so the gauge-sensitive fraction of the
terminal voltage is small and almost entirely **inductive** — and it is the
inductance that moves while the resistance does not.

### Fields, which are gauge invariant by theorem

| field | base→permA | base→permB |
|---|---|---|
| E | 1.45e‑01 | 1.20e‑01 |
| B | 4.95e‑01 | 7.32e‑01 |
| H | 4.95e‑01 | 7.32e‑01 |
| **J** | **7.72e‑01** | 6.37e‑01 |

`E`, `B`, `H` and `J` **cannot** depend on the gauge. Case 07 measures them
invariant to 1e‑12 with both ports on the outer boundary. Here they move by tens
of percent, which means **these are not two gauges of one solution — they are two
different discrete solutions.**

### It is localized at the port

`J`, by distance from the port face, against the band's own maximum as well as
the global one (the global max sits at the port, so normalizing by it alone
would flatter the far field):

| band | nodes | vs global max | vs this band's max |
|---|---|---|---|
| 0–2 mm | 10559 | **7.72e‑01** | **7.72e‑01** |
| 2–5 mm | 11282 | 1.76e‑03 | 4.71e‑03 |
| 5–10 mm | 18980 | 1.09e‑04 | 2.92e‑04 |
| 10 mm+ | 69215 | 5.31e‑05 | 1.42e‑04 |

**77 % inside 2 mm, 0.47 % by 5 mm, 0.014 % far away** — a factor of 160 across
one band boundary. Restricting to nodes inside the conductor gives identical
numbers, so this is not an artefact of reporting `J` where `σ = 0`.

That 69215 far-field nodes agree to 1.4e‑04 is also the **internal control on the
measurement itself**: positions are matched exactly, and a broken comparison
would be wrong everywhere rather than within 2 mm of one face.

## Where this sits against the two existing baselines

| port configuration | local `J` | terminal `L` |
|---|---|---|
| both ports on the outer boundary (case 07) | 1e‑12 | 1e‑15 |
| cut inside a conductor (case 05, §1.3) | 8.3 % | 0.018 % |
| **interior conductor/air port (this case)** | **77 %** | **3.1 %** |

**About ten times worse in the field and about 170 times worse in the
inductance** than an interior cut through a conductor — which was itself already
outside what case 07's result covers.

3.1 % on an inductance is a number that fails an engineering tolerance, and it is
set by an arbitrary graph traversal.

## What this does and does not establish

**Establishes:** the mixed-material interior port is expressible, well-posed and
gives the right DC resistance, *and* its inductance and near fields are
controlled by the spanning tree rather than by the physics. The failure the
project has been reasoning about from the literature is now measured in our own
solver, with a quantitative baseline a new gauge has to beat.

**Does not establish a mechanism.** The most likely one is the field singularity
at the port rim, where the cap meets the cylinder's lateral surface: the solver
already warns, for internal cuts, that a rim interior to Φ's support makes the
field singular there and the gap capacitance mesh-dependent. A singular field has
no determinate discrete value, and the tree selects among the possibilities. That
is consistent with the sharp localization but is **not proven here**.

## The control, and it is decisive

```bash
pvpython measure_gauge.py control
```

The same wire, the same 1 A drive, the **same permutation seeds**, and the same
measurement script — with both caps on the outer boundary. That is case 03, run
three times through this harness. If the 77 % were an artefact of the script,
the matching, or the stub-rod mesh, it would appear here too.

| | R (µΩ) | L (nH) | ΔR | ΔL |
|---|---|---|---|---|
| base | 97.969950 | 22.595612 | — | — |
| permA | 97.969950 | 22.595612 | 4.2e‑12 | 9.9e‑12 |
| permB | 97.969950 | 22.595612 | 6.5e‑12 | 1.3e‑11 |

Identical to ten significant figures, and `Z = 9.796995e‑05 + j7.0986e‑06`
reproduces case 02's published impedance exactly — so the harness is also
measuring the right physics, not merely measuring it consistently.

| field | base→permA | base→permB |
|---|---|---|
| E | 9.7e‑12 | 1.0e‑11 |
| B | 3.9e‑11 | 3.6e‑11 |
| H | 3.9e‑11 | 3.6e‑11 |
| J | 1.5e‑11 | 1.7e‑11 |

And `J` by band: 7e‑12, 1.3e‑11, 1.5e‑11, 1.5e‑11 — **flat**, with none of the
localization the mixed port shows.

**The gauge really did change in the control.** `Φ` moves 0.79–0.95 % and `A`
moves 66–100 %, comparable to the mixed-port run. So the three spanning trees are
genuinely different in both cases; the only difference is whether the
**observables** notice.

### The two side by side

| | boundary ports (control) | interior mixed port | ratio |
|---|---|---|---|
| terminal `L` | 9.9e‑12 | **3.1e‑02** | **3e9×** |
| `J` | 1.5e‑11 | **7.7e‑01** | **5e10×** |
| `Φ` moves (the gauge changing) | 7.9e‑03 | 1.8e‑02 | ~2× |
| `A` moves | 1.0e+00 | 6.4e‑01 | ~1× |

The gauge transformation is the same size in both. The response to it differs by
ten orders of magnitude, and the only thing changed is **where the port sits**.

That closes the alternative explanations. The 77 % is the solver's answer.

**Nothing here is asserted as a regression yet** beyond the invariants:
`expected.txt` pins the single-valuedness of the port Φ and the DC resistance,
which must hold under any gauge. The gauge-movement numbers are *recorded*, not
asserted — asserting them would be asserting the defect.

## Practitioner corroboration

Reported independently from Ansys Maxwell: a mixed-material port there ran to
completion without crashing and returned wrong solutions, and its A‑Φ matrices
were badly conditioned, more so under a Coulomb gauge. This case's factorization
reports `d max/min = 7.05e+12`, identically across all three trees.
