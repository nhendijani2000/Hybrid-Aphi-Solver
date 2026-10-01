# User guide — writing an input file and getting results

How to set up and run a case, and how to get pictures out of it.

This is the **how-to**. `INPUT_FILE_PLAN.md` is the design record for the same
format — it explains why the format is the way it is and what was deliberately
left out, which is what you want when changing the parser and not what you want
when using it.

**Contents**

1. [A complete file, explained](#1-a-complete-file-explained)
2. [How a run works](#2-how-a-run-works)
3. [The sections](#3-the-sections)
4. [**Post-processing — getting pictures**](#4-post-processing--getting-pictures)
5. [Running it](#5-running-it)
6. [When it refuses](#6-when-it-refuses)
7. [Everything at a glance](#7-everything-at-a-glance)

---

## 1. A complete file, explained

A copper ring in air, driven with 1 A at 50 Hz, with one picture asked for.
Everything in it is explained in §3 and §4; read it once to see the shape.

```ini
# Lines starting with # are comments. Keys are case-insensitive;
# names of volumes and surfaces are not -- they must match the mesh.

[mesh]
file        = loop.msh          # relative to THIS file, not your shell
length_unit = mm                # Gmsh files carry no units. You must say.

[analysis]
type        = frequency         # or `dc`
frequencies = 50                # Hz

[Body B1]                       # one section per material region
volume      = ring              # a Physical Volume in the mesh
sigma       = 5.8e7             # S/m -- copper

[Body B2]
volume      = air
sigma       = 0                 # an insulator

[port P1]                       # how the problem is driven
type              = internal_current
surface           = loop_cut    # a Physical Surface cutting the ring
current           = 1.0         # A
current_direction = +y          # which way round the loop

[postprocess PP1]               # optional -- ask for a picture
geometry = body
body     = ring
field    = J
display  = vector
```

Two rules that catch people out:

- **Every volume in the mesh needs a `[Body]`.** Air is a body with `sigma = 0`,
  not something you can leave out.
- **`sigma` is never defaulted.** A conductor cannot silently become an
  insulator because you forgot a line.

---

## 2. How a run works

```
    your .geo  --gmsh-->  mesh.msh  ─┐
                                     ├─>  solve_mesh.exe  ──>  output/
    your .aphi  ─────────────────────┘                          ├─ potential.vtk, {A,B,E,H,J}_field.vtk
                                                                ├─ *.out  (the same, as text)
                                                                ├─ run_summary.txt
                                                                └─ postprocess.json
                                                                       │
                                                 pvbatch tools/postprocess.py output
                                                                       │
                                                                       └─>  PP1_J_vector.png, ...
```

The input file is checked **before** the mesh is opened and long before anything
is solved, so a typo costs a second rather than a minute. Names of volumes and
surfaces are the exception: the parser cannot know them, so a misspelled
`volume = rng` is caught when the mesh is read, a step later.

---

## 3. The sections

Order does not matter. Sections named `[Body B1]`, `[port P1]`,
`[postprocess PP1]` can repeat; the rest appear at most once.

### `[mesh]` — required

| key | required | what |
|---|---|---|
| `file` | yes | path to the `.msh`, **relative to the input file** |
| `length_unit` | yes | `m`, `mm`, `um`, `nm` |

`length_unit` is required because a Gmsh file does not record its units. Getting
it wrong scales your whole problem and the results will look plausible, so it is
worth a second look.

### `[analysis]` — required

| key | required | what |
|---|---|---|
| `type` | yes | `dc` or `frequency` |
| `frequencies` | one of these two | a list in Hz, e.g. `50 1e3 1e6` |
| `f_start`, `f_stop`, `f_points` | one of these two | a sweep; `f_points` counts **both** endpoints |
| `sweep` | no | `log` (default) or `linear` |
| `formulation` | no | `full_wave` (default) or `reduced` |

Give either an explicit `frequencies` list **or** the three sweep keys, not
both. At `type = dc` neither is allowed — there is no frequency to give.

```ini
[analysis]                      [analysis]
type        = frequency         type    = frequency
frequencies = 50 500 5000       f_start = 1e3
                                f_stop  = 1e7
                                f_points = 41        # log-spaced, hits 1e7 exactly
```

### `[Body NAME]` — one per material region, at least one

| key | required | what |
|---|---|---|
| `volume` | yes | one Physical Volume name from the mesh |
| `sigma` | yes | conductivity, S/m. `0` for an insulator |
| `eps_r` | no | relative permittivity, default 1 |
| `mu_r` | no | relative permeability, default 1 |

One volume per section. For two volumes of the same material, write two
sections.

### `[port NAME]` — how the problem is driven

| key | required | what |
|---|---|---|
| `type` | yes | `boundary_current`, `boundary_voltage`, `internal_current`, `internal_voltage` |
| `surface` | yes | one or more Physical Surface names |
| `current` | current types | amplitude in A |
| `voltage` | voltage types | amplitude in V |
| `phase_deg` | no | the drive's phase; `A·cos(ωt + phase_deg)`. Forbidden at DC |
| `current_direction` | internal only | `+x`, `-y`, `z`, … or three numbers |

**`boundary_*` versus `internal_*`** is the distinction that matters:

- A **boundary** port is a terminal on the outside of the domain — the flat end
  of a wire that runs to the edge of the box. Current enters there and leaves
  at another one. You need **two**: one driven, one as the voltage reference.
- An **internal** port is a *cut* through a conductor, used when the conductor
  closes on itself. A ring has no ends to attach to, so the only way to drive it
  is to cut it and inject across the cut. One port is enough: the cut is its own
  reference.

`current_direction` only applies to an internal port, and it is a **hint, not a
normal**. The cut's own faces fix the axis; the hint only resolves which of the
two signs you meant. For a ring in the x-y plane cut at `y = 0, x > 0`, `+y`
reads as "counterclockwise seen from +z".

### `[boundary]` — optional

| key | required | what |
|---|---|---|
| `outer` | no | `flux_tangential` (default and currently the only one) |

`flux_tangential` sets `n × A = 0` on the domain wall. It confines the return
flux, so it behaves like a flux-excluding shell. **This makes the size of your
air box matter**: too small and the inductance reads low and the far field is
wrong. See `ComplexVectorPhasorConcept.md`'s companion case study in
`regression_tests/05_Loop_1A_50Hz` — moving the wall from 2.3 loop radii to 7
changed `L` by 2 %.

### `[solver]` — optional

| key | required | what |
|---|---|---|
| `conditioning` | no | `natural` (default), `row_scaled`, `scaled_phi` |
| `backend` | no | `internal` (default) or `mumps` |

`backend = mumps` needs a build configured with it, and is refused at parse time
rather than silently falling back — a run that quietly used a different
factorization would invalidate any timing you measured with it.

### `[output]` — optional

| key | required | what |
|---|---|---|
| `directory` | yes if the section is present | where results go, relative to the input file |

Relative to the **input file**, not your shell, so a case runs the same from
anywhere.

---

## 4. Post-processing — getting pictures

Add a `[postprocess]` section for each picture you want. Several are fine; they
are rendered in the order written.

Each one answers **three questions**, and that is all the section is:

```
    1. WHERE do I want to look?      geometry  =  body | plane | points
    2. WHICH field?                  field     =  phi | a | e | b | h | j
    3. WHAT about it?                display   =  ... see 4.3
                                   (+ phase_deg, when the answer is "at an instant")
```

**You only write the keys your answers need.** Choosing `geometry = plane` means
you write `plane`, and `body` would be an error — not ignored, an error with a
line number. That is deliberate: a key silently dropped looks exactly like a
solver ignoring you.

### 4.1 Where to look — `geometry`

| `geometry` | what you get | keys it brings |
|---|---|---|
| `body` | the whole of one material region, in 3-D | `body = <volume name>` |
| `plane` | a flat cut through everything | `plane = xy\|yz\|zx`, `offset` |
| `points` | values at coordinates you list → a **line plot**, not a render | `points = x y z  x y z  …` |

```ini
geometry = body            geometry = plane        geometry = points
body     = ring            plane    = xy           points   = 5.7 0 0  6.5 0 0  7.3 0 0
                           offset   = 0.1
```

`offset` is how far the cut sits along its own normal (so `plane = xy` with
`offset = 0.1` cuts at `z = 0.1`). Both `offset` and `points` are in your
file's `length_unit`.

> **Tip.** Slicing exactly at `0` can land on a plane of mesh faces and come out
> ragged. A small offset — a fraction of the smallest feature — avoids it.

### 4.2 Which field — `field`

| `field` | what it is | scalar or vector |
|---|---|---|
| `phi` | the electric scalar potential Φ | **scalar** |
| `a` | the magnetic vector potential **A** | vector |
| `e` | electric field **E**, V/m | vector |
| `b` | magnetic flux density **B**, T | vector |
| `h` | magnetic field **H**, A/m | vector |
| `j` | current density **J**, A/m² | vector |

> **Φ and A are gauge-dependent.** Their values depend on a bookkeeping choice
> inside the solver, not only on your problem — Φ was measured changing by 254×
> between 1 Hz and 50 Hz on an identical mesh. **E, B, H and J are not**; they
> are what you should quote. Φ is still worth plotting to see current paths, and
> the *jump* in Φ across a port is a real voltage; the renderer stamps the
> warning on the figure so a screenshot cannot lose it.

### 4.3 What about it — `display`

**This is the part that needs explaining, so here it is slowly.**

Everything in a frequency-domain run oscillates. The solver does not store a
value; it stores an **amplitude and a phase**, as a complex number `P + jQ`,
meaning the real field is

```
    value(t)  =  P·cos(ωt)  −  Q·sin(ωt)
```

For a **scalar** like Φ that is simple: it swings between `+amplitude` and
`−amplitude`, and "the magnitude" is one number.

For a **vector** it is not simple, and this is the whole reason there are
several options. A vector has three components, each with its own phase. As the
cycle runs, the arrow does not just grow and shrink along a fixed line — **its
tip travels around an ellipse**. An ellipse has a long radius and a short one,
so "how big is E here?" has more than one honest answer:

```
                      ___
                   .-'   `-.          a = the long radius: the biggest the
                 .'    _    `.            arrow ever gets          <- `peak`
                /    .' `.    \
               |    |  +  |    |      b = the short radius: the smallest
                \    `._.'    /
                 `.         .'        the arrow's tip runs round this path
                   `-.___.-'          once per cycle
```

The four magnitude-ish options:

| `display` | what it means | when you want it |
|---|---|---|
| **`complex_magnitude`** | `√(\|P\|²+\|Q\|²)`, which is `√(a²+b²)`, equal to `√2 × RMS` | **the default.** A steady picture that does not flicker, and the one every published figure in this project uses |
| `peak` | `a` — the largest the field ever actually reaches | worst case: insulation, saturation, "does it exceed X anywhere" |
| `magnitude_at_phase` | the length of the arrow **at one instant** `ωt = phase_deg` | snapshots, animations, seeing how a wave moves |
| `axial_ratio` | `b/a`, between 0 and 1 | a check: is `complex_magnitude` safe to read as a peak here? |

And three that are not magnitudes:

| `display` | what it means |
|---|---|
| `vector` | arrows showing the actual field **at the instant** `phase_deg` |
| `phase` | the timing, in degrees. Needs `component` for a vector field (see below) |
| `real`, `imag` | `P` alone, `Q` alone — the raw stored halves |

**Which should you pick?**

```
    Just show me the field.                        ->  complex_magnitude  (the default)
    How big does it get, worst case?               ->  peak
    I want an animation / a snapshot in time.      ->  magnitude_at_phase, or vector
    Which way is the current flowing?              ->  vector
    Is this field leading or lagging?              ->  phase  + component
    Can I trust complex_magnitude as a peak?       ->  axial_ratio
```

**The honest practical note.** `complex_magnitude` and `peak` are the same number
whenever the field is *linearly polarized* — when the arrow swings along a line
rather than circling, so `b = 0`. In ordinary low-frequency work that is nearly
always true: measured across a whole torus case, the axial ratio's median was
`9e-06` and `complex_magnitude` overstated the peak by less than 0.001 %. The
two differ only where the components fall out of phase — near corners, in
strong skin effect, where two current paths interfere. **If you are unsure, plot
`axial_ratio` once.** If it is near zero everywhere, the distinction does not
affect you. If it approaches 1, `complex_magnitude` can overstate the true peak
by up to 41 %.

`ComplexVectorPhasorConcept.md` has the full derivation if you want it.

### 4.4 Which instant — `phase_deg`

Only two displays use it: `magnitude_at_phase` and `vector`. For those,
`phase_deg` picks the moment in the cycle, with the solver's convention
`A·cos(ωt + φ)` — so `phase_deg = 0` is the instant the drive is at its
positive peak, `90` is a quarter cycle later.

**For every other display, `phase_deg` is an error**, because those quantities
do not depend on it and a key that did nothing would be a lie. If you set a
phase and the picture does not move, you want to know, not wonder.

```ini
display   = vector
phase_deg = 90
```

At `type = dc` there is nothing to be a phase of, so `phase_deg` and any
instant-asking display are both refused.

### `component` — only for the phase of a vector

`display = phase` on `e`, `b`, `h`, `j` or `a` needs `component = x`, `y` or `z`.

That is not fussiness. **A vector does not have "a phase."** Each component has
its own, and their disagreement is exactly what makes the ellipse an ellipse.
Picking one for you would be inventing an answer. `field = phi` is a scalar, so
it has only one phase and needs no component.

### 4.5 Worked examples

```ini
# Current flow through a conductor, as arrows, at the drive's peak
[postprocess PP1]
geometry  = body
body      = ring
field     = J
display   = vector
phase_deg = 0

# |B| over a horizontal cut, slightly off centre to miss the mesh plane
[postprocess PP2]
geometry = plane
plane    = xy
offset   = 0.1
field    = B
# display omitted -> complex_magnitude

# E across the conductor's width, as a graph
[postprocess PP3]
geometry = points
points   = 5.7 0 0  6.1 0 0  6.5 0 0  6.9 0 0  7.3 0 0
field    = E

# Is the complex magnitude trustworthy as a peak in this problem?
[postprocess PP4]
geometry = body
body     = ring
field    = B
display  = axial_ratio

# How much does the current lag the drive, deep in the conductor?
[postprocess PP5]
geometry  = plane
plane     = zx
field     = J
display   = phase
component = z
```

---

## 5. Running it

```bash
# 1. mesh it
gmsh loop.geo -3 -o loop.msh

# 2. solve it   -- solve_mesh.exe, NOT aphi_solver.exe
build\solve_mesh.exe loop.aphi

# 3. draw it    -- only if the file has [postprocess] sections
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" tools/postprocess.py output
```

> **`aphi_solver.exe` does not solve.** It reads, validates, binds and
> assembles, then stops — it is a check, useful for confirming a file is sound
> without waiting for a factorization. **`solve_mesh.exe` is the solver.** If you
> ran the wrong one you will get no `output/` and no error, which is confusing
> exactly once.

`solve_mesh.exe` takes two optional extras:

```bash
build\solve_mesh.exe loop.aphi amd              # ordering: natural | rcm | amd
build\solve_mesh.exe loop.aphi --compare-plain  # also solve unequilibrated, to
                                                # compare residuals (2x the time)
```

What lands in `output/`:

| file | what |
|---|---|
| `run_summary.txt` | mesh size, unknowns, fill-in, where the time went |
| `potential.vtk` / `.out` | Φ |
| `{A,B,E,H,J}_field.vtk` / `.out` | one pair per field |
| `postprocess.json` | your requests, for the renderer — only if you asked for pictures |
| `PP1_J_vector.png`, … | the pictures, after step 3 |

Open the `.vtk` files directly in ParaView if you want to explore by hand; the
`[postprocess]` route is for results you want reproducibly, without clicking.

---

## 6. When it refuses

The parser is strict on purpose, and every message carries a line number. The
ones you are most likely to meet:

| message | what happened |
|---|---|
| `unknown key 'x' in [...]` | a typo, or a key that belongs in a different section |
| `'phase_deg' has no meaning for this display` | the display is phase-independent — see §4.4 |
| `'display = phase' on a vector field needs 'component = x, y or z'` | see §4.4 — a vector has no single phase |
| `'display = peak' describes a vector; 'field = phi' is a scalar` | a scalar has no ellipse to have a peak of |
| `[postprocess PP1] is missing the required key 'body'` | `geometry = body` needs to know *which* body |
| `a 'phase_deg' has no meaning when type = dc` | nothing oscillates at DC |
| `'backend = mumps' needs a build configured with …` | that build has no MUMPS; drop the line |
| `body 'B1': volume 'rng' is not a Physical Volume in the mesh. Available: 'ring', 'air'` | caught a step later, when the mesh is read. **It lists the names the mesh does have** — usually the typo is obvious from that line alone |
| `no potential reference: every port is a boundary current source, so Phi is …` | a current-driven boundary port needs a second, voltage-driven port to measure against. An internal cut does not — it is its own reference |

A rule that explains most of the surprising ones: **if a key cannot affect the
answer, writing it is an error.** The parser would rather stop than accept
something you clearly did not mean.

---

## 7. Everything at a glance

| section | repeats | keys |
|---|---|---|
| `[mesh]` | no | `file`, `length_unit` |
| `[analysis]` | no | `type`, `frequencies`, `sweep`, `f_start`, `f_stop`, `f_points`, `formulation` |
| `[boundary]` | no | `outer` |
| `[solver]` | no | `conditioning`, `backend` |
| `[output]` | no | `directory` |
| `[Body NAME]` | yes | `volume`, `sigma`, `eps_r`, `mu_r` |
| `[port NAME]` | yes | `type`, `surface`, `current`, `voltage`, `phase_deg`, `current_direction` |
| `[postprocess NAME]` | yes | `geometry`, `field`, `display`, `body`, `plane`, `offset`, `points`, `phase_deg`, `component` |

**Values**

| key | allowed |
|---|---|
| `length_unit` | `m`, `mm`, `um`, `nm` |
| `type` (analysis) | `dc`, `frequency` |
| `sweep` | `log`, `linear` |
| `formulation` | `full_wave`, `reduced` |
| `outer` | `flux_tangential` |
| `conditioning` | `natural`, `row_scaled`, `scaled_phi` |
| `backend` | `internal`, `mumps` |
| `type` (port) | `boundary_current`, `boundary_voltage`, `internal_current`, `internal_voltage` |
| `current_direction` | `+x`, `-x`, `+y`, `-y`, `+z`, `-z`, or three numbers |
| `geometry` | `body`, `plane`, `points` |
| `field` | `phi`, `a`, `e`, `b`, `h`, `j` |
| `display` | `complex_magnitude`, `magnitude_at_phase`, `peak`, `axial_ratio`, `phase`, `vector`, `real`, `imag` |
| `plane` | `xy`, `yz`, `zx` |
| `component` | `x`, `y`, `z` |

---

## See also

- `INPUT_FILE_PLAN.md` — the design record: why the format is this shape
- `ComplexVectorPhasorConcept.md` — the three magnitudes, derived and measured
- `FIELD_POSTPROCESSING.md` — what is in each output file, array by array
- `examples/` — runnable input files
- `regression_tests/*/` — complete worked cases, each with a validation report
