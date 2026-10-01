# 06 — Two wires in quadrature: a tour of the `[postprocess]` options

Two parallel copper wires driven with equal currents **a quarter cycle apart**.
Thirteen `[postprocess]` requests covering all three geometries, every display,
`phase_deg`, `colormap` and `data`.

The point of the geometry is that it is the only case in the suite where the
field is **elliptically polarized**. `docs/USER_GUIDE.md` §4 explains the
options; this case is where you can see them do something.

## Run it

```bash
gmsh two_wires.geo -3 -o two_wires.msh          # already committed; only if you edit it
..\run_case.bat two_wires.aphi                  # ~8 s
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" ..\..\tools\postprocess.py output
```

Thirteen PNGs land in `output/`. Edit the `[postprocess]` sections in
`two_wires.aphi` and re-run the last line — you do **not** need to re-solve
unless you change the physics.

## What each picture shows

| | shows | look for |
|---|---|---|
| **PP1, PP2** | `J` as arrows in each wire | the current direction |
| **PP3** | `B`, `complex_magnitude` | the phase-independent default — never blinks |
| **PP4** | `B`, `magnitude_at_phase`, 0° | one instant: wire 1 at its peak |
| **PP5** | `B`, `magnitude_at_phase`, 90° | a quarter cycle later: **the bright region has moved to the other wire** |
| **PP6** | `B`, `peak` | the largest `\|B\|` ever reaches — compare with PP3 |
| **PP7** | `B`, **`axial_ratio`** | **the payoff.** Two lobes at `(0, ±5 mm)` reaching 0.98, a null along the line of centres |
| **PP8** | `J`, `phase`, component z | the two wires about 90° apart; inside each, phase varies with radius |
| **PP9** | `E` at ten points | a **line plot**, not a render |
| **PP10** | `B` as arrows on the `zx` cut | a different plane |
| **PP11** | `phi` | the gauge warning the renderer stamps on it |
| **PP12** | PP3 again with `colormap = jet` | the colour maps |
| **PP13** | PP3 again with `data = per_tet` | nodal vs per-tet — PP13 is faceted, PP3 is smooth |

**PP3, PP4 and PP5 are the ones to look at first.** Three pictures of the same
field that look completely different, and nothing is wrong with any of them —
they answer three different questions.

**PP7 against PP3 and PP6 is the second thing.** Where `axial_ratio` reaches 1
the field traces a circle, and there `complex_magnitude` (PP3) exceeds the true
peak (PP6) by √2 = 41 %. Where it is 0 they are equal. That is the whole of
`docs/ComplexVectorPhasorConcept.md` in two pictures.

## Why it is built this way

**Why two wires.** `axial_ratio`, `peak` vs `complex_magnitude`, and `phase` are
all about elliptical polarization. A single conductor has none — `B` is purely
azimuthal, the ellipse collapses to a line, all three go trivial. Case 05's
median axial ratio is `9e-06`.

**Why the quarter cycle.** Two sources 90° apart make `B` rotate. The two-wire
superposition says exactly where, and the verifier asserts both:

```
b/a = 0   on the line of centres — both contributions point the same way
          there, so the sum stays linear however it is phased
b/a = 1   at (0, ±d/2) — perfectly circular
```

**Why 2.5 kHz.** `δ = 1.32 mm` against `a = 2 mm`, so `a/δ = 1.5`. Enough skin
effect that the phase varies with radius *inside* each conductor, which is what
makes PP8 worth looking at.

**Why the mesh is coarse.** About 2 elements per skin depth, where case 04 uses
8. This case is meant to be edited and re-run, so it is sized for an eight-second
solve. **Do not read a skin-depth profile off it** — case 04 is the one for that.

## What it asserts

`../check.bat 06_TwoWires_Quadrature` runs `verify_two_wires.py`:

| | measured | target |
|---|---|---|
| `\|I\|` in each wire | 1.0041, 0.9970 A | 1 A ± 2 % |
| phase difference | 90.09° | 90° ± 1° |
| max `b/a` near `(0, ±d/2)` | 0.975 | ≥ 0.90 |
| `N/a` at the most circular cell | 1.3968 | √2 ± 0.05 |
| median `b/a` over `\|y\| < 0.08d` | 0.0658 | 0.0722 ± 0.02 |

The last one is the analytic value **for that window**, not zero: `b/a` is zero
*on* the axis and grows linearly away from it, so a band median is not a line
value. The first draft of this case asserted zero and was wrong.

## It found a bug on its first run

`axial_ratio` read `2.7e-05` where theory says exactly 1.0. The fault was in
`tools/postprocess.py`, which divided `|P×Q|` by `a` once — giving the
semi-minor axis `b`, not the ratio `b/a`. It needed `a` twice. The wrong
quantity is in the field's own units and looks entirely plausible on a colour
bar, so nothing but a case with a known analytic answer would have caught it.
