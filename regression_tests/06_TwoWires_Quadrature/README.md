# 06 — Two wires in quadrature: a tour of the `[postprocess]` options

Two parallel copper wires driven with equal currents **a quarter cycle apart**.
Sixteen `[postprocess]` requests covering all three geometries, every display,
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

Sixteen PNGs land in `output/`. Edit the `[postprocess]` sections in
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
| **PP14, PP15, PP16** | the PP3/4/5 comparison again, for `J` | the current is what is *driven*, so this is the most direct view of what a port `phase_deg` does |

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

**The mesh resolves the skin depth**, at about 4.2 elements per `δ` — median
element 0.317 mm against `δ = 1.32 mm`. 192k tets, ~20 s.

The refinement is spent **inside** the conductors, which needs a little care: a
distance-to-the-surface field cannot express it, because distance is positive on
both sides, so asking for 0.25 mm at the surface also refines a 2 mm shell of
*air*. That came to 660k elements, more than half of them air, which varies on
no such scale. Two gmsh `Cylinder` fields instead.

### What the skin effect actually comes out as

`pvbatch skin_profile.py output` — azimuthal mean of `|J|` against Kelvin's
`J₀(kr)/J₀(ka)` at `a/δ = 1.51`:

| `r/a` | measured | Kelvin | err | lag | Kelvin | err |
|---|---|---|---|---|---|---|
| 0.00–0.20 | 0.8035 | 0.8003 | **+0.39 %** | −53.14° | −53.20° | **+0.06°** |
| 0.20–0.40 | 0.8047 | 0.8027 | +0.26 % | −47.90° | −48.07° | +0.16° |
| 0.40–0.60 | 0.8183 | 0.8173 | +0.12 % | −37.56° | −37.75° | +0.19° |
| 0.60–0.80 | 0.8629 | 0.8629 | −0.00 % | −22.65° | −22.89° | +0.24° |
| 0.80–0.92 | 0.9345 | 0.9357 | −0.12 % | −8.70° | −8.93° | +0.23° |

A 53° lag from surface to core, reproduced to **0.06°**, and the magnitude
profile to better than 0.4 % everywhere.

> **BOTH SIDES MUST BE ANCHORED AT THE SAME RADIUS**, and an earlier version of
> this table did not. Kelvin's profile is normalised at `r = a`; a cell band "at
> the surface" is `r > 0.92a`, whose **mean radius is 0.9566a**. Dividing the
> measurement by `J(0.9566a)` while dividing the theory by `J(a)` inflates every
> measured ratio by `1/0.9660 = +3.5 %` — an offset with **no `h` in it**.
>
> It told a convincing false story: errors of 7.2 % that "improved" to 3.9 %
> under refinement and then stalled, with the apparent convergence order
> collapsing to 0.37. All of that was a constant offset dominating a real error
> ten times smaller. What finally placed it was two controls — a **single
> isolated wire** showed the identical error, ruling out the neighbour, and
> narrowing the axial sample band changed nothing, ruling out end effects. An
> error that survives removing the neighbour, ignores position and ignores mesh
> size is in the instrument, not the model.

### Does the refinement earn its cost?

Measured properly, at three levels:

| median `h` | per `δ` | wire tets | core `\|J\|` err | core lag err |
|---|---|---|---|---|
| 0.724 mm | 1.8 | 1,604 | +1.483 % | +0.304° |
| 0.502 mm | 2.6 | 4,977 | +0.661 % | +0.265° |
| 0.317 mm | 4.2 | 19,915 | **+0.390 %** | **+0.060°** |

Yes — 1.48 % to 0.39 % for 8 s to 20 s. The pairwise order in `h` comes out
between 1 and 2 for the magnitude and is too noisy to read for the lag, whose
errors are already down at a tenth of a degree. **Three points are not enough to
claim an order**, and this report does not claim one.

### And the proximity effect, which is not an error

| | `\|J\|` facing the other wire / away |
|---|---|
| wire 1 | **1.137** |
| wire 2 | **0.766** |

The neighbour's field redistributes current around the circumference — ±14 % and
−23 %. It is what makes the arrows in PP1 asymmetric. The two wires are pushed
**opposite** ways because they are 90° apart: at any instant one is near its peak
while the other is near its zero crossing.

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

## A note on PP8, and on phase plots generally

`display = phase` needs masking wherever the field is zero, and this case is
what showed it. In the air `sigma = 0`, so `J` is exactly zero and has **no
phase** — a zero vector points nowhere. But `atan2` returns an angle regardless,
and which angle depends on the **sign bit of the zero**: of 37,706 air nodes,
21,727 hold *negative* zero, and `atan2(+0.0, -0.0)` is `pi`. So the first
version of PP8 painted half the frame a confident, uniform **+180 degrees**,
indistinguishable from a measurement.

The renderer now writes NaN where the component's magnitude is below `1e-9` of
its maximum, and draws NaN in white — a colour the scale cannot produce. Drawing
it at the *bottom* of the map instead would read as a real value at the low end
of the range, which is the same mistake in a nicer colour. **Absent data must
never be drawn in a colour the data itself can take.**

Nothing is wrong in the solver: a scan of all six output files finds no NaN and
no Inf anywhere, and `J` in air is a clean `0.0` exactly as `J = sigma E`
requires. The NaN is a display instruction, added by the renderer.

## It found a bug on its first run

`axial_ratio` read `2.7e-05` where theory says exactly 1.0. The fault was in
`tools/postprocess.py`, which divided `|P×Q|` by `a` once — giving the
semi-minor axis `b`, not the ratio `b/a`. It needed `a` twice. The wrong
quantity is in the field's own units and looks entirely plausible on a colour
bar, so nothing but a case with a known analytic answer would have caught it.
