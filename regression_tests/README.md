# Regression tests

Each folder is one self-contained case: its geometry, its input files, and its
own `output/`. Nothing here writes outside its own folder.

```
regression_tests/
  run_case.bat            the runner
  check.bat               solve every case and assert its physics
  verify.py               the assertions, run by check.bat
  01_OneCylinder/
    cylinder.geo          geometry
    cylinder.msh          generated from the .geo
    cylinder_1hz.aphi     1 Hz
    cylinder_50hz.aphi    50 Hz -- the one check.bat runs
    cylinder_sweep.aphi   50 Hz to 500 Hz
    expected.txt          the known-good physics, and which .aphi to run
    output/               everything the solver writes
```


## Checking the physics before and after a change

`run-tests.bat` is unit tests and finishes in seconds. It does **not** solve
anything, so it cannot tell you that a change moved `R`, `L`, `Phi` or `J`.
That is what this is for:

```
regression_tests\check.bat                  solve every case and assert its physics
regression_tests\check.bat all verify-only  re-check existing output, seconds
regression_tests\check.bat 03_Cylinder_1A_50Hz    just one
```

It runs every case folder containing an `expected.txt`, and exits non-zero
naming the ones that drifted. About 6 minutes per case with the internal
solver; seconds in `verify-only`.

**Run all three, not one.** Each covers something the others cannot.

| case | drive | `a/delta` | skin effect in `R` | what only this case can catch |
|---|---|---|---|---|
| `01_OneCylinder` | 1 V | 1.07 | 2.6 % | a mild skin effect against Kelvin; the **gauge study** (1 Hz vs 50 Hz) |
| `02_Ansys_Cylinder_50Hz` | 1 V | 0.16 | ~0 | `R` against the exact DC value; `Phi` as a clean `z/L` gradient |
| `03_Cylinder_1A_50Hz` | **1 A** | 0.16 | ~0 | port current self-consistency; the same mesh as 02 |
| `04_Cylinder_SkinDepth` | 1 V | **3.00** | **77 %** | the skin effect where it **dominates** — `\|J\|` falls 3.7x to the axis |

02 and 03 are **duals** — 02 drives 1 V and reads the current out, 03 drives
1 A on the same mesh and reads the voltage out — and they must report the same
`R` and `L`, because those belong to the geometry and the material rather than
to how the thing is driven. **A change that breaks that duality shows up as the
two disagreeing, which neither case alone can see.**

01 and 04 sit in the other regime, and that is why they are worth the runtime:
02 and 03 are at `a/delta = 0.16`, where the profile is flat to five digits and
a bug in the skin-effect physics would be invisible.

**01 and 04 are not redundant with each other.** At `a/delta = 1.07` the skin
effect is only 2.6 % of `R`, so matching Kelvin to 0.038 % validates *the effect
itself* to about 1.4 %. At `a/delta = 3.00` it is **77 %** of `R`, so matching
to 0.22 % validates the effect to 0.29 %. 04 is the stronger statement; 01 is
cheap, carries the gauge study, and covers the weak-effect end where a
formulation could get the limit wrong.

**A case may omit a check, and the omission is a claim.** `verify.py` skips any
check whose key is absent from `expected.txt` and prints `SKIPPED` with the
reason. Case 01 omits both `Phi` checks because at `omega*L/R = 2.89` the
tree-cotree gauge dominates `Phi` and it is *not* a `z/L` gradient — asserting
one there would be asserting something false. It omits `R_rel_error_max` for
the same kind of reason: `R` is 2.6 % above `R_dc` because of the skin effect,
so it is checked against Kelvin's `R_ac/R_dc` instead. Each omission is
documented in the `expected.txt` that makes it.

**`expected.txt` may also name the input to run,** with an `input` line. Case 01
needs this because its folder holds three `.aphi` files — 1 Hz, 50 Hz and a
4-frequency sweep — and without it `check.bat` would pick the sweep, which
writes a different frequency into the output directory `verify.py` reads.

The suggested sequence around any change to the solver, the assembly or the
post-processing:

```
build.bat  &&  run-tests.bat  &&  regression_tests\check.bat
```

## Running one

From inside the case folder:

```
cd regression_tests\01_OneCylinder
..\run_case.bat cylinder_50hz.aphi
```

That builds the solver if needed, meshes the `.geo` if `cylinder.msh` is
missing or older than it, runs the case, and lists what it produced. Set
`SKIP_BUILD=1` to skip the build step when iterating.

Outputs land in `output/` because the input file says so:

```
[output]
directory   = output
```

That path is resolved **relative to the input file, not the working
directory**, so a case writes to the same place however it is invoked. Without
that, running from the repository root and from inside the case folder would
scatter results in two places, and the second run would look as though it had
produced nothing.

### What you need installed

| | why | checked how |
|---|---|---|
| Visual Studio (any recent) | to compile the solver | `build.bat` locates `vcvars64.bat` itself — no developer prompt needed |
| gmsh | to turn `cylinder.geo` into `cylinder.msh` | `run_case.bat` finds the winget install automatically, or put `gmsh` on PATH |
| ParaView | figures (`pvbatch`) | only needed for the plots, not for the solve |
| Git Bash | `build_doc.sh` | only needed for the HTML/PDF documents |

### After the solve: figures, then documents

`run_case.bat` stops once the fields are written. Two more steps produce the
pictures and the report:

```
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" make_plots.py
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" pv_bcell.py
sh build_doc.sh
```

`make_plots.py` writes the report figures, `pv_bcell.py` the per-cell-vs-nodal
`B` study, and `build_doc.sh` rebuilds `CopperRodValidation.html` and `.pdf`
with the new figures inlined. **Run `build_doc.sh` after regenerating figures**
or the documents keep showing the previous run's pictures — the images are
baked into them.

The PDF step launches a headless browser. Some sandboxed shells block that; the
HTML is still written and the script says it skipped the PDF.

### Cost, and how to tell it worked

By default each frequency is factorized **once**, equilibrated. Pass
`--compare-plain` to also solve it unequilibrated and print the two residuals
side by side -- that is a second full factorization and roughly doubles the
runtime, since factorization is ~97 % of a large run. The answer written to
disk is the equilibrated one either way.

Every run writes `output/run_summary.txt` — the same timing table the console
shows, plus the mesh, unknown count and fill-in. That file is the record of
what a run cost; the per-stage timings are also stamped into the header of
every `.out`. Nothing else has to be captured from the console.


`02_Ansys_Cylinder_50Hz` is the heavier of the two cases:

    mesh          18994 nodes, 101873 tets, a few seconds
    solve         240990 unknowns, about 8 min, 4.03 GB peak
    figures       about 2 min

The solve prints a line per frequency. These are the numbers to compare:

    unknowns      240990
    backward err  4.8e-21
    gauge         complete

Then `R` from `tools/pv_extract_rl.py` should read `9.796995e-05` ohm against
`9.796863e-05` exact. If those match, the run reproduced.

### Forcing a re-mesh

`run_case.bat` only re-meshes when `cylinder.msh` is missing or older than
`cylinder.geo`. Editing the `.geo` is therefore enough; to force it without an
edit, delete `cylinder.msh`.

## 02_Ansys_Cylinder_50Hz

Full write-up with figures: [02_Ansys_Cylinder_50Hz/REPORT.md](02_Ansys_Cylinder_50Hz/REPORT.md)

The same solver on a SLENDER conductor: `a = 1.5 mm`, `L = 40 mm`, box 40 mm, a
40-gon, at the same 50 Hz. Built to match the regime of the Ansys Maxwell A-Phi
voltage example, whose dimensions are inferred in the .geo header.

`a/delta = 0.1605` and `omega*L/R = 0.0725`, so this case is resistance dominated
and `Phi` is readable -- the opposite branch from 01. Measured:

| | |
|---|---|
| `Phi` vs the exact `z/L` | worst 8.26e-04, mean 2.51e-04, over all 95149 wire nodes |
| `Phi` at mid height | 84 nodes within 1 um, mean 0.499996, spread 6.02e-04 |
| the same spread on case 01 | **0.31** |
| `R` | 9.796995e-05 ohm against 9.796863e-05 DC exact, **1.3e-05** |
| `L` from `Im(Z)/omega` | 22.60 nH |
| `J(0)/J(a)` | 0.999968 against the Bessel 0.999959 |
| `B` vs exact, inside the conductor | every band within **0.5 %** |
| `B` azimuthal scatter at the peak | 0.0138, coherent content 0.31 % at `m = N` |
| nodes / unknowns | 18994 / 240990, 503 s factorise, 4.03 GB |

**`a/delta` and `omega*L/R` are the same parameter**, both scaling as
`omega a^2`, so a case cannot show a strong skin effect and a readable potential
at once. 01 takes one branch, 02 the other. Together they show the solver
tracking the Bessel profile across both: 1.062 at `a/delta = 1.07`, 0.99997 at
0.16, each matching its own exact value.

---

## 01_OneCylinder

One solid copper cylinder spanning a box, driven by 1 V across its end caps.
`a = 10 mm`, `L = 40 mm`, box 200 mm, a 36-gon cross-section.

The radius is the design decision. Skin depth in copper at 50 Hz is 9.346 mm,
so `a/delta = 1.07` and the current is **already** crowding toward the surface
at the lowest frequency of interest. The earlier 0.2 mm geometry had
`a/delta = 0.02` and could not show a skin effect at any frequency this solver
targets.

The mesh is graded from the conductor surface outward -- 2.0 mm at the surface,
16 mm in the far air -- and **capped at 2.0 mm inside the conductor**, which is a
separate field and the thing that was missing. See "The core sizing, and why it
had to be capped" below; briefly, the grading field is an unsigned distance from
the lateral surface, so without the cap it coarsens inward as well as outward and
leaves the axis at 4.4 mm elements.

| | |
|---|---|
| nodes | 4755 |
| unknowns | 55652 |
| nnz(L) after AMD | 38.43 M |
| factorization, internal solver | 90.0 s per frequency |
| factorization, `backend = mumps` | **2.14 s** -- 42x |
| one frequency end to end | 91.5 s internal, **3.1 s** MUMPS |

Element size inside the conductor, measured as the equivalent regular-tet edge:

| `r/a` | 0.0-0.1 | 0.1-0.2 | 0.3-0.5 | 0.5-0.7 | 0.7-0.9 | 0.9-1.0 |
|---|---|---|---|---|---|---|
| cells | 53 | 184 | 853 | 1350 | 2560 | 4492 |
| `<h>` mm | 2.60 | 2.68 | 2.70 | 2.63 | 2.30 | 1.65 |

**The direct solver, not gmsh, is what limits the mesh** -- or it was. Fill grows
roughly as `n^1.5`: an 8057-node version of this same geometry needed 4.3 GB and
12 min 19 s per frequency with the internal solver. MUMPS factorizes the present
matrix in 2.14 s, so that curve is no longer the binding constraint on
refinement; **the facet rule is.** `lc_core` cannot usefully go below about 2 mm
at `N = 36`, because the polygon facet is `2a sin(pi/N) = 1.743 mm` and a volume
size well under the facet width produces slivers. A finer interior needs a larger
`N` first.

### What is known about the answer

- `R = L/(sigma*A_poly) = 2.2114e-06 ohm`, against the **polygon** area
  `311.87 mm^2`, not `pi a^2 = 314.16`.
- `|E|` rises toward the surface. Measured **6.09 %** from axis to rim at 50 Hz
  against the exact **6.70 %** from `J0(kr)/J0(ka)`, so 91 % of the crowding is
  captured. It was 5.76 % (86 %) before the core was capped.
- `|B|` is zero on the axis, linear in `r` inside, `1/r` outside, peaking at the
  conductor surface.

### Viewing

```
output/
  potential.out  potential.vtk     Phi (about 67 MB per frequency in total)
  A_field.out    A_field.vtk       A   (gauge dependent -- see below)
  B_field.out    B_field.vtk       B,  plus the exact per-cell B_cell
  H_field.out    H_field.vtk       H
  E_field.out    E_field.vtk       E,  plus the exact per-cell E_cell
  J_field.out    J_field.vtk       J = sigma E, nodal and per-cell
```

One file per field. There is no combined file: to see two fields together,
open two of these in the same ParaView session -- they share a mesh, so the
views line up, and nothing has to be written twice to allow it.

Each `.out` is a plain table: index, position, that field's six real/imaginary
components, and the interface flag. Each `.vtk` carries the nodal field, its
magnitude, `material_interface`, and the per-cell array where one exists.

Colour `E_field.vtk` by `E_magnitude` or `B_field.vtk` by `B_magnitude`, and
**threshold `material_interface` to 0 first** -- the nodal average straddles
the conductor surface and is meaningless there. The per-cell arrays
(`E_cell_real`, `B_cell_real`, `J_real`) need no such care.

Use a **linear** colour scale for `|B|`: a log ramp compresses the `1/r` decay
and hides the ring entirely. See `docs/POSTPROCESSING_PLAN.md`.

`tools/pv_slice.py`, `tools/pv_bmag.py` and `tools/pv_skin.py` render the usual
cross-sections with pvbatch.

## Adding a case

Copy the shape of `01_OneCylinder`: a `.geo`, a `.aphi` naming it, an
`[output] directory = output`, and a header comment saying what the answer
should be and how you know. A case whose expected answer is not written down
is not a regression test.

Then add an `expected.txt` so `check.bat` picks the case up — it runs every
folder that has one and ignores every folder that does not. Copy the shape of
`02_Ansys_Cylinder_50Hz/expected.txt` for a resistive case or
`01_OneCylinder/expected.txt` for one with a real skin effect, and **say what
sets each tolerance**. A tolerance loose enough never to fail is not a test; the
way to confirm one is a test is to perturb the expected value, watch the check
fail, and put it back.

### Measured, on the mesh in this folder

`tools/pv_extract_rl.py output` and `tools/pv_ampere.py output`:

    conductor volume     1.250267e-05 m3   exactly the 36-gon
    terminal current     47199.07 - 136437.18 j A   |I| = 144370.55 A
    R                    2.264520e-06 ohm
    R (DC, exact)        2.206425e-06 ohm
    R_ac / R_dc          1.026330

`R` is **2.63 % above DC**, against the exact Kelvin-function result
`R_ac/R_dc = (u/2)[ber bei' - bei ber']/(ber'^2 + bei'^2)` at
`u = sqrt(2) a/delta = 1.513191`:

| | |
|---|---|
| `ber`, `bei` | 0.918265305, 0.567230807 |
| `ber'`, `bei'` | -0.215566682, 0.735963532 |
| Kelvin `R_ac/R_dc` | **1.0267245** |
| measured | **1.026330**, low by **0.038 %** |

Take those from the **defining series**, not from differencing `J0`: a central
difference put `bei'` out in the 6th digit, which moved this reference in the
5th, and it is the number the others are judged against.

That agreement is the point of the geometry: on the old 0.2 mm wire this rise was
about 1e-7 and indistinguishable from nothing. `check.bat` asserts it as
`R_over_Rdc`, and it is the **only** place in the suite where the skin effect
itself is checked against theory.

The `J` profile is checked too, and the reference needs care. `verify.py`
measures a mesh average over `r < 0.25a` divided by one over `r > 0.90a`. Where
the profile is flat that is the same as `J(0)/J(a)`; here it is not, because the
core cells sit well off the axis:

| | |
|---|---|
| Bessel at `r = 0` over `r = a` | 0.926497 &nbsp; ← **not** the right comparison |
| Bessel averaged over **the radii the mesh actually sampled** | **0.93993** |
| measured | 0.94229, high by **0.250 %** |

Comparing against the `r = 0` value would read as 1.7 % solver error when the
solver is within a quarter of a percent. **Because that reference is an average
over sampled radii it moves when the mesh moves** -- it was 0.94089 before the
core was capped -- so it has to be recomputed, not carried over, whenever
`cylinder.geo` changes.

### The core sizing, and why it had to be capped

`cylinder.geo` sized the whole mesh from one field: an **unsigned** `Distance`
from the wire's lateral surface, mapped through a `Threshold`. That distance
reads 10 mm on the axis -- the same as a point 10 mm out in the air -- so the
field coarsened **inward** as well as outward. An earlier version of the `.geo`
called that "correctly so, since nothing happens there", which is true at case
02's `a/delta = 0.16` and false here.

The fix is the same `MathEval` + `Restrict` + `Min` construction case 02 uses:
`Restrict` returns a huge size outside its volume, which is exactly what `Min`
wants, so `lc_core` binds only inside the wire and the air grading is untouched.
`Mesh.OptimizeNetgen` was switched on at the same time.

| | before | after | exact | |
|---|---|---|---|---|
| `R/R_dc` | 1.025955 | **1.026330** | 1.0267245 | 0.075 % -> **0.038 %** |
| `<J>core/<J>surf` vs own ref | +0.483 % | **+0.250 %** | -- | halved |
| axis-to-rim &#124;J&#124; | 5.76 % | **6.09 %** | 6.70 % | 86 % -> **91 %** captured |
| core cells, mid band | 28 / 2658 | **130 / 3686** | -- | none at all inside `r<0.1a` before |
| `<h>` in the core | 4.41 mm | **2.60 mm** | -- | |
| nodes | 4041 | 4755 | -- | +17.7 % |
| `Phi/V` vs `z/L` | 5.42e-01 | 4.88e-01 | -- | **barely moves -- see below** |

Every gauge-**independent** quantity roughly halved its error. `Phi` moved 10 %,
which is the point of the next section.

### Phi is gauge dominated here, and refining the mesh does not fix it

At `omega*L/R = 2.89` the reactance dominates. `E = -j*omega*A - grad Phi`, and
how the axial `E` splits between those two terms is set by that ratio: here `A`
carries most of it, `grad Phi` carries little, and **the continuum answer for
`Phi` is not `V z/L`.** So `Phi/V` deviating from `z/L` by 0.49 is not an error to
be refined away.

That is a claim worth testing rather than asserting, and the test is cheap: run
the **same mesh** at 1 Hz instead of 50 Hz, so the elements are identical and
only `omega` changes.

| same mesh | `omega*L/R` | `a/delta` | `Phi/V` vs `z/L` |
|---|---|---|---|
| 50 Hz | 2.891 | 1.070 | **4.88e-01** |
| 1 Hz | 0.059 | 0.151 | **1.93e-03** |

**253x**, on identical elements. Discretisation error cannot do that -- it would
be about the same at both frequencies. The deviation scales with `omega`, so it
is the gauge. The residual 1.9e-03 at 1 Hz *is* the mesh-limited part, and
capping the core improved that from 2.1e-03 -- about 10 %, which at 50 Hz is
invisible under 0.49.

This is why `01_OneCylinder/expected.txt` omits both `Phi` checks. `Phi` is gauge
dependent; `E`, `B`, `H` and `J` are not, and those carry the validation.

`cylinder_1hz_mumps.aphi` regenerates `output_1hz/` in about 3 s, so the
comparison can be redone whenever the mesh changes. **Both directories must come
from the same mesh for it to mean anything**, which is why that variant exists.

`|B|` against Ampere's law, per cell, ratio across the whole domain. Re-measured
on the present mesh -- the table that used to sit here had `r` running to 0.4,
which cannot be this geometry's millimetres when the conductor surface is at
`r = 10 mm`, so it was carried over from the old 0.2 mm wire:

     r (mm)        |B| meas   |B| Ampere    ratio   cells
      1.0 -  2.0     0.4518      0.4535     0.996      50
      3.0 -  4.0     1.0094      1.0173     0.992      96
      5.0 -  6.0     1.5670      1.5766     0.994     169
      8.0 -  9.0     2.5012      2.4892     1.005     351
      9.5 - 10.0     2.7727      2.7712     1.001     797   <- conductor surface
     10.0 - 11.0     2.7053      2.7511     0.983    1277
     12.0 - 14.0     2.2391      2.2570     0.992     530
     18.0 - 22.0     1.4481      1.4512     0.998     201
     35.0 - 45.0     0.7371      0.7316     1.008      99
     70.0 - 90.0     0.3681      0.3654     1.007     103

Linear in `r` inside, `1/r` outside, **within 0.4-1.7 % from `r = 1 mm` to
`r = 80 mm`**. Inside the conductor the comparison uses
`mu0 I r / (2 pi a^2)`, which treats the 36-gon as a circle, so the 1.7 % just
outside `r = a` is the polygon rather than solver error.

**Both scripts compare phasor magnitudes.** At 50 Hz here the current is
`47270 - 136555j`, so `Re(I)` is only 0.327 of `|I|`: comparing `|Re(B)|`
against `mu0 |I| / (2 pi r)` gives a ratio of 0.33 at every radius and looks
like the field is three times too small. It is not -- the two sides are
different quantities. This bit once and is worth not repeating.

## 03_Cylinder_1A_50Hz

The **dual** of case 02: the same conductor, the same mesh, the same frequency,
driven with **1 A** instead of 1 V. It is the first current-driven case this
project has solved -- the port machinery for it had existed since the DOF map
was written but nothing had ever exercised it.

It answers two things a voltage-driven case cannot.

**Port current self-consistency** (`SOLVER_PLAN.md` §8, item 7). Inject 1 A at
the top cap and exactly 1 A must be collected; charge does not accumulate in a
conductor. Measured **1.000000 A**.

**Does the dual agree?** `R` and `L` belong to the geometry and the material,
not to how it is driven, and the problem is linear:

| | case 02, 1 V | case 03, 1 A |
|---|---|---|
| `I` | 10180.56 A | **1.000000 A** |
| `V` | 1 V | **9.796995e-05 +7.098621e-06 j V** |
| `Z` | 9.796995e-05 +7.098621e-06 j | **9.796995e-05 +7.098621e-06 j** |
| `R` | 9.796995e-05 ohm | identical |
| `L` | 22.5956 nH | identical |

Identical to every digit printed.

**Two things worth knowing about it.**

`Phi` is *more* linear here than in case 02 -- worst deviation from `z/L` of
3.9e-06 against 8.3e-04, two hundred times better. A voltage port forces the
whole cap to one potential; a current port constrains only the total current
and lets the cap relax to whatever the field wants. The prescribed equipotential
is the thing bending `Phi` in case 02, not an error.

The **relative residual is 1.4e-09** where case 02 gives 6e-13, while the
backward error is 3.8e-21 against 3.7e-21 -- as good. That is exactly the case
`factorization.hpp` warns about: `||Ax-b||/||b||` is misleading when `||b||` is
small, and a 1 A drive makes a far smaller right-hand side than a 1 V one.
Judge this case by the backward error.

The mesh is **referenced, not copied**: `file = ../02_Ansys_Cylinder_50Hz/cylinder.msh`.
"Exactly the same mesh" is the point of the comparison, and a copy can drift.

---

## 04_Cylinder_SkinDepth

The `01_OneCylinder` rod -- copper, `a = 10 mm`, `L = 40 mm`, box 200 mm -- run
at **393 Hz** instead of 50, where `delta = 3.334 mm` and `a/delta = 3.000`.

![current density over the cross-section](04_Cylinder_SkinDepth/fig/j_cross_section.png)

That is the point of the case in one picture: the current is confined to a
surface layer about one skin depth thick, and the core carries almost none of
it. Quantitatively, `|J|` falls **3.7x** from surface to axis.

### Why a second high-`a/delta` case at all

| | 01 @ 50 Hz | **04 @ 393 Hz** |
|---|---|---|
| `a/delta` | 1.070 | **3.000** |
| `\|J(0)/J(a)\|` exact | 0.9265 | **0.2522** |
| `R_ac/R_dc` exact | 1.0267 | **1.7680** |
| the skin effect is | 2.6 % of `R` | **77 % of `R`** |
| measured `R_ac/R_dc` | 1.02633 | **1.77191** |
| error vs Kelvin | −0.038 % | **+0.220 %** |
| so the EFFECT is validated to | ~1.4 % | **~0.29 %** |

Matching Kelvin to 0.038 % sounds better than 0.220 % until you notice what
fraction of the answer the skin effect actually is. At `a/delta = 1.07` almost
all of `R` is just `R_dc`, which any solver that integrates `sigma` correctly
will get right; the thing under test is a 2.6 % correction. At `a/delta = 3` the
correction is most of the answer.

### The radial profile

![radial profile against Bessel](04_Cylinder_SkinDepth/fig/j_radial_profile.png)

38294 cells, one point each, against the exact `|J0(kr)/J0(ka)|` evaluated at the
radii actually sampled. **Worst band error 1.27 %**, and the error is monotonic
-- +1.27 % at the axis falling to 0.00 % at the surface. That shape is the
signature of an `h`-limited solution rather than a wrong one: the error is
largest exactly where `delta/h` is worst and vanishes where the mesh is finest.

![where each run sits on the Kelvin curve](04_Cylinder_SkinDepth/fig/kelvin_curve.png)

### How the mesh was sized, and the measurement it came from

The sizing was not guessed. `01_OneCylinder`'s mesh was run at three frequencies
first, and the error was read off against Bessel:

| f | `delta` | `a/delta` | `delta/h` | worst band err | `R` vs Kelvin |
|---|---|---|---|---|---|
| 50 Hz | 9.346 mm | 1.070 | 5.66 | 0.26 % | −0.038 % |
| 175 Hz | 4.996 mm | 2.002 | 3.03 | 2.10 % | −0.050 % |
| 300 Hz | 3.815 mm | 2.621 | 2.31 | 4.73 % | +0.386 % |

So the error is set by `delta/h`, and `delta/h >= 3` is what holds it near 1 %.

**Where the refinement had to go was the non-obvious part.** At 300 Hz the
per-band error reads 0.1 % in the annulus (`r > 0.8a`) and 4.7 % in the core
(`r < 0.3a`), which looks like "refine the core". That reading is wrong: `|J|` in
the core is *flat* there -- 0.3934, 0.3929, 0.3949 across `r/a` 0 to 0.3 -- so
the core is not failing to resolve local variation. It is inheriting accumulated
error from the region where the field actually decays, and `delta/h` was
**1.45 to 2.31 across the whole conductor**. The annulus's own relative error
only looks small because `|J|` is large there.

At `a/delta = 3` the decay region `r > a − 2*delta` is 89 % of the cross-section,
so there is nothing to gain by grading inside the conductor: `lc_core = lc_skin`
and the wire is meshed uniformly.

**`N` had to rise with it.** The facet is `2a sin(pi/N)`, and the volume element
size has to stay comparable to it or the mesh fills with slivers. At `N = 36` the
facet is 1.743 mm, so `h = 0.833 mm` would be less than half of it. `N = 76`
gives 0.827 mm, matched. This is the same rule that killed the `N = 96`
experiment in case 02 from the other direction -- there the facet was 0.098 mm
and the volume size 0.25 mm. The ratio is what matters, in either direction.

One thing that did **not** come out as designed: `lc_core = 0.833 mm` was meant
to give `delta/h = 4`, and the realised `delta/h` is **3.0–3.1**. gmsh's
characteristic length is not the equivalent-tet edge -- the realised `h` runs
about 1.3x the target, which case 01 also showed (target 2.0 mm, measured
2.60 mm). Reaching `delta/h = 4` would need `N = 96` and roughly twice the
elements. The 1.27 % that remains is mostly that.

### Cost, and why this case runs MUMPS by default

| | |
|---|---|
| nodes | 31922 |
| tets | 179477, of which 98204 in the wire |
| unknowns | **418318** |
| mesh | 24 s |
| factorization, MUMPS | **267 s** |
| one frequency end to end | 293 s |

**This is the first case that cannot practically run on the internal solver**, so
its `expected.txt` names `cylinder_393hz_mumps.aphi` rather than the plain input
-- the opposite of cases 01-03, which name the plain one on purpose so
`check.bat` exercises the default backend. `run_case.bat` reads `backend = mumps`
out of the input file and switches to the `build-mumps-omp` binary with the
oneAPI environment loaded; see `docs/MUMPS_SETUP.md`.

Note that 267 s is well above what linear fill extrapolation predicted (~60 s)
from the 55652- and 240990-unknown data points. Fill is not linear this far out.

### A caveat that refinement will not remove

**Kelvin's formula is for a circular conductor and this is a 76-gon.** At
`a/delta = 3` the current lives in a thin surface layer, so what matters is the
perimeter rather than the area: 62.8140 mm against the circle's 62.8319 mm, or
−0.028 %, which raises `R_ac` by about +0.028 %. That is a real part of the
+0.220 % and it will not go away with a finer mesh -- only with a larger `N`,
which also changes `R_dc`. A tolerance tight enough to exclude it would be
asserting something false.

### Regenerating

```
gmsh cylinder_skin.geo -3 -o cylinder_skin.msh
..\run_case.bat cylinder_393hz_mumps.aphi
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" make_plots.py output
```
