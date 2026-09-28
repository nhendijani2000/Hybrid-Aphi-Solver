# Regression tests

Each folder is one self-contained case: its geometry, its input files, and its
own `output/`. Nothing here writes outside its own folder.

```
regression_tests/
  run_case.bat            the runner
  01_OneCylinder/
    cylinder.geo          geometry
    cylinder.msh          generated from the .geo
    cylinder_50hz.aphi    one frequency
    cylinder_sweep.aphi   50 Hz to 500 Hz
    output/               everything the solver writes
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

## 01_OneCylinder

One solid copper cylinder spanning a box, driven by 1 V across its end caps.
`a = 10 mm`, `L = 40 mm`, box 200 mm, a 36-gon cross-section.

The radius is the design decision. Skin depth in copper at 50 Hz is 9.346 mm,
so `a/delta = 1.07` and the current is **already** crowding toward the surface
at the lowest frequency of interest. The earlier 0.2 mm geometry had
`a/delta = 0.02` and could not show a skin effect at any frequency this solver
targets.

The mesh is graded from the conductor surface outward: 2.0 mm at the surface,
16 mm in the far air, with 60 % of the nodes in the two radial bands straddling
`r = a`. Node density at the surface is 600x that at the box wall.

| | |
|---|---|
| nodes | 4041 |
| unknowns | 49276 |
| nnz(L) after AMD | 39.4 M, 751 MB |
| factorization | ~93 s per frequency |
| one frequency, end to end | ~3 min |

**The direct solver, not gmsh, is what limits the mesh.** Fill grows roughly as
`n^1.5`: an 8057-node version of this same geometry needed 4.3 GB and took
12 min 19 s for a single frequency. If you refine `lc_skin`, expect that curve.

### What is known about the answer

- `R = L/(sigma*A_poly) = 2.2114e-06 ohm`, against the **polygon** area
  `311.87 mm^2`, not `pi a^2 = 314.16`.
- `|E|` rises toward the surface. Measured 5.5 % from axis to rim at 50 Hz,
  against 0.9414 from the exact `J0(kr)/J0(ka)` for the same radial bin
  centres -- 0.7 %.
- `|B|` is zero on the axis, linear in `r` inside, `1/r` outside, peaking at the
  conductor surface.

### Viewing

```
output/
  potential.out  potential.vtk     Phi (about 67 MB per frequency in total)
  A_field.out    A_field.vtk       A   (gauge dependent -- see below)
  B_field.out    B_field.vtk       B,  plus the exact per-cell B_cell
  H_field.out    H_field.vtk       H
  E_field.out    E_field.vtk       E,  plus E_cell and J = sigma E
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

### Measured, on the mesh in this folder

`tools/pv_extract_rl.py output` and `tools/pv_ampere.py output`:

    conductor volume     1.250267e-05 m3   exactly the 36-gon
    terminal current     47269.85 - 136555.16 j A
    R                    2.263692e-06 ohm
    R (DC, exact)        2.206425e-06 ohm
    R_ac / R_dc          1.0260

`R` is **2.60 % above DC**, against **2.67 %** from the exact Kelvin-function
result `R_ac/R_dc = (u/2)[ber bei' - bei ber']/(ber'^2 + bei'^2)` at
`u = sqrt(2) a/delta`. That agreement is the point of the geometry: on the old
0.2 mm wire this rise was about 1e-7 and indistinguishable from nothing.

`|B|` against Ampere's law, per cell, ratio across the whole domain:

     r (mm)        |B| meas   |B| Ampere    ratio
     0.05 - 0.075    0.9298      0.9319     0.998
     0.175 - 0.200   2.6619      2.7001     0.986     <- conductor surface
     0.250 - 0.275   2.1875      2.2022     0.993
     0.375 - 0.400   1.4919      1.4926     1.000

Linear in `r` inside, `1/r` outside, within 1-3 % throughout.

**Both scripts compare phasor magnitudes.** At 50 Hz here the current is
`47270 - 136555j`, so `Re(I)` is only 0.327 of `|I|`: comparing `|Re(B)|`
against `mu0 |I| / (2 pi r)` gives a ratio of 0.33 at every radius and looks
like the field is three times too small. It is not -- the two sides are
different quantities. This bit once and is worth not repeating.
