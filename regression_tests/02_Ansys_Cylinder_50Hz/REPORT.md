# 02_Ansys_Cylinder_50Hz — a slender copper cylinder at 50 Hz

A single copper rod driven by a 1 V terminal voltage, solved with the hybrid
A-Φ formulation and the tree-cotree gauge. Sized to match the regime of the
Ansys Maxwell A-Φ voltage example.

Regenerate every figure here with:

```bash
"C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" make_plots.py
```

---

## 1. Case

| | |
|---|---|
| conductor | copper, `sigma = 5.8e7 S/m`, `a = 1.5 mm`, `L = 40 mm`, a 40-gon |
| domain | 40 mm cube of air, `n x A = 0` on the outer boundary (flux tangential) |
| excitation | 1.0 V on the top face, 0.0 V on the bottom |
| frequency | 50 Hz |
| skin depth | `delta = 9.346 mm`, so `a/delta = 0.1605` — **no skin effect** |
| regime | `omega*L/R = 0.0725`, resistance dominated, so `Phi` is readable |
| formulation | first-order Whitney edge `A`, second-order P2 nodal `Phi` |
| mesh | 18994 nodes, 101873 tets, **Netgen-optimised**, 0.25 mm at the surface, 0.4 mm in the core |
| solve | 240990 unknowns, 503 s, 4.03 GB, backward error 4.8e-21 |

`a/delta` and `omega*L/R` are the same parameter, both scaling as `omega*a^2`.
A case cannot show a strong skin effect and a readable potential at once; this
one takes the resistive branch.

---

## 2. Mesh

| whole domain | the wire alone |
|---|---|
| ![](output/plots/01_mesh_domain.png) | ![](output/plots/02_mesh_wire.png) |

Graded from 0.25 mm at the conductor surface to 6 mm in the far air, then
optimised with Netgen. The element budget is heavily concentrated where the
field is: **6.7 %** inside `r = 1 mm`, **84.8 %** in the band `1.0–2.0 mm`, and
**3.8 %** beyond `r = 3 mm`.

`Mesh.OptimizeNetgen = 1` is the single most cost-effective setting in this
case. It improves element *shape* at fixed element size, and measured against
the same mesh without it, it reduced the tet count (its `CombineImprove` merges
elements), cut the peak-band azimuthal scatter by **24 %**, improved the error
in six of eight radial bands, and made the solve **26 % faster in 12 % less
memory** — 105715 → 101873 tets, 248294 → 240990 unknowns, 681 s / 4.56 GB →
503 s / 4.03 GB. Better element shape gives the sparse ordering less fill-in.

---

## 3. Potential

`Phi` over the lateral surface of the wire, seen side on — not a cross-section,
because the potential is driven along the axis.

| `Re(Phi)` | magnitude |
|---|---|
| ![](output/plots/03_phi_on_surface.png) | ![](output/plots/03b_phi_magnitude_surface.png) |

A clean `z/L` gradient from 0 to 1 V. Against the exact `z/L` over all 95149
wire nodes: worst deviation **8.26e-04**, mean **2.51e-04**.

`Phi` is gauge dependent — a different spanning tree shifts it by `-j*omega*psi`,
almost purely imaginary. `E`, `B`, `H` and `J` are not.

---

## 4. Magnetic flux density

Mid-length cross-section. Linear colour scale: a log ramp compresses the `1/r`
decay into the top few colours and hides the ring entirely.

| whole domain | zoomed to 3a |
|---|---|
| ![](output/plots/04_B_magnitude.png) | ![](output/plots/04b_B_magnitude_zoom.png) |

Zero on the axis, rising linearly inside the conductor, falling as `1/r`
outside — the peak sits at `r = a`.

### 4.1 Smoothed variant, for comparison only

Commercial tools smooth plotted fields by default, which is the leading
explanation for why their `|B|` cross-sections look cleaner. These apply one
smoothing pass so a like-for-like comparison can be made. **They are not the
solver's output** and each carries a banner saying so.

| honest default | smoothed, 1 pass |
|---|---|
| ![](output/plots/04b_B_magnitude_zoom.png) | ![](output/plots/08b_B_magnitude_zoom_SMOOTHED.png) |
| ![](output/plots/04_B_magnitude.png) | ![](output/plots/08_B_magnitude_SMOOTHED.png) |

Measured cost of that one pass:

| | unsmoothed | 1 pass | 2 passes | 4 passes |
|---|---|---|---|---|
| scatter, `1.20–1.49 mm` | 0.0296 | 0.0183 | 0.0179 | 0.0225 |
| error, `1.20–1.49 mm` | −0.39 % | −1.84 % | −2.43 % | −3.48 % |
| displayed peak (T) | 1.3385 | 1.2954 | 1.2716 | 1.2387 |
| peak change | — | −3.2 % | −5.0 % | −7.5 % |

38 % less scatter for 3.2 % of the peak — a better trade than on the pre-Netgen
mesh, where it cost 6.2 % for 14 %. The bias still grows while the scatter stops
improving past one pass. **Read numbers off the unsmoothed figures.** Smoothing
is not a route to a better answer — it is the control needed to compare like
with like if the reference picture is itself smoothed.

---

## 5. Electric field and current density

| magnitude, whole domain | magnitude inside the wire |
|---|---|
| ![](output/plots/05_E_magnitude.png) | ![](output/plots/06_E_in_wire.png) |

`E` is uniform inside the conductor at `V/L = 25 V/m` — there is no skin effect
at `a/delta = 0.16`. It is also one order smoother than `B`: `E = -j*omega*A -
grad(Phi)` is *linear* per tet because `grad(Phi)` comes from the P2 space,
while `B = curl A` is *constant* per tet. That is the whole reason `B` plots
rougher than `E`, in this code and in any other first-order edge-element code.

![](output/plots/07_J_vectors.png)

`J = sigma*E` inside the conductor, drawn on a longitudinal slice — on a
mid-length cut every arrow points at the viewer and nothing is visible.

---

## 6. Validation

Nothing below fits a free parameter. The current follows from `R_dc` under the
1 V drive: `I = 10207.3 A`, `B(a) = 1.3610 T`.

| quantity | result |
|---|---|
| `R` | 9.796995e-05 ohm against 9.796863e-05 DC exact, **1.3e-05** |
| `L` from `Im(Z)/omega` | 22.60 nH |
| `J(0)/J(a)` | 0.999968 against the Bessel 0.999959 |
| `Phi` vs exact `z/L` | worst 8.26e-04, mean 2.51e-04, over 95149 wire nodes |
| `Phi` at mid height | 84 nodes within 1 um, mean 0.499996, spread 6.02e-04 |
| `B` direction at the peak | azimuthal component 0.9999 of the total — **0.01 %** off |
| `B` vs exact, inside the conductor | every band within **0.5 %** |

`|B|` against the exact axisymmetric solution, by radius:

| band (mm) | mean, T | exact, T | error | azimuthal scatter |
|---|---|---|---|---|
| 0.30 – 0.60 | 0.4167 | 0.4163 | +0.24 % | 0.0477 |
| 0.60 – 0.90 | 0.6914 | 0.6884 | +0.40 % | 0.0309 |
| 0.90 – 1.20 | 0.9906 | 0.9929 | −0.22 % | 0.0210 |
| 1.20 – 1.49 | 1.2734 | 1.2722 | +0.10 % | 0.0138 |
| 1.49 – 2.00 | 1.2380 | 1.2633 | −1.99 % | 0.0307 |
| 2.00 – 3.00 | 0.8776 | 0.8869 | −1.15 % | 0.0559 |
| 3.00 – 5.00 | 0.5490 | 0.5601 | −2.01 % | 0.0768 |
| 5.00 – 8.00 | 0.3267 | 0.3349 | −2.37 % | 0.0868 |

Linear rise inside, `1/r` outside, correct absolute scale, over a 25x span in
radius. Inside the conductor every band is within 0.5 %.

---

## 7. Known limitations

**Azimuthal scatter in `B` of 1.4–3.1 % near the conductor.** The problem is
axisymmetric, so this is error. It is `O(h)` element noise from `B` being
constant per tetrahedron — the lowest-order quantity in the formulation.

**A coherent polygon harmonic is now visible, at 0.31 %.** Fourier analysis of
the azimuthal profile at `r = a` finds a peak at `m = 40`, the polygon order.
On every earlier mesh no such peak existed: it was buried under element-quality
noise. Netgen removed enough of that noise to expose it. At 0.31 % of the mean
against 1.4 % total scatter it is not what the eye sees, but it means that if
quality noise falls further, raising `N` finally becomes worthwhile — with
`lc_skin` lowered in step, so the facet stays comparable to the volume size.

**Far-field accuracy was traded for near-field resolution.** Beyond `r = 3 mm`
the error is ~2 % against ~1 % on a uniformly-graded mesh. That band holds
3.8 % of the elements and does not affect `R`, `L` or the conductor fields, but
it does degrade `tools/pv_ampere.py`.

**What was tried and did not work**, all measured and documented in
`docs/FIELD_POSTPROCESSING.md`: superconvergent patch recovery (made `B` more
than twice as rough), plot smoothing (costs peak amplitude for modest gain),
raising the polygon order alone (no effect on scatter, and slivers above
`N = 48`), and coarsening the far field further (at most a few per cent of the
mesh to reclaim).

**Uniform refinement is reaching its limit.** `lc_skin` 0.7 → 0.35 → 0.25 gave
peak-band scatter 0.029 → 0.0215 → 0.0181, ratios x0.74 and x0.84 where `O(h)`
allows x0.50 and x0.71. Netgen then took it to 0.0138 for free, which is a
larger gain than the last refinement step bought — element *shape* was the
cheaper lever all along. Going further by refinement alone needs a solver
holding more than the ~250k unknowns this direct factorisation manages.
