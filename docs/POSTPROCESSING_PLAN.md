# Post-processing, probes and validation — plan

Phase 05b. The solve is verified (`docs/SOLVER_PLAN.md` §16–§19: reconstruction
to 1.4e-16, ordering-independent, agreeing with a dense LU, backward errors at
1e-25). **The physics is not verified at all.** Nothing has been extracted from a
solution, so every claim about ports, the gauge, the weak form and the sign
conventions rests on structural and algebraic checks alone.

This is the phase that can find out. It comes before hardening the solver
(`SOLVER_PLAN` §20) for a plain reason: a solver that is numerically exact on a
problem whose answer is wrong is not worth hardening.

---

## 1. What a solution vector contains

`x` is indexed by `DofMap`'s global numbering, `[a | phi | V]`:

| block | indices | meaning |
|---|---|---|
| `a` | `0 … num_a-1` | one coefficient per **free** edge, for the **global** edge direction |
| `phi` | `num_a … num_a+num_phi-1` | one per **free** P2 node — vertices *and* edge midpoints |
| `V` | the last `num_ports` | one per port |

Three things a reader of `x` must not forget, each of which is a way to get a
plausible wrong field:

- **Eliminated DOFs are not in `x`.** A Dirichlet or tree edge has `a = 0`; an
  absent Φ node has no value; a prescribed Φ node carries `phi_fixed`, not an
  entry. `DofMap::local_dofs` already returns all of this per tet, and
  reconstruction must go through it rather than indexing `x` directly.
- **The edge sign is not optional.** `local_dofs` returns `coeff = ±1` from
  `tet_edge_signs`, and roughly 58 % of (tet, local edge) pairs are reversed.
  Dropping it flips more than half the field.
- **Under `ScaledPhi` the unknown is `Φ' = Φ/(jω)`, not `Φ`.** Extraction must
  multiply by `jω`. This is listed as a negative control in `SOLVER_PLAN` §8 and
  is the single likeliest source of a wrong answer that still looks smooth.

---

## 2. Fields, and where each one naturally lives

From `docs/FORMULATION.md` §1, under `e^{+jωt}`:

    A    = Σ_e a_e W_e(x)                    Whitney edge functions
    B    = curl A = Σ_e a_e curl W_e         **constant per tet**
    H    = B / μ
    E    = -jω A - grad(Phi)                 FORMULATION.md line 50
    J    = sigma E                           in conductors only

The natural home of each quantity is not the same, and pretending otherwise is
how post-processing loses accuracy silently:

| quantity | natural representation | at a vertex |
|---|---|---|
| `Φ` | **nodal** — a P2 DOF at every vertex and every edge midpoint | **exact, no averaging** |
| `A` | linear within a tet, tangentially continuous across faces | evaluate per tet, then average |
| `B`, `H` | **constant per tet** | inherently an average |
| `E` | linear within a tet (`grad Φ` of a quadratic) | evaluate per tet, then average |

### Φ needs no averaging, and averaging it would be wrong

This is a deliberate difference from the common practice of computing a field at
the vertices and taking the mean of the two endpoint values for a mid-edge node.
That is the right thing for `B`, `H` and `E`, which have no mid-edge degree of
freedom. **It is the wrong thing for Φ here**, because Φ lives in a *P2* space:
the mid-edge value is a genuine unknown that the solve determined. Averaging the
two endpoints would replace the quadratic with its linear interpolant and discard
the very term that makes the space second order.

Solvers that average Φ at mid-edge nodes do so because their Φ is P1 and the
mid-edge node is only there to draw a curved element. Ours is not, and the plan
is explicit about it so nobody "fixes" it later.

### Averaging for `B`, `H`, `E` is smoothing, not accuracy

`curl W_e` is constant per tet, so `B` is genuinely piecewise constant. A
per-vertex `B` is a **post-processing choice**, volume-weighted over the incident
tets:

    B_vertex = ( Σ_t vol_t * B_t ) / ( Σ_t vol_t )

It is smoother to look at and no more accurate than the per-tet values, and the
output will carry both so nobody mistakes the smoothing for convergence. For a
mid-edge node, the mean of its two endpoint values — which is what the
conventional recipe does, and is right for these quantities.

---

## 3. The API

```cpp
/// Everything reconstructed at one point, in one tet.
struct FieldSample {
    Vec3 position;
    std::complex<double> phi;
    std::array<std::complex<double>, 3> a, b, h, e;   // A, B, H, E
    int tet = -1;                                     // -1 if the point is outside
};

/// Reconstruct at a point. Locates the containing tet, then evaluates.
FieldSample sample_fields(const Mesh&, const BoundProblem&, const DofMap&,
                          const Solution&, const Vec3& point);

/// A probe along a straight line: `count` samples from `from` to `to`.
/// The workhorse for validation -- Phi down the cylinder's axis must be linear,
/// and B around the loop must fall as 1/r outside it.
std::vector<FieldSample> probe_line(..., const Vec3& from, const Vec3& to, int count);

/// Per-vertex and per-tet fields for the whole mesh.
struct FieldOutput {
    std::vector<std::complex<double>> phi_vertex;   // exact P2 values
    std::vector<std::complex<double>> phi_edge;     // exact P2 values, mid-edge
    std::vector<std::array<std::complex<double>, 3>> b_tet, e_tet;   // as computed
    std::vector<std::array<std::complex<double>, 3>> b_vertex, e_vertex, h_vertex, a_vertex;
};

FieldOutput compute_fields(const Mesh&, const BoundProblem&, const DofMap&, const Solution&);

/// Writes a Gmsh post-processing file the mesh can be viewed with directly --
/// `$NodeData` for the per-vertex quantities and `$ElementData` for the per-tet
/// ones, so both are visible and can be compared.
void WriteSolution(const std::string& path, const Mesh&, const FieldOutput&);
```

`Solution` bundles what a solve produced — the vector, the conditioning it was
obtained under (so `ScaledPhi`'s `jω` is applied once, in one place, rather than
at every call site), and the frequency.

Gmsh's own format is chosen over VTK because the mesh came from Gmsh, which is
already installed and can open the result without a second tool.

---

## 4. Extraction: R, L and C

Two independent routes wherever one exists, because agreement between them is
worth more than either number alone.

### The cylinder at DC — `examples/cylinder_box.aphi`

| quantity | route A | route B |
|---|---|---|
| `R` | `V/I` from the port unknown | `P/I²` with `P = ∫ σ\|E\|² dV` |
| `L` | `2·W_m/I²` with `W_m = ∫ \|B\|²/(2μ) dV` | `(1/I²) ∫ A·J dV` |

Targets, from `examples/cylinder_box.aphi` and its derivation:

    R = 0.1388 mOhm      exact to round-off against the MESHED cross-section
    L = 0.3870 nH        within ~0.5 %

`R` being exact is the sharpest check in the project: the DC solution `Φ = V·z/l`
is linear, hence lies in the P2 space, and the 24-gon's planar vertical lateral
faces make `dΦ/dn = 0` hold there exactly. If `R` is not right to several digits,
something upstream is wrong — and the polygon area, not `πa²`, is the denominator.

Also free, and worth asserting: **the current at P2 must be −1 A.** A voltage
port's current is no longer a row residual (`SOLVER_PLAN` §10), so this requires
re-forming the terminal's current balance — which is work this phase owes anyway.

### C is *not* extractable from this configuration, and that needs saying

The cylinder is driven as a **series** element: 1 A in at the bottom cap, 0 V at
the top. That measures a series impedance `R + jωL`. The capacitance between wire
and box wall is a **shunt**, and with `n×A = 0` on the walls and no port on them
there is no shunt path in the model at all. No series measurement can produce it.

Three ways to get `C`, in increasing honesty of effort:

1. **From electric energy on an AC run**: `C = 2·W_e/ΔV²` with
   `W_e = ∫ (ε/2)\|E\|² dV` and `ΔV` the wire-to-wall potential difference taken
   from the solution. Workable, but `ΔV` is a derived quantity here because Φ is
   free on the box.
2. **A dedicated electrostatic problem**: σ = 0 everywhere, wire at 1 V, walls at
   0 V. Clean and unambiguous, and needs a `boundary_voltage` port on the box
   wall, which the input format already supports.
3. **A genuine two-port**, which is what a real extraction flow does.

The target, from `L·C = μ₀ε₀l²` for a coaxial geometry using the **external**
part of `L` only (0.337 nH of the 0.387 nH; the remaining 0.050 nH is internal to
the conductor and pairs with nothing):

    C = 2π ε₀ l / ln(1.0787 W / 2a) = 33.0 fF

Checked: `0.337 nH × 33.0 fF = 1.112e-23` against `μ₀ε₀l² = 1.113e-23`.

**Recommendation: route 2.** It is the only one that measures `C` rather than
inferring it, and it costs one small input file.

### The loop — `examples/loop_internal_port.aphi`

No closed form as tidy as the cylinder's, so the checks are relational rather
than absolute:

- `L = 2·W_m/I²` compared against the standard thin-ring formula
  `L ≈ μ₀R[ln(8R/a) − 2]`, which for a square cross-section is only indicative —
  expect tens of percent, and say so rather than pretending otherwise.
- **`Φ` is single-valued and jumps only across the cut.** Walk a circle inside
  the ring: Φ must be continuous everywhere except one step of exactly `V`.
- **`B` circulates.** On a loop threading the ring's hole, `∮H·dl = I` to
  discretisation error. This is Ampère's law as a direct check on the gauge and
  the sign conventions, and it is the loop fixture's real purpose.
- Current continuity: `∫ J·dS` over any cross-section of the ring equals `I`.

---

## 5. What can go wrong, and the controls for it

Every one of these produces a smooth, plausible field:

| mistake | control |
|---|---|
| drop the edge sign from `local_dofs` | `∮H·dl` collapses; `B` reverses in ~58 % of tets |
| forget `Φ = jω·Φ'` under `ScaledPhi` | the three conditionings disagree on `Φ` by a factor of `ω` |
| average Φ at mid-edge nodes | `R` loses accuracy; the linear-Φ exactness argument breaks |
| use `E = +jωA − ∇Φ` | `R` from `P/I²` disagrees with `V/I` |
| treat `B` as per-vertex without volume weighting | `L` from energy drifts with mesh grading |
| index `x` directly instead of via `local_dofs` | eliminated DOFs read as garbage |

**The strongest single check is the three conditionings.** `Natural`,
`RowScaled` and `ScaledPhi` describe one problem and `test_assembly` already
proves their matrices are related exactly. So the extracted `R`, `L` and the
fields must agree across all three — which exercises assembly, solve *and*
extraction end to end, and no single-formulation test can.

---

## 6. Order

1. `Solution`, `sample_fields`, and the tet-location it needs. Tested on a hand-
   built mesh where the fields are known.
2. `compute_fields`, with the per-vertex averaging and the exact Φ handling.
3. `WriteSolution` to Gmsh format; open the cylinder's Φ and check it by eye.
4. `probe_line`; assert Φ is linear down the cylinder's axis.
5. Extraction of `R` and `L`, both routes, and the port current.
6. **The cylinder DC milestone**: `R = 0.1388 mΩ`, `L = 0.3870 nH`, `I(P2) = −1 A`.
7. The electrostatic `C` run, against 33.0 fF.
8. The loop: Ampère's law, the Φ jump across the cut, current continuity.
9. The three-conditioning agreement, on both fixtures.

Step 6 is the milestone. Everything before it is machinery; everything after it
extends the validation to a second geometry and a second quantity.

---

## 7. An external reference: Ansys Maxwell, 1 A loop at 50 Hz

Supplied Sept 2026. A ring driven by a 1 A current source at 50 Hz; the reported
terminal voltage, with `I = 1 A` so `Z = V`:

    V = 0.00034575462656563323623 + j 1.2395463798949440417e-5  V

| | |
|---|---|
| **R** = Re(Z) | **345.75 µΩ** |
| X = Im(Z) | 1.2395e-5 Ω |
| **L** = X/ω | **39.456 nH** |
| Q = X/R | **0.0359** |

This is worth more than the analytic loop formula in §4, because it is an
independent solver's answer to the same physics rather than a closed form with
its own approximations.

### The dimensions were lost, and are recoverable only up to σ

Inverting `R = 2r/(σa²)` and `L = μ₀r[ln(8r/a) − 7/4]` — the uniform-current form,
justified below — gives one equation short of a unique answer:

| σ | major radius `r` | wire radius `a` | outer diameter |
|---|---|---|---|
| **5.8e7 (copper)** | **11.60 mm** | **1.075 mm** | **25.3 mm** |
| 3.5e7 (aluminium) | 12.58 mm | 1.442 mm | 28.0 mm |
| 1.0e7 | 15.83 mm | 3.03 mm | 37.7 mm |

All three reproduce `R` and `L` to the digits given. Copper is the likeliest —
its 25.3 mm outer diameter matches the scale bar in the supplied figure, where
the ring sits just inside the 3 cm mark. **Confirming the material pins the
geometry completely**, and is the one piece of information worth recovering.

### Why this is a good validation case

**It is effectively a DC problem.** Skin depth in copper at 50 Hz is 9.35 mm,
**8.7× the wire radius**, so the current is essentially uniform. Nothing about
the answer depends on resolving a skin layer — and it probes the ω → 0 limit,
which is the regime `docs/CONDITIONING.md` exists for and the one where an A-Φ
formulation is hardest.

**The imaginary part is the demanding half.** `Q = 0.036` means the inductive
term is 28× smaller than the resistive one, so matching `L` to a few percent
requires `Im(V)` to be right when it is 3.6 % of `Re(V)`. A formulation can get a
dominant quantity right by accident; it cannot do that to a term this small.

### What it would take to use it

- **A new fixture.** The inversion assumes a **circular** cross-section;
  `loop_cut.msh` is an extruded annulus with a **rectangular** one, so it cannot
  reproduce this. A revolved polygonal cross-section is needed — the same
  two-half-annuli construction as `tools/loop_cut.geo`, with the cross-section
  revolved rather than extruded.
- **The material confirmed**, or the case treated as three candidate geometries.

### The caveat on the recovered numbers

`L = μ₀r[ln(8r/a) − 7/4]` is the thin-ring approximation and here `a/r = 0.093`,
so it is good to roughly 1 %. That uncertainty attaches to the **dimensions**,
not to `R` and `L`, which come straight from the reported voltage and are exact
to the digits supplied.

### Copper confirmed (Sept 2026) — and what is still uncertain

The material is **copper**, so:

| | |
|---|---|
| major radius `r` | **11.60 mm** |
| wire radius `a` | **1.075 mm** |
| outer diameter | **25.3 mm** |

The two equations are **not equally trustworthy**, and it matters which is
leaned on:

**`R` is domain-independent and exact.** `R = 2r/(σa²)` with σ = 5.8e7 and
R = 345.75 µΩ gives

    r / a² = 10027 per metre

whatever the outer boundary does. Nothing about the air region, the truncation
or the boundary condition can change a DC resistance. This is the equation to
trust and the one a fixture must reproduce first.

**`L` carries a domain caveat.** The inversion assumed the reported 39.456 nH is
the **free-space** thin-ring value. The supplied figure shows the air region is a
**cylinder roughly 3 cm across against the ring's 2.5 cm** — only about 1.2× the
ring diameter. A boundary that close generally suppresses `L` relative to free
space, so the true ring may be slightly larger than the numbers above. If the
free-space `L` were 5 % higher, the pair moves along `r = 10027 a²` to
r = 12.08 mm, a = 1.098 mm — 4 % and 2 %.

So the honest statement of the fixture's geometry is: **`r/a² = 10027 m⁻¹`
exactly, with `(r, a) = (11.60, 1.075) mm` as the best point estimate**, carrying
a few percent of uncertainty that only Maxwell's domain size and boundary
condition would remove.

**The air region is a cylinder, not a box.** A fixture built to match this should
use a cylindrical air domain rather than the square box of `loop_cut.geo`, since
with the boundary this close its shape is not a detail.

### What to compare, in order of confidence

1. **`R` = 345.75 µΩ.** Domain-independent, skin-effect-free at 50 Hz
   (δ = 9.35 mm = 8.7 a), and the DC solution is the one our formulation should
   reproduce most accurately. **Expect agreement to a fraction of a percent**, and
   treat a larger discrepancy as a defect rather than a modelling difference.
2. **`L` = 39.456 nH.** Expect a few percent, limited by the domain question above
   rather than by either solver. Worth reporting the discrepancy rather than
   tuning the geometry until it vanishes.
3. **`Q` = 0.0359.** The ratio is a cleaner target than `L` alone, since both
   parts come from the same solve and the same geometry.
