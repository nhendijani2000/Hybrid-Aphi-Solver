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
