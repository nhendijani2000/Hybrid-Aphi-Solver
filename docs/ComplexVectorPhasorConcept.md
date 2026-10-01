# Complex vector phasors — magnitude, peak, and instantaneous value

A scalar phasor has one amplitude. **A vector phasor has two**, because its tip
traces an ellipse rather than oscillating along a line. Three different numbers
are all reasonably called "the magnitude of `E`", they are not equal, and a
post-processing UI that offers them has to say which is which.

This document derives all three from scratch, proves the relations between
them, measures how far apart they actually are in this project's own output,
and records which one every published figure and number in the suite currently
uses.

Companion to `FIELD_POSTPROCESSING.md`, which defines the field arrays this
operates on, and prerequisite for the planned `[postprocess]` input section.

---

## 1. The phasor and the time convention

A phasor vector is three complex numbers, which is the same thing as **two real
3-vectors**:

```
    Ê  =  P + jQ ,        P = Re Ê ,   Q = Im Ê ,    P, Q ∈ ℝ³
```

This solver uses `e^{+jωt}`. That is not a choice made in this document; it is
already in the assembly, where `src/element_matrix.cpp:22` builds

```
    α  =  jωσ − ω²ε           as  complex(−ω²ε, +ωσ)
```

with a **positive** imaginary part on the conduction term, and in
`src/problem.cpp:27`, where a port drive of amplitude `A` and `phase_deg = φ`
becomes `A(cos φ + j sin φ) = A e^{jφ}`. So a driven quantity is

```
    a(t)  =  Re{ A e^{jφ} e^{jωt} }  =  A cos(ωt + φ)
```

which is the `A·cos(ωt + φ)` convention the user-facing tools in this field
normally assume. Nothing below changes if the opposite convention `e^{−jωt}` is
ever adopted except the sign of every `Q`; the geometry is identical.

The physical field recovered from a stored phasor is therefore

```
    E(t)  =  Re{ Ê e^{jωt} }  =  P cos ωt  −  Q sin ωt
```

Write `u = ωt` throughout. The **phase** a user asks for is a value of `u`.

---

## 2. The trajectory is an ellipse

As `u` runs over one cycle, `E(u) = P cos u − Q sin u` is the image of the unit
circle `(cos u, sin u)` under the linear map `ℝ² → ℝ³` whose columns are `P` and
`−Q`. A linear image of a circle is an ellipse. It lives in the plane
`span{P, Q}`, and it degenerates:

| condition | trajectory |
|---|---|
| `P × Q ≠ 0` | a proper ellipse |
| `P ∥ Q` (includes `Q = 0`) | a **line segment** through the origin |
| `P = Q = 0` | a point |

This is called the **polarization ellipse**. Everything that follows is a
property of it.

A scalar phasor is the special case where the "vector" has one component: the
ellipse is always degenerate, always a segment, and the ambiguity below never
arises. That is why `Φ` needs none of this and `E`, `B`, `H`, `J`, `A` all do.

---

## 3. The magnitude over a cycle

```
    |E(u)|²  =  (P cos u − Q sin u) · (P cos u − Q sin u)

             =  |P|² cos²u  −  2(P·Q) cos u sin u  +  |Q|² sin²u
```

Apply `cos²u = (1 + cos 2u)/2`, `sin²u = (1 − cos 2u)/2`, `2 cos u sin u = sin 2u`:

```
    |E(u)|²  =  S  +  C cos 2u  −  D sin 2u
```

with the three scalars

```
    S  =  ( |P|² + |Q|² ) / 2            the cycle mean
    C  =  ( |P|² − |Q|² ) / 2
    D  =  P · Q
```

Collapse the oscillating part into a single cosine, `C cos 2u − D sin 2u =
R cos(2u + ψ)`, which requires `R cos ψ = C` and `R sin ψ = D`:

```
    R  =  √( C² + D² )           ψ  =  atan2( D , C )

    |E(u)|²  =  S  +  R cos(2u + ψ)
```

**`|E|²` oscillates at twice the field frequency** — which is why the extremes
come a *quarter* cycle apart, not a half.

---

## 4. The semi-axes

The extremes of `|E(u)|²` are `S ± R`:

```
    a  =  √( S + R )        semi-major,   attained at  u* = −ψ/2
    b  =  √( S − R )        semi-minor,   attained at  u* + 90°
```

`b` is real because

```
    S² − R²  =  ( (|P|²+|Q|²)/2 )²  −  ( (|P|²−|Q|²)/2 )²  −  (P·Q)²
             =  |P|²|Q|²  −  (P·Q)²
             =  |P × Q|²                                           ≥ 0
```

by Cauchy–Schwarz, with equality exactly when `P ∥ Q` — i.e. exactly the
degenerate, linearly polarized case of §2, where `b = 0`.

---

## 5. The two invariants

The step above gives `a²b² = S² − R² = |P × Q|²`, and `a² + b² = 2S`
immediately. So:

```
    ┌─────────────────────────────────────────────┐
    │   a²  +  b²   =   |P|²  +  |Q|²   =   N²    │
    │                                             │
    │   a   ·  b    =   |P × Q|                   │
    └─────────────────────────────────────────────┘
```

where `N` is defined in §6. These two symmetric functions determine `a` and `b`
completely: they are the roots of

```
    t²  −  N² t  +  |P × Q|²  =  0
```

The first invariant is the one that matters for the UI. It says

```
    N  =  √( a² + b² )
```

— the stored magnitude is the **hypotenuse of the two semi-axes**, never the
semi-major alone.

The second has a pleasant geometric reading: the ellipse's area is
`π a b = π |P × Q|`.

---

## 6. The three quantities, explicitly

Let `Ê = (Êx, Êy, Êz)` with `Êk = Pk + jQk`.

### Option 1 — magnitude at a phase  `|E(θ)|`

Take the instantaneous *real* vector first, then its ordinary Euclidean norm.

```
    Ek(θ)   =  Pk cos θ  −  Qk sin θ

    |E(θ)|  =  √( Σk ( Pk cos θ − Qk sin θ )² )
```

**Depends on θ.** This is what a field plot shows once a phase is set, and it is
the only one of the three that animates.

### Option 2 — complex magnitude  `N`

Take the modulus of each component *first*, then combine. The phase cancels.

```
    |Êk|  =  √( Pk² + Qk² )

    N  =  √( Σk |Êk|² )  =  √( Σk ( Pk² + Qk² ) )  =  √( |P|² + |Q|² )
```

**Independent of θ.** This is what `src/postprocess.cpp:766` writes as
`E_magnitude`, `B_magnitude`, … and what every figure in the regression suite
is coloured by.

### Option 3 — peak over the cycle  `a`

```
    a  =  max over θ of |E(θ)|  =  √( S + R )

        S = (|P|² + |Q|²)/2 ,    R = √( ((|P|²−|Q|²)/2)² + (P·Q)² )
```

**Independent of θ.** The largest value the field ever actually reaches.

Equivalently, from the invariants and the axial ratio `AR = b/a ∈ [0,1]`:

```
    a  =  N / √( 1 + AR² )
```

---

## 7. Ordering, bounds, and the energy reading

```
    0  ≤  b  ≤  |E(θ)|  ≤  a  ≤  N  ≤  √2 · a
```

The right-hand bound follows from `N² = a² + b² ≤ 2a²` since `b ≤ a`, with
equality iff `b = a` — circular polarization. So:

```
    1  ≤  N / a  =  √(1 + AR²)  ≤  √2
```

**The complex magnitude never understates the peak, and overstates it by at most
41.4 %.**

Averaging `|E(u)|²` over a cycle kills the `cos(2u + ψ)` term and leaves `S`:

```
    ⟨ |E|² ⟩  =  S  =  N² / 2        ⟹        N  =  √2 · RMS
```

**for any polarization whatsoever.** That is the cleanest statement of what `N`
*is*: it is an energy measure, not a geometric one — the amplitude a *linearly*
polarized field carrying the same mean-square would have. `a` is the geometric
peak. Both are legitimate; they answer different questions, which is precisely
why a UI must name which it is reporting.

---

## 8. Worked examples

| | `Ê` | `P` | `Q` | `N` | `a` | `b` | `b/a` | `N/a` |
|---|---|---|---|---|---|---|---|---|
| linear | `(3, 0, 0)` | `3x̂` | `0` | 3.0000 | 3 | 0 | 0.0000 | **1.0000** |
| elliptical | `(3, 2j, 0)` | `3x̂` | `2ŷ` | 3.6056 | 3 | 2 | 0.6667 | **1.2019** |
| circular | `3(x̂ + jŷ)` | `3x̂` | `3ŷ` | 4.2426 | 3 | 3 | 1.0000 | **1.4142** |

Check the middle row by hand:

```
    E(θ)  =  ( 3 cos θ , −2 sin θ , 0 )
    |E(θ)|²  =  9cos²θ + 4sin²θ  =  4 + 5cos²θ
```

| θ | 0° | 30° | 45° | 60° | 90° |
|---|---|---|---|---|---|
| `\|E(θ)\|` | 3.000 | 2.784 | 2.550 | 2.291 | 2.000 |

Option 1 sweeps 2 → 3. Option 3 is 3. Option 2 is 3.606, **above everything the
field ever reaches**. All three are correct answers to three different
questions. And `a² + b² = 9 + 4 = 13 = N²` ✓, `ab = 6 = |P × Q| = |(0,0,6)|` ✓.

**The circular row is the one to remember.** There `|E(t)| = 3` *constantly* —
the vector merely rotates — yet the complex magnitude reports 4.243. Nothing is
wrong: the field really does carry that mean-square energy. It is simply not a
peak.

### Why both of the first two exist

Take the linear row, `Ê = (E₀, 0, 0)`. Then

```
    |E(θ)|  =  E₀ |cos θ|
```

which is **zero at θ = 90°**. An animation of option 1 blinks to black twice per
cycle, everywhere at once. That is physically correct and visually useless.
Option 2 sits at `E₀` and never blinks — that is what it is for. It buys the
steadiness by being an energy measure, and the circular row is the price.

---

## 9. The ellipse as vectors, not just lengths

Lengths alone do not give the ellipse's orientation. To get the axes as
vectors, rotate the phasor by `e^{jδ}` until its real and imaginary parts become
orthogonal. With `Ê e^{jδ} = P' + jQ'`:

```
    P'  =  P cos δ  −  Q sin δ
    Q'  =  P sin δ  +  Q cos δ

    P'·Q'  =  C sin 2δ  +  D cos 2δ  =  0     ⟹     tan 2δ = −D/C
```

so `δ = −ψ/2`, the **same angle as the peak instant** `u*` of §4. Then `P'` and
`Q'` are the two semi-axis vectors, the longer being the semi-major:

```
    a = max(|P'|, |Q'|) ,     b = min(|P'|, |Q'|)
```

This also gives the handedness, via the sign of `P' × Q'` relative to the
propagation or symmetry direction, should it ever be wanted.

---

## 10. "Phase of a vector" does not exist

Each component has its own phase `atan2(Qk, Pk)`, and in general they differ —
that is exactly what makes the polarization elliptical. There is no single
number that is "the phase of `Ê`". A UI must therefore offer phase only for

- a **scalar** field (`Φ`), or
- a **named component** (`E_z`), or
- the **projection onto a user-given direction** `n̂`:
  `∠(n̂·Ê) = atan2(n̂·Q, n̂·P)`.

When the field *is* linearly polarized the component phases agree modulo 180°
and a single phase is meaningful — but the UI cannot assume that, and should
report the axial ratio alongside so the user can see whether the assumption
holds.

A related trap: `Φ` is **gauge-dependent** in this formulation. Case 01 measured
it changing by 254× between 1 Hz and 50 Hz on identical elements, purely from
the gauge. `E`, `B`, `H`, `J` are gauge-invariant; `Φ` and `A` are not. Any plot
of `Φ` should carry the gauge in its caption.

---

## 11. Numerical notes

**Computing `b` by the obvious formula loses half the digits.** `b = √(S − R)`
suffers catastrophic cancellation precisely when `S ≈ R`, which is the nearly
linear case — and that is most of the domain in every case in this suite. Use
the product invariant instead:

```
    a  =  √( S + R )                 well conditioned, no cancellation
    b  =  |P × Q| / a                stable
```

Measured on case 05, the two agree wherever `AR ≳ 1e-7`; below that the
cancelling form is reporting its own roundoff. Since `√ε ≈ 1.5e-8`, the
cancelling form can never resolve an axial ratio below about `1e-8`, so a field
that is *exactly* linear will show a spurious `AR` at that level.

**Clamp the discriminant.** `N⁴ − 4|P×Q|²` is non-negative in exact arithmetic
but can go slightly negative in floating point; take `max(·, 0)` before the
square root.

**Cost.** Both `a` and `b` are a handful of flops per cell on data already in
memory. There is no performance argument for not offering them.

---

## 12. What this project currently uses

**Option 2, everywhere, without exception.**

| | |
|---|---|
| `src/postprocess.cpp:766` | writes `√(Σ \|Êk\|²)` as `*_magnitude` |
| every regression-suite figure | ParaView Calculator `sqrt(X_real_X² + … + X_imag_Z²)` |
| every analysis script | `np.sqrt((re**2).sum(1) + (im**2).sum(1))` |

Consequences worth stating plainly:

**No figure in any case report is at a phase.** Not phase 0 — no phase at all.
They are phase-independent envelopes. A phase-0 plot would show `|P|`, which is
a different picture.

**The published analytic comparisons remain valid.** The analytic references are
themselves phasor magnitudes of linearly polarized fields — e.g. case 05's
`E(ρ) = I/(σρK)` descends from `Φ = Vθ/2π`, a single complex scalar times a real
unit vector, so `b = 0` and `a = N` identically. Comparing measured `N` against
`|analytic phasor|` was like for like. Where phase carries information it was
checked separately: case 04 asserts the core-to-surface phase lag at
**−134.24°** against Kelvin's **−133.85°**.

---

## 13. How much it matters here

Axial ratio over case 05's shipped mesh (90 mm domain), computed with the stable
form of §11:

| field | region | cells | median `b/a` | 90th | max `b/a` | max `N/a` |
|---|---|---|---|---|---|---|
| `E` | conductor | 23,409 | 9.40e-06 | 2.48e-05 | 1.23e-01 | 1.00758 |
| `E` | air | 194,310 | 4.20e-03 | 8.78e-03 | 7.84e-02 | 1.00307 |
| `B` | conductor | 23,409 | 4.74e-04 | 1.82e-03 | 2.44e-01 | 1.02929 |
| `B` | air | 194,310 | 6.57e-04 | 8.94e-03 | **9.21e-01** | **1.35953** |
| `J` | conductor | 23,409 | 9.40e-06 | 2.48e-05 | 1.23e-01 | 1.00758 |

`J` tracks `E` exactly, as it must where `J = σE`.

The 36 % outlier looks alarming until you ask where it sits:

| `\|B\|` band, fraction of max | cells | median `b/a` | max `N/a` |
|---|---|---|---|
| ≥ 10 % | 167,010 | 4.91e-04 | 1.0104 |
| 1 – 10 % | 30,764 | 1.85e-03 | 1.0475 |
| 0.1 – 1 % | 11,130 | 7.67e-03 | 1.2321 |
| 0.01 – 0.1 % | 7,468 | 2.21e-02 | **1.3595** |

The most nearly circular cell carries `|B| = 2.578e-07 T`, which is `8.1e-04` of
the peak — three decades below anything a plot renders. **Restricted to the
197,774 cells carrying ≥ 1 % of peak field, the worst overstatement anywhere in
the domain is 4.7 %, and the median is 1.00000.**

So for the cases in this suite the distinction is invisible wherever anyone
looks. It is not structurally zero, though, and it grows whenever a vector's
components acquire different phases: two current paths out of phase, conduction
mixing with displacement, strong skin effect near a corner.

### Phase-0 is not the same picture as the complex magnitude

| case | `\|E(0)\|/N` median | range across conductor cells |
|---|---|---|
| 05, 50 Hz | **1.0000** | 0.9922 – 1.0000 |
| 04, 393 Hz | **0.5938** | **0.0001** – 1.0000 |

Case 05's field is almost purely real (`ωL/R = 0.017`, so `|Q|` is 0.09 % of
`|P|`), so there a phase-0 plot and the complex magnitude coincide to within
0.8 %. That is luck.

Case 04 is the warning. With `a/δ = 3` the phase lags ~134° from surface to
core, so at `θ = 0` different radii are caught at different points in their
cycles and some cells read **one ten-thousandth** of their complex magnitude. A
phase-0 image of that case shows a dark ring where the field is crossing zero at
that instant. It is correct, and it looks nothing like the figure in the report.

---

## 14. Consequences for the `[postprocess]` design

1. **Default to option 2.** It reproduces every existing figure, it never
   blinks, and it is the only one of the three that needs no phase input.
2. **Offer all three under names that cannot be confused**, and grey out the
   phase box for options 2 and 3 — that teaches the distinction better than any
   tooltip.
3. **Expose `b/a` as a display option.** One plot then answers "is `magnitude`
   safe to read as a peak in my problem?" — cheaper than making users learn
   §§3–7.
4. **Do not silently redefine `magnitude` to mean the peak.** It would shift
   numbers in four committed case reports by amounts too small to notice and
   large enough to break their cross-checks.
5. **Document that setting a phase changes the picture**, with case 04 as the
   example, so it does not look like a regression to anyone who has read the
   reports.

| option | formula | phase box | reproduces existing figures |
|---|---|---|---|
| `complex_magnitude` | `√(\|P\|² + \|Q\|²)` | greyed | **yes** (default) |
| `magnitude_at_phase` | `\|P cos θ − Q sin θ\|` | active | no |
| `peak` | `√(S + R)` | greyed | within 0.001 % here |
| `axial_ratio` | `b/a`, via `b = \|P×Q\|/a` | greyed | diagnostic |
| `vector` | `P cos θ − Q sin θ` | active | — |

---

## 15. Verification

Every closed form above is checked against a brute-force sweep of 200,001
phases, on random `(P, Q)` pairs spanning six decades of scale, with a quarter
of the draws forced to exact linear polarization. The sweep uses no closed form,
so a disagreement would mean the derivation is wrong.

```bash
pvpython tools/verify_polarization_ellipse.py
```

| # | claim | worst error |
|---|---|---|
| 1 | `a` vs sweep maximum | 1.2e-10 |
| 2 | `b` vs sweep minimum | 1.6e-05 |
| 3 | `a² + b² = N²` | 6.6e-16 |
| 4 | `a·b = \|P × Q\|` | 1.4e-08 |
| 5 | the two formulas for `b`, against each other | 1.4e-08 |
| 6 | `⟨\|E\|²⟩ = N²/2` | 8.0e-16 |
| 7 | `δ = −ψ/2` orthogonalises, and `\|P'\|,\|Q'\|` are the semi-axes | 3.3e-16 |
| 8 | `a ≤ N ≤ √2 a` and `b ≤ \|E(θ)\| ≤ a` | 0 violations in 4,000 draws |

Three entries need their numbers explained rather than tolerated.

**Row 2 is the sweep's limitation, not the formula's.** `|E(u)|` is quadratic
near its minimum, so a grid of 200,001 points locates that minimum to about one
step — which is the 1.6e-05.

**Rows 4 and 5 are the `√ε ≈ 1.5e-8` cancellation of §11**, and finding it is
how that section came to exist. Row 4 deliberately uses the *independently*
derived `b = √(S − R)`; using the stable `b = |P×Q|/a` would make the identity
true by construction and the test would check nothing. Row 5 measures the two
formulas against each other directly and finds the same 1.4e-08.

**Row 3 is genuine despite using the stable `b`.** With `a = √(S+R)` and
`b = |P×Q|/a`, the claim `a² + b² = N²` is equivalent to
`a⁴ − N²a² + |P×Q|² = 0` — that is, it tests whether `a²` really is a root of
the quadratic of §5. It is, to machine precision.

An earlier version of this check reported "1.4e+08" for row 4 and three
violations in row 8. Both were faults in the test, not the mathematics: the
first divided by `|P × Q|`, which is *exactly zero* for the linear draws, and
the second compared the sweep minimum against a `b` that was pure roundoff at
the `1.5e-8` level. They are recorded here because the same two mistakes are
easy to repeat in an implementation.

---

## See also

- `docs/FIELD_POSTPROCESSING.md` — the field arrays, their definitions, and
  §6's note that `*_magnitude` is not the physical peak
- `docs/FORMULATION.md` §1 — the `e^{+jωt}` convention
- `docs/POSTPROCESSING_PLAN.md` — the staged plan this feeds
- `tools/verify_polarization_ellipse.py` — the checks of §15
- `src/postprocess.cpp:766` — `put_magnitude`, which writes option 2
- `regression_tests/04_Cylinder_SkinDepth/` — the case where phase actually
  matters, and the only one asserting a phase
