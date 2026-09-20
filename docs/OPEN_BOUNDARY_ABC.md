# Phase 07 open boundary: a first-order ABC for A-Phi

## Decision

Phase 07's interim, approximate open-boundary treatment is a **first-order
absorbing boundary condition (ABC)**, not a PML. This document records the
literature check behind that choice and the from-scratch derivation of the
ABC written directly in terms of the potentials **A** and Phi, since no such
formula was found already published.

## Why ABC over PML (recap)

PML is the more common default in general full-wave E/H-field codes, but two
things push against it here specifically:

1. **PML has known low-frequency degeneracy** (the complex coordinate-stretching
   parameters become ill-conditioned as omega -> 0). That would put a second
   low-frequency failure mode right next to a formulation whose entire premise
   is all-frequency stability from DC through THz -- a direct conflict with
   the reasoning behind every other choice in this project (keeping the
   epsilon terms, preferring ACA over MLFMA, etc.).
2. **PML is the larger implementation lift for this specific solver.** It
   needs an added volumetric absorbing layer (a mesh change), anisotropic
   complex-tensor assembly inside that layer, and empirical
   thickness/grading-profile tuning -- on top of translating a technique
   whose published form is for E/H fields onto the coupled A-Phi system while
   keeping the chosen gauge (Sec. below) valid inside the stretched region.
   An ABC costs one short field-to-potential substitution and drops in as a
   single boundary term added to the existing weak form: no new DOFs, no mesh
   change.

An ABC's trade-off -- accuracy that degrades away from near-normal incidence,
and a need for reasonable standoff distance from the scatterer -- is a
bounded, documentable limitation, not a new failure regime. Phase 12 replaces
this whole boundary treatment with an exact FEM-BI coupling once the core
solver is proven, so the interim treatment's job is just to be honest about
its error budget, per the roadmap's own requirement.

## Literature check: does a published A-Phi-native ABC exist?

Searched specifically for an absorbing or radiation boundary condition
written directly in terms of **A** and Phi, rather than applied to derived E/H
fields. Result: **no such formula was found.** What exists instead:

- **Zhang, Na, Jiao & Chew, "An A-Phi Formulation Solver in Electromagnetics,"
  arXiv:2207.02260 (2022)** -- the closest existing broadband A-Phi solver
  (already cited in `docs/REFERENCES.md` for its DEC discretization, with no
  boundary-integral coupling). Its own numerical validation does not use an
  A-Phi-native ABC: a rod-antenna example uses a simple impedance boundary
  condition (IBC) described as "a simple absorbing boundary condition," and a
  nano-scale rod-antenna example uses plain PEC truncation "for simplicity."
  PML is discussed narratively but its implementation is deferred to a
  separate discrete-exterior-calculus reference, with no A-Phi equations
  given. This is consistent with the gap already noted in
  `docs/FULL_WAVE_SCOPE.md`'s novelty note for the BEM-coupling question --
  the open-boundary treatment for A-Phi solvers generally falls back to
  simplified field-based conditions or is left unspecified.
- **Ciuprina & Sabriego, "Electric circuit element boundary conditions for
  electromagneto-quasistatic and full wave models in A, phi potentials,"
  J. Math. Industry, vol. 14, art. 27 (2024)** -- covers port/circuit boundary
  conditions (equipotential terminals, voltage/current excitation) for A-Phi
  FEM, not open/radiating boundaries. Not applicable here, but worth knowing
  about for Phase 06's port abstraction.
- **Chervyakov, "On the use of mixed potential formulation for finite-element
  analysis of large-scale magnetization problems," arXiv:2307.12308 (2023)**
  -- magnetostatic only; truncates with PMC (n x H = 0) and magnetic
  insulation (n . B = 0) on a large sphere, not a radiation condition. Not
  applicable to a full-wave radiating problem.

**Conclusion:** the ABC below is a from-scratch derivation applying a
standard, widely-published, peer-reviewed field-based result (the
Silver-Muller/Sommerfeld first-order ABC, verified against Chen & Chew 2017,
Eq. 83, above) to this project's own A-Phi system, via the standard potential
definitions already in use elsewhere in this roadmap. No employer-authored or
otherwise non-public material entered this derivation.

## Literature check: PML for A-Phi specifically (Sept 2026)

A separate, dedicated search for a PML formulated on A-Phi potentials (as
opposed to E/H fields) turned up nothing published either -- if anything, a
narrower gap than the ABC case:

- **Chen & Chew, J. Comput. Phys. 350 (2017)** (above) explicitly implements
  PEC, PMC, Dirichlet, periodic, and ABC boundary treatments in a DEC/E-H
  setting, but **does not implement PML** -- confirmed directly from the
  paper's boundary-condition section.
- **Zhang, Na, Jiao & Chew, arXiv:2207.02260 (2022)** discusses PML narratively
  but defers its actual implementation to reference [9] in that paper's own
  bibliography -- J. Rabina, "On a numerical solution of the Maxwell equations
  by discrete exterior calculus," Ph.D. dissertation, University of
  Jyvaskyla, 2014. That dissertation was not directly accessible to verify,
  but everything else found about this DEC lineage (Chen & Chew 2017, and
  Zhang/Na/Jiao/Chew's own paper) formulates its Maxwell solver in E/H
  differential forms, not A-Phi -- there is no indication anywhere in this
  citation chain that a PML was ever written directly on A-Phi potentials,
  only that a DEC-framework PML exists for the field-based version of the
  same broader numerical lineage.
- General searches for "vector potential" + PML in the wider literature
  return eddy-current/magnetostatic vector-potential-only formulations (no
  coupled scalar potential, no radiation, e.g. A-V or T-Omega eddy-current
  methods) -- a different, non-radiating problem class, not evidence of an
  A-Phi-plus-PML radiation treatment.

**Conclusion:** no A-Phi-native PML exists in the literature either. This
doesn't change the Phase 07 decision (ABC over PML) -- if anything it
reinforces it: adopting PML here would mean adapting an *unpublished-for-A-Phi*
technique with a harder physical translation (anisotropic complex-stretched
material tensors threaded through the coupled curl-curl/gradient system while
preserving the gauge) and a known low-frequency failure mode, versus the ABC's
short, verified field-to-potential substitution.

## Why this gap exists (Sept 2026)

Worth recording deliberately, since it bears directly on the Phase 12
"candidate paper" angle already noted in `docs/FULL_WAVE_SCOPE.md`: the two
literature checks above found no A-Phi-native ABC or PML, and there's a
reasoned explanation for that gap -- part evidenced, part this project's own
inference, kept clearly separated below.

**Evidenced: two communities with non-overlapping needs.** Every source
found points the same direction. A-Phi's home has historically been the
magnetoquasistatic/eddy-current world (motors, transformers, power
electronics, geophysics), where fields don't radiate by construction, so
there was never a reason to develop a radiation boundary for potentials.
Meanwhile, every ABC/PML result found stays in E/H fields even when the
paper's whole point is DC-to-high-frequency validity: Chen & Chew (2017)
implements ABC on E/H inside a DEC framework; Zhu & Jiao's "theoretically
rigorous full-wave FEM from DC to high frequencies" (IEEE Trans. Adv.
Packaging, vol. 33, no. 3, 2010) is *also* E-field, and even it doesn't
implement a true open/radiating boundary -- it truncates a packaging geometry
with a Neumann-type condition on an air layer, because its target problems
are closed packaging structures, not radiating ones. So the DC-to-high-
frequency crowd and the radiating-antenna crowd are two different
literatures, and A-Phi historically sits with neither: it wasn't built by
the packaging crowd (who stayed in E-field and didn't need a radiation
boundary anyway) or the antenna crowd (who never needed DC stability).
Ansys's own product split -- Maxwell (A-Phi, magnetoquasistatic, no
radiation) versus HFSS (E-field, full-wave, proper radiation boundaries) --
is the commercial expression of the same non-overlap: nobody who had A-Phi
also needed radiation, and nobody who needed radiation had a reason to touch
potentials.

**This project's own inference, not sourced to a citation:** two secondary
factors plausibly compound the gap, offered here as reasoning rather than a
literature finding.

1. *Gauge non-uniqueness makes a "generic" A-Phi boundary condition less
   publishable than the E-field version.* E is physical and gauge-invariant,
   so one Silver-Muller formula covers every case. **A** alone is not
   physical -- it only means something once a gauge is fixed. The condition
   derived above assumes no gauge (it uses only the potential definitions),
   but a researcher wanting a single, general-purpose closed-form result
   would need to either commit to a specific gauge or carry it as a
   parameter, which is a less clean, less citable contribution than the
   field-based version.
2. *The derivation is arguably too mechanical, on its own, to be
   paper-worthy.* What this document does -- substitute known potential
   definitions into an already-published equation -- is a page of algebra
   applied to existing results, not new physics. Anyone who needed this
   before likely did the same substitution privately (as this project just
   did) without publishing it separately, since it doesn't read as a novel
   enough contribution on its own to submit anywhere.

**Ties back to this project's own novelty case.** The same underlying cause
already identified in `docs/REFERENCES.md`'s novelty note for the BEM-coupling
question applies here too: nobody needed an all-frequency-stable A-Phi system
with a real open/radiating boundary until a problem class required both ends
of the spectrum in the *same* solve. Sub-micron THz phased-array feeds
co-designed with electrically large radiating apertures is exactly that kind
of problem. The ABC gap and the BEM-coupling gap are plausibly the same gap,
found twice, from two different angles -- worth keeping in mind together when
scoping the Phase 12 standalone-paper angle, rather than treating them as two
unrelated findings.

## Derivation

### Starting point: the standard field-based first-order ABC

The standard first-order (Silver-Muller-type) absorbing boundary condition,
applied on a truncation surface Gamma_inf enclosing all sources/scatterers, is

```
n x (curl E) - j k [n x (n x E)] = 0      on Gamma_inf
```

with `n` the outward unit normal, `k = omega * sqrt(mu * eps)` the free-space
wavenumber, and the `e^{+j omega t}` time convention used throughout this
project's equations (matching `docs/ROADMAP.md` and `docs/FULL_WAVE_SCOPE.md`).
This condition is exact for an outgoing spherical wave in the far field and
approximate at finite radius; its accuracy degrades for non-normal incidence,
which is the ABC's known, bounded limitation referenced above.

**Verified against a specific peer-reviewed source (Sept 2026 update):** this
exact 3-D first-order condition, `n x (curl E) = j*omega*mu0*(n x H) ~=
-j*k0*[n x (n x E)]`, appears as Eq. (83) in Chen & Chew, "Numerical
Electromagnetic Frequency Domain Analysis with Discrete Exterior Calculus,"
J. Comput. Phys., vol. 350, pp. 668-689, 2017 (arXiv:1704.05145) -- quoted
there in the opposite (physics, `e^{-i omega t}`) time convention, which is
why their sign on the `k` term is `+i` where this document's `e^{+j omega t}`
convention gives `-j`; the physical content is identical. That paper (by W.
C. Chew, the same author behind the generalized Lorenz gauge chosen in Phase
03) implements exactly this ABC in a differential-forms/DEC setting, but
formulated on **E** and **H**, not on **A** and Phi -- confirming the gap
this document's derivation fills. It explicitly does not implement PML; see
the PML literature check below.

### Substituting the potential definitions

This project already defines the fields in terms of the potentials as

```
E = -j*omega*A - grad(Phi)
B = curl A          =>   curl E = curl(-j*omega*A - grad(Phi)) = -j*omega*(curl A)
```

(the `curl(grad Phi) = 0` identity removes Phi from `curl E` directly, which
is exactly why a field-based ABC's Phi-dependence has been hiding inside E
all along, not because Phi has no boundary role here.)

Substitute both into the field-based ABC:

```
n x [-j*omega*(curl A)] - j*k*{ n x [n x (-j*omega*A - grad(Phi))] } = 0
```

Expand the second bracket:

```
n x [n x (-j*omega*A - grad(Phi))] = -j*omega*[n x (n x A)] - [n x (n x grad(Phi))]
```

so the equation becomes

```
-j*omega*[n x (curl A)] - j*k*{ -j*omega*[n x (n x A)] - [n x (n x grad(Phi))] } = 0
```

```
-j*omega*[n x (curl A)] - k*omega*[n x (n x A)] + j*k*[n x (n x grad(Phi))] = 0
```

Divide through by `-j*omega` (nonzero for omega != 0; see the DC caveat below):

```
n x (curl A) - j*k*[n x (n x A)] - (k/omega)*[n x (n x grad(Phi))] = 0
```

### Result

```
+-----------------------------------------------------------------------+
|  n x (curl A) - j*k*[n x (n x A)] - (k/omega)*[n x (n x grad Phi)] = 0 |
|                                                on Gamma_inf            |
|              k = omega * sqrt(mu * eps)                                |
+-----------------------------------------------------------------------+
```

This is a genuinely coupled boundary condition: even though the field-based
ABC has no explicit Phi term (it was absorbed into E), decomposing back into
the two potential unknowns produces a condition that mixes the tangential
curl of **A**, the tangential component of **A** itself, and the tangential
component of grad(Phi). All three terms need to be assembled at the
truncation boundary in Phase 07 -- dropping the grad(Phi) term would silently
reduce this back to an incomplete condition.

### Weak-form boundary term

For assembly, this enters the weak form the same way any Robin-type natural
boundary condition does: test against a vector test function `w` in the
**A**-space and integrate over Gamma_inf,

```
integral_{Gamma_inf} w . [n x (curl A)] dS
   = integral_{Gamma_inf} w . { j*k*[n x (n x A)] + (k/omega)*[n x (n x grad Phi)] } dS
```

with the left-hand side being exactly the natural boundary term that already
appears from integrating the curl-curl term by parts in Phase 04's weak form
-- so this ABC replaces the zero-Neumann (or PEC) boundary term with this
Robin-type expression, rather than requiring a new integration-by-parts step.

### Refined assembly: a symmetric surface-gradient coupling block (Sept 2026 update)

An independent review of this derivation (requested from Google Gemini, then
checked against this document before being accepted -- same practice already
used for the earlier Gemini review noted in `docs/REFERENCES.md`'s terahertz
phased-array section) identified a cleaner way to assemble the third,
Phi-coupling term above than leaving it as a raw `grad Phi` surface integral.

The tangential double-cross-product of a surface gradient reduces to the
negative surface (tangential) gradient itself:

```
n x (n x grad Phi) = -grad_Gamma(Phi)
```

where `grad_Gamma` denotes the gradient restricted to the tangent plane of
Gamma_inf. Substituting this into the right-hand side's third term and
integrating that piece by parts over the closed surface Gamma_inf (which has
no boundary of its own, so no extra edge term appears) moves the derivative
off Phi and onto the test function:

```
integral_{Gamma_inf} w . [ -(k/omega) * grad_Gamma(Phi) ] dS
   = (k/omega) * integral_{Gamma_inf} [ div_Gamma(w) ] * Phi  dS
```

with `div_Gamma` the surface divergence of the tangential test function `w`.

This is a genuine implementation improvement, not just an algebraic
rearrangement: it turns the Phi-coupling term into a symmetric bilinear block
between the discrete surface divergence of the Whitney (curl-conforming) edge
elements already used for **A** and the nodal (H1) basis functions already
used for Phi, evaluated with the *same* Gamma_inf surface mesh and the *same*
quadrature rule used for the rest of the boundary term -- rather than
assembling a raw tangential-gradient-of-Phi term against a vector test
function, which needs care to keep consistent with how grad(Phi) is
discretized elsewhere in the volume. Because both sides reduce to the same
kind of edge-to-node coupling integral, this keeps the boundary matrix
exactly consistent with the interior variational formulation and its
existing edge/nodal basis pairing, with no new element type or extra
consistency check required at the truncation surface.

The full boxed result from the earlier section is unchanged by this -- this
is a weak-form assembly refinement of the third term only, applied after
testing and integrating, not a change to the strong-form boundary condition
itself.

### DC caveat

Dividing by `-j*omega` above means this condition is undefined at omega = 0,
exactly like Chew's generalized Lorenz gauge chosen for Phase 03 (see
`docs/REFERENCES.md`, "Gauge for mixed conductor/dielectric ports"). This is
expected and consistent with the rest of the project: Phase 01, step 3
already commits to a dedicated magnetostatic solve at omega = 0 rather than a
naive limit of the frequency-domain system, and an open-boundary radiation
condition is not physically meaningful at DC anyway (there is no radiating
wave to absorb).

## Implementation note

Document, per the roadmap's own requirement, the standoff distance from the
scatterer/source region to Gamma_inf used in the Phase 07 test problem, and
report the observed error against the analytical solution as a function of
incidence angle -- this is what turns "an ABC has angle-dependent accuracy"
from a caveat into an actual, citable number for this project's own
benchmark suite (Phase 08).

## References

- S. Chen, W. C. Chew, "Numerical Electromagnetic Frequency Domain Analysis
  with Discrete Exterior Calculus," J. Comput. Phys., vol. 350, pp. 668-689,
  2017 (arXiv:1704.05145) -- Eq. (83) is the specific, verified 3-D
  first-order Silver-Muller ABC (on E, H) this derivation starts from.
  Implements PEC/PMC/Dirichlet/periodic/ABC; explicitly does not implement
  PML.
- A. F. Peterson, S. L. Ray, R. Mittra, *Computational Methods for
  Electromagnetics*, IEEE Press -- general Silver-Muller/Sommerfeld ABC
  background (already cited in `docs/REFERENCES.md` for EFIE/MFIE/CFIE
  background).
- J. Rabina, "On a numerical solution of the Maxwell equations by discrete
  exterior calculus," Ph.D. dissertation, University of Jyvaskyla, 2014 --
  the PML-in-DEC reference that Zhang, Na, Jiao & Chew (2022) defer to; not
  directly verified (dissertation not accessible), and in any case part of
  the same E/H-based DEC lineage, not A-Phi.
- B. Zhang, D.-Y. Na, D. Jiao, W. C. Chew, "An A-Phi Formulation Solver in
  Electromagnetics," arXiv:2207.02260, 2022 -- literature-check reference;
  uses IBC/PEC truncation, not an A-Phi-native ABC.
- G. Ciuprina, R. V. Sabriego, "Electric circuit element boundary conditions
  for electromagneto-quasistatic and full wave models in A, phi potentials
  and their finite element implementation," J. Math. Industry, vol. 14,
  art. 27, 2024 -- literature-check reference; port/circuit BCs only.
- A. Chervyakov, "On the use of mixed potential formulation for
  finite-element analysis of large-scale magnetization problems with large
  memory demand," arXiv:2307.12308, 2023 -- literature-check reference;
  magnetostatic PMC/insulation truncation only.
- The E = -dA/dt - grad(Phi) and B = curl A potential definitions used above
  are the same ones already in use throughout this project's formulation
  (see `docs/ROADMAP.md`, Phase 01).
- J. Zhu, D. Jiao, "A Theoretically Rigorous Full-Wave Finite-Element-Based
  Solution of Maxwell's Equations from dc to High Frequencies," IEEE Trans.
  Adv. Packaging, vol. 33, no. 3, 2010 -- "why this gap exists" reference;
  E-field (not A-Phi), DC-to-high-frequency, and truncates a packaging
  geometry with a Neumann-type condition rather than a true radiation
  boundary -- evidence that even the DC-to-high-frequency literature has
  stayed in closed (non-radiating) problem classes.
- Independent review conducted with Google Gemini,
  `Gemini/Open_Boundary_ABC_vs_PML_Evaluation.tex` (and compiled PDF) --
  a secondary synthesis of this document's own derivation and the public
  literature above, not an independent source; checked against this
  document's derivation (Sept 2026) before folding in the surface-gradient
  weak-form refinement above, per this project's standard practice for
  Gemini-sourced reviews (see `docs/REFERENCES.md`, terahertz phased array
  scope note).
