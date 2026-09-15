# Scope decision: full-wave A-Phi now, FEM-BI hybrid deferred

## What changed

The original plan specialized the A-Phi formulation to low-frequency (dropped the
epsilon / displacement-current terms) from Phase 01 onward. That's been revised:
the project now keeps the full form of the coupled system --

```
curl(1/mu curl A) + (j*omega*sigma - omega^2*eps) A + (sigma + j*omega*eps) grad(Phi) = J_imp
-div[(j*omega*sigma - omega^2*eps) A + (sigma + j*omega*eps) grad(Phi)] = 0
```

-- rather than dropping the eps terms, so the same potential system is valid from
DC through full-wave/radiating regimes. This is the "all-frequency stable"
formulation described in Zhao & Fu (2017) and Yan (2021), both already in
`APhi_Papers/` and cited in `docs/REFERENCES.md`. The gauge work (tree-cotree now,
Coulomb later) and the DOF layout carry over unchanged -- this is a change to
which coefficients get zeroed out, not a different codebase.

The low-frequency-reduced system (drop epsilon, confine Phi to the conductive
region) stays available as a cheaper optional path for pure eddy-current/EDA
problems that don't need full displacement current -- a runtime/config choice
sharing the same operator blocks, not a fork.

## Why this doesn't require abandoning A-Phi

The original motivation for reconsidering the formulation was a real problem:
a zero-tangential-field (PEC) truncation boundary can't represent radiation to
infinity -- a dipole or scattering problem needs fields that genuinely extend to
infinity, not vanish at an artificial wall. But that's a *boundary-condition*
problem, not a *field-variable* problem -- an E-field full-wave FEM with the same
PEC wall would have the identical failure. So the fix belongs in the open-boundary
treatment (Phase 07: ABC/PML now, Phase 12: exact FEM-BI coupling later), not in
abandoning the A-Phi potential formulation.

## FEM-BI hybrid: deferred, not abandoned

A hybrid finite-element / boundary-integral method (coupling the interior FEM
solve to a surface integral equation on the truncation boundary, giving an exact
rather than approximate radiation condition) is the right long-term answer for
open-region problems -- this is the well-precedented approach in the literature
(Jin & Volakis) and is what distinguishes the most advanced tier of commercial
full-wave solvers. It is also a substantial, separate subsystem: choice of
boundary formulation, singular-integral quadrature, the FEM-BI interface
coupling, and -- critical for any real problem size -- a fast method (FMM or
ACA) to avoid the dense BI matrix's O(N^2)/O(N^3) cost.

Decision: defer this to Phase 12, after Phases 00-11 produce a working, validated
FEM-only solver (full-wave A-Phi + ABC/PML open boundary). Don't let the two
undertakings force each other's timeline -- see the roadmap artifact for the full
phase breakdown.

### Updated (Sept 2026): prefer a potential-based BEM over classical EFIE/MFIE/CFIE

Classical field-based surface integral equations (EFIE/MFIE, combined as CFIE to
avoid interior-resonance artifacts) are the default choice in most FEM-BI
literature, but EFIE has its own well-known low-frequency breakdown -- a
different failure mode from the one A-Phi already solves on the FEM side. A
classical-SIE coupling can therefore still degrade at low frequency even with a
perfect interior solve.

Pairing the all-frequency-stable A-Phi FEM instead with a *potential-based* BEM
(Sharma & Triverio, arXiv:2108.02764 and arXiv:2112.07360) keeps the whole
FEM-BI hybrid uniformly stable from DC through full-wave, and couples more
naturally at the interface, since both sides are already expressed in A and Phi
rather than needing E/H reconstructed just to hand quantities across the
boundary. CFIE stays as a documented fallback if the potential-based route runs
into trouble specific to a given geometry class.

A literature check at the time of this decision found closely related but
distinct prior work -- a broadband A-Phi solver via discrete exterior calculus
with no BI coupling (Zhang, Na, Jiao & Chew, arXiv:2207.02260); the
potential-based BEM above, which is pure BEM with no FEM coupling; and older
A-V/BEM-FEM couplings limited to magnetostatic/eddy-current problems -- but not
an existing all-frequency-stable A-Phi FEM coupled to a potential-based BEM that
stays uniformly low-frequency-stable end to end. That gap, if it holds up under
a fuller citation-checked review, is also a plausible angle for a standalone
paper alongside the product work -- see `docs/REFERENCES.md`.

## References for this decision

- Y. Zhao, W. N. Fu, "A New Stable Full-Wave Maxwell Solver for All Frequencies,"
  IEEE Trans. Magn., vol. 53, no. 6, 2017.
- S. Yan, "Continuous Discontinuous Galerkin Method for Electromagnetic
  Simulations Based on an All Frequency Stable Formulation," Progress In
  Electromagnetics Research M, vol. 106, 2021.
- S. Sharma, P. Triverio, "Electromagnetic Modeling of Lossy Materials with a
  Potential-Based Boundary Element Method," arXiv:2108.02764, 2021.
- S. Sharma, P. Triverio, "Electromagnetic Modeling of Lossy Interconnects From
  DC to High Frequencies With a Potential-Based Boundary Element Formulation,"
  arXiv:2112.07360, 2021.
- B. Zhang, D.-Y. Na, D. Jiao, W. C. Chew, "An A-Phi Formulation Solver in
  Electromagnetics," arXiv:2207.02260, 2022 -- cited for contrast (no BI
  coupling).
- J.-M. Jin, J. L. Volakis, and collaborators -- hybrid finite-element /
  boundary-integral method for scattering and radiation (general method,
  standard in the computational electromagnetics literature).
- A. F. Peterson, S. L. Ray, R. Mittra, *Computational Methods for
  Electromagnetics*, IEEE Press -- general EFIE/MFIE/CFIE background.
