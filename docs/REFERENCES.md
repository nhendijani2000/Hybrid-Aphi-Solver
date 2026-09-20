# Reference literature

Public, published sources this project's formulation and numerics are built from. Keep this list
current as new papers are added to the project's literature folder.

## A-Phi formulation, general

- J. D. Jackson, *Classical Electrodynamics*, Wiley -- standard textbook source for Maxwell's
  equations and the general A/Phi potential-definition derivation used in `docs/FORMULATION.md`,
  Sec. 1 (general background, not specific to any one implementation).
- P. Dular et al., "Dual magnetodynamic formulations and their source fields associated with
  massive and stranded inductors," IEEE Trans. Magn., vol. 36, no. 4, 2000.
- Y. Zhao and W. N. Fu, "A New Stable Full-Wave Maxwell Solver for All Frequencies," IEEE Trans.
  Magn., vol. 53, no. 6, 2017. -- basis for the project's all-frequency-stable formulation
  decision; see `docs/FULL_WAVE_SCOPE.md`.
- S. Yan, "Continuous Discontinuous Galerkin Method for Electromagnetic Simulations Based on an
  All Frequency Stable Formulation," Progress In Electromagnetics Research M, vol. 106, 2021. --
  same role as above (also listed under Coulomb gauge below for its gauge treatment).

## Full-wave open boundary (Phase 12, deferred)

- S. Sharma, P. Triverio, "Electromagnetic Modeling of Lossy Materials with a Potential-Based
  Boundary Element Method," arXiv:2108.02764, 2021 -- full-wave, DC-to-high-frequency
  potential-based BEM (pure BEM, no FEM coupling); basis for the phase's preference for a
  potential-based BEM over classical EFIE/MFIE when coupling to the interior A-Phi FEM.
- S. Sharma, P. Triverio, "Electromagnetic Modeling of Lossy Interconnects From DC to High
  Frequencies With a Potential-Based Boundary Element Formulation," arXiv:2112.07360, 2021 --
  companion/extension of the above, targeted at lossy interconnects (directly EDA-relevant).
- B. Zhang, D.-Y. Na, D. Jiao, W. C. Chew, "An A-Phi Formulation Solver in Electromagnetics,"
  arXiv:2207.02260, 2022 -- broadband A-Phi solver via discrete exterior calculus (DEC); cited for
  contrast -- no boundary-integral coupling, domain truncation only.
- J.-M. Jin, J. L. Volakis, and collaborators -- hybrid finite-element / boundary-integral (FE-BI)
  method for scattering and radiation by complex objects (general method; add the specific paper
  once selected for detailed reference).
- A. F. Peterson, S. L. Ray, R. Mittra, *Computational Methods for Electromagnetics*, IEEE Press --
  general EFIE/MFIE/CFIE background.
- (for context only, not directly used) older A-V/BEM-FEM couplings limited to magnetostatic/eddy
  current problems, e.g. IEEE Xplore doc. 9096899, "BEM-FEM formulation based on magnetic vector
  and scalar potentials for eddy current problems," 2020.

**Novelty note (Sept 2026):** a focused literature search did not turn up an existing
all-frequency-stable A-Phi FEM coupled to a potential-based BEM (as opposed to a classical
field-based SIE) that stays uniformly low-frequency-stable end to end. Treat this as a lead worth a
fuller citation-checked review before claiming it as a novel contribution -- see
`docs/FULL_WAVE_SCOPE.md` and the roadmap's Phase 12.

## Basis functions / element order (Phase 02)

**Decision (Sept 2026, final):** after two revisions this month (recorded in full, with history,
in `docs/FORMULATION.md` Sec. 5.1), settled on a **mixed-order pairing**: first-order (lowest-order)
Whitney/Nedelec A together with second-order (P2, 10-node) Phi -- not the matched-second-order
plan this section originally recorded, and not the original matched-first-order plan either.

- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed., Wiley, 2014, Ch. 5
  "Three-Dimensional Finite Element Analysis" -- verified, direct source for the second-order
  (10-node) nodal tetrahedral shape functions used for Phi (Eq. 5.52, Fig. 5.3; chapter and figure
  read directly from the shared PDF, Sept 2026, not taken from memory). Nodal/Lagrange elements
  only -- no edge, Whitney, vector-basis, H(curl), or Nedelec material in this chapter.
- Bossavit, *Computational Electromagnetism*, 1998, Ch. 5 "Whitney Elements," Prop. 5.4 --
  verified, direct source for the first-order Whitney/Nedelec A element and the discrete complex
  identity (`grad(W^0) subset of W^1`, etc.) that grounds Phase 02's `CG = 0` exit test.
- Nastaran Hendijani, `BasisFunction.docx` (shared Sept 2026), summarizing `EdModel.cpp` from her
  university "3dedyaphi" multi-author eddy-current A-Phi research codebase (confirmed not
  Ansys-derived) -- confirms the first-order Nedelec A decision against a working prior
  implementation: edge-DOF storage, the `mEdgeSign` orientation convention, and the
  `CurlFromEdges` routine all match the textbook element exactly. The same document's separate
  10-point-per-tet reconstruction of A (vertex values from edge DOFs, mid-edge values by linear
  averaging) is a post-processing/visualization utility, not part of the Phase 02 discretization
  -- see `docs/FORMULATION.md` Sec. 5.1 for the distinction.
- R. D. Graglia, D. R. Wilton, A. F. Peterson, "Higher Order Interpolatory Vector Bases for
  Computational Electromagnetics," IEEE Trans. Antennas Propag., vol. 45, no. 3, 1997 -- candidate
  primary source for a possible *future* matched second-order tetrahedral edge element for A.
  Confirmed to exist (IEEE Xplore document 558649) via search but **not yet read in full**; no
  longer a Phase 02 dependency now that A stayed first-order.
- L. E. Garcia-Castillo et al., "Second-order Nedelec tetrahedral element for computational
  electromagnetics," Int. J. Numer. Model., vol. 13, no. 2-3, 2000 -- second candidate primary
  source for the same possible future A-side upgrade; title is an exact match for what's needed.
  Confirmed to exist (Wiley Online Library) via search but **not yet read in full**; also no
  longer a Phase 02 dependency.

**Known, accepted consequence of the mixed-order pairing (Sept 2026):** the gradient of a P2
nodal (Phi) shape function has a linear part with a symmetric Jacobian, while the first-order
Nedelec A space's linear part is restricted to skew-symmetric-Jacobian fields -- so
`grad(P2 Phi)` is not exactly representable in the first-order A edge space (a from-scratch
linear-algebra check, not sourced from any paper; see `docs/FORMULATION.md` Sec. 5.1 for the
derivation). This caps A-Phi coupling accuracy at first order, a known and accepted trade-off.
Tree-cotree gauging (Phase 03) is unaffected -- it operates on A's own edge space and the mesh
graph, independent of Phi's order -- so the classical treatment applies without qualification and
the two A-side candidates above are deferred to a possible Phase 10 upgrade, not a current
dependency. See `docs/ROADMAP.md`, Phase 02 step 3 and Phase 03's now-resolved note.

## Tree-cotree gauge

- S.-C. Lee, J.-F. Lee, R. Lee, "Hierarchical Vector Finite Elements for Analyzing Waveguiding
  Structures," IEEE Trans. Microwave Theory Tech., vol. 51, no. 8, 2003. -- Sec. II states directly
  that their A-V formulation "can be better conditioned than the one obtained by E- or H-field
  formulation" (verified against the actual paper text, Sept 2026); this is the direct source for
  the Phase 01 "why A-Phi" pitch in `docs/FORMULATION.md`, Sec. 4. **Scope correction (Sept
  2026):** this paper's Table I hierarchical basis set is a 2-D triangular-element scheme for
  waveguide cross-sections (Fig. 1), not 3-D tetrahedra -- it does not source the A-side
  second-order element decision above (see the Basis functions section) and is retained here only
  for the A-V-vs-E/H conditioning claim, which is unaffected by the correction.
- I. Munteanu, "Tree-cotree condensation properties" (COMPUMAG technical article, PDF at
  compumag.org/jsite/images/stories/newsletter/ICS-02-09-1-Munteanu.pdf) -- projection-based
  interpretation of tree-cotree gauging techniques (Albanese-Rubinacci, Manges, Ticar, Munteanu
  sym/unsym) with condition-number comparison; fetched and checked directly (Sept 2026) before
  implementing Phase 03 steps 2-3 (`src/gauge_variants.cpp`). Verified with high confidence,
  directly quoted: the general oblique/orthogonal projection framework ("W^T A V y = W^T f" for
  test/trial subspace bases W, V); the edge vector **a** = [a_c; a_t] after cotree-then-tree
  reordering; Method A (Albanese-Rubinacci) sets a_t = 0 exactly and reduces to M's cotree
  principal submatrix; the essential incidence matrix F = G_c * G_t^{-1} and a_t = -F^T a_c;
  Method D (Munteanu unsymmetric) uses trial basis V = L^T = [I | -F]^T with the plain
  cotree-row test basis W; and the reported ordering kappa_D < kappa_C < kappa_A < kappa_E <
  kappa_B. **One caveat, disclosed rather than silently resolved:** the specific closed-form
  reduced-matrix equations (her eqs. 19-20) came back from the PDF extraction with what looks
  like a symbol collision -- "F" used both for "number of mesh faces" (an earlier section) and
  for "the essential incidence matrix" (this section) -- so copying that equation text verbatim
  was judged too risky. `build_munteanu_unsymmetric_gauge` instead derives the reduced matrix
  from the paper's own general projection framework plus her explicitly stated choice of test/
  trial bases for this method, and cross-checks the result two ways: the resulting matrix is
  generically non-symmetric (matching the variant's own name), and on both of this project's test
  meshes the derived Method D gives a *lower* condition-number estimate than Method A, matching
  the direction of her reported kappa_D < kappa_A ordering even at this tiny test-mesh scale
  (see `tests/test_gauge_variants.cpp`). **The full independent derivation itself** (not just
  this summary of what was and wasn't trusted from the paper) is written up step-by-step in
  `docs/TREE_COTREE_GAUGE.md` -- it re-derives `G_t` invertibility, `F = G_c*G_t^{-1}`, the
  `G^T a = 0 <=> a_t = -F^T a_c` equivalence, a from-scratch transversality proof that both
  Method A's and Method D's gauge choices are valid, and the block-matrix reduction for both
  reduced systems, all without depending on trusting the paper's specific OCR'd equations.
- F. Rapetti, A. Alonso Rodriguez, E. De los Santos, "On the Tree Gauge in Magnetostatics," J, vol.
  5, no. 1, 2022 -- clarifies that the tree gauge is not a discretization of the Coulomb (or any
  orthogonality) condition; relevant to why tree-cotree has no variational fallback at a mixed
  conductor/dielectric port boundary (Phase 03).

## Coulomb gauge

- Y.-L. Li, S. Sun, Q. I. Dai, W. C. Chew, "Vectorial Solution to Double Curl Equation With
  Generalized Coulomb Gauge for Magnetostatic Problems," IEEE Trans. Magn., vol. 51, no. 8, 2015.
- S. Yan, "Continuous Discontinuous Galerkin Method for Electromagnetic Simulations Based on an
  All Frequency Stable Formulation," Progress In Electromagnetics Research M, vol. 106, 2021.
- (folder) AphiMstaticColoumbgauge_yanpu2017, AphiTdomainColoumbgauge_yanpu2017,
  AphiFreqDomainCoulombgauge_yanpu2017, AphifreqdomainGauge_guanghuaJFlee1999 -- add full
  citations here as each is reviewed in detail.
- N. A. Demerdash, R. Wang, "Theoretical and Numerical Difficulties in 3-D Vector Potential Methods
  in Finite Element Magnetostatic Computations," IEEE Trans. Magn., vol. 26, no. 5, 1990 --
  documents Coulomb-gauge breakdown in regions mixing high- and low-permeability materials.

## Gauge for mixed conductor/dielectric ports (Phase 03, benchmark candidates)

- S. M. Ansari, C. G. Farquharson, S. P. MacLachlan, "A Gauged Finite-Element Potential Formulation
  for Accurate Inductive and Galvanic Modelling of 3-D Electromagnetic Problems," Geophysical
  Journal International, vol. 210, no. 1, 2017 -- explicitly-constrained (Lagrange-multiplier)
  Coulomb gauge for combined galvanic (contact/conduction) + inductive problems; shows an
  implicitly- or softly-enforced Coulomb gauge can give non-unique potentials at material
  interfaces where the normal component of A is discontinuous -- the likely mechanism behind
  Coulomb-gauge inaccuracy at a conductor/dielectric port.
- W. C. Chew, "Vector Potential Electromagnetic Theory with Generalized Gauge for Inhomogeneous
  Anisotropic Media," arXiv:1406.4780, 2014; companion formulation paper, Progress In
  Electromagnetics Research, 2014 -- generalized Lorenz gauge (nabla.(eps*A) = -chi*dPhi/dt, chi =
  alpha*eps^2*mu) built for spatially-varying eps/sigma, with conductor/dielectric interface
  conditions derived as part of the gauge itself; a candidate to benchmark against the explicit
  Coulomb gauge above for the mixed-port case.

**Note (Sept 2026):** these two are candidates to benchmark against tree-cotree and the existing
Phase 09 penalty-Coulomb gauge on a synthetic mixed conductor/dielectric port geometry (Phase 03) --
not yet a confirmed fix. See `docs/ROADMAP.md`, Phase 03.

**Update (Sept 2026): Chew's generalized Lorenz gauge is prioritized for implementation first.** It
adds no new unknowns (no Lagrange-multiplier field, no saddle-point solve or inf-sup pairing to get
right, unlike Ansari's explicit-constraint method) and its frequency-dependent gauge condition fits
a full-wave/radiating solver more naturally than a method built for galvanic/inductive geophysical
problems. This is a prioritization for implementation order, not a decision made from the
literature alone -- Ansari, Farquharson & MacLachlan's method remains the required benchmark
comparison on the same synthetic mixed-port geometry before either is called the default, per the
note above. Like the classical Lorenz gauge, Chew's generalized version is frequency-dependent and
should be expected to degenerate at omega = 0 the same way; the dedicated DC/magnetostatic solve
already planned (Phase 01, step 3) covers that case regardless of which AC gauge is chosen. See
`docs/ROADMAP.md`, Phase 03, step 5.

## Open boundary / ABC for A-Phi (Phase 07)

**Literature check (Sept 2026):** searched specifically for an absorbing or radiation boundary
condition written directly in terms of the A and Phi potentials (as opposed to derived E/H fields).
None was found. What exists instead:

- S. Chen, W. C. Chew, "Numerical Electromagnetic Frequency Domain Analysis with Discrete Exterior
  Calculus," J. Comput. Phys., vol. 350, pp. 668-689, 2017 (arXiv:1704.05145) -- implements PEC,
  PMC, Dirichlet, periodic, and a first-order/second-order ABC (Eq. 83 for the 3-D first-order
  case, the equation this project's ABC derivation starts from), all on E/H fields, not A-Phi.
  Explicitly does not implement PML.
- B. Zhang, D.-Y. Na, D. Jiao, W. C. Chew, "An A-Phi Formulation Solver in Electromagnetics,"
  arXiv:2207.02260, 2022 (already listed above) -- the closest existing broadband A-Phi solver.
  Its own numerical examples use a simple impedance boundary condition (IBC) "as a simple absorbing
  boundary condition" for a rod antenna case, and plain PEC truncation "for simplicity" for a
  nano-scale rod antenna case; PML is discussed in the text but its implementation is deferred to
  J. Rabina's 2014 Ph.D. dissertation (Univ. of Jyvaskyla, not directly verified), with no explicit
  A-Phi equations given anywhere in this citation chain.
- G. Ciuprina, R. V. Sabriego, "Electric circuit element boundary conditions for
  electromagneto-quasistatic and full wave models in A, phi potentials and their finite element
  implementation," Journal of Mathematics in Industry, vol. 14, article 27, 2024 -- covers port/circuit
  (ECE1-ECE5) boundary conditions only; does not address open/radiating boundaries.
- A. Chervyakov, "On the use of mixed potential formulation for finite-element analysis of
  large-scale magnetization problems with large memory demand," arXiv:2307.12308, 2023 --
  magnetostatic only; boundary treatment is PMC/magnetic-insulation domain truncation, not a
  radiation condition.

**A separate, dedicated PML-for-A-Phi search (Sept 2026) also found nothing published.** The
DEC/E-H lineage above (Chen & Chew 2017; Zhang, Na, Jiao & Chew 2022) is the closest anyone gets,
and it implements PML (per Rabina's dissertation, unverified) only on E/H fields, never on A-Phi.
General "vector potential + PML" searches return eddy-current/magnetostatic vector-potential-only
formulations (A-V, T-Omega; no coupled scalar potential, no radiation) -- a different problem
class, not evidence of an A-Phi-plus-PML radiation treatment. This reinforces rather than
undercuts the Phase 07 ABC choice above.

**Decision:** since no A-Phi-native ABC (or PML) exists in the literature, Phase 07's first-order
ABC is a from-scratch derivation, not a citation -- substituting E = -dA/dt - grad(Phi) into the
verified Chen & Chew (2017) Eq. (83) first-order Silver-Muller/Sommerfeld condition (also standard
background in A. F. Peterson, S. L. Ray, R. Mittra, *Computational Methods for Electromagnetics*,
IEEE Press). The full derivation, including the algebraic verification of the starting equation
against a specific peer-reviewed source, is in `docs/OPEN_BOUNDARY_ABC.md`.

## Matrix conditioning

- J. Jin et al., "Improving matrix conditioning" (tree-cotree context), 2008 (add full citation).
- M. R. Hestenes and E. Stiefel, "Methods of Conjugate Gradients for Solving
  Linear Systems," J. Research Nat. Bur. Standards, vol. 49, 1952 -- the
  Conjugate Gradient method itself, used (Sept 2026) as the inner solve of
  `estimate_condition_number`'s inverse power iteration (`cg_solve_ata` in
  `src/gauge_variants.cpp`), replacing a dense Gaussian-elimination inner
  solve so the estimator scales to production mesh sizes. Standard,
  well-established numerical linear algebra, not specific to this project's
  EM formulation; cited here rather than derived, unlike the EM-specific
  material elsewhere in this document. G. H. Golub and C. F. Van Loan,
  "Matrix Computations" (any edition) is the standard reference for the
  surrounding technique -- inverse power iteration tolerating an inexact
  (CG-based, not-fully-converged) inner solve at each outer step.

## Terahertz phased array scope (Sept 2026)

No new literature enters the project from this scope addition -- every technical
claim in `docs/THZ_PHASED_ARRAY_SCOPE.md` (Floquet boundary conditions,
multi-port active impedance, the ACA-over-MLFMA choice for Phase 12) traces
back to papers already listed above (Zhao & Fu 2017; Yan 2021; Sharma &
Triverio 2021, 2022; Lee, Lee & Lee 2003; Ansari, Farquharson & MacLachlan
2017; Chew, generalized Lorenz gauge). The independent review conducted with
Google Gemini (`Gemini/APhi_Solver_Proposal_Gemini.pdf`,
`Gemini/THz_Phased_Array_APhi_Evaluation.md`) is a secondary synthesis of this
same public literature, not an independent source, and was checked against it
before folding into the roadmap.

---

**Note:** this project intentionally does not draw on any employer-authored internal material
(slide decks, internal reports, or internal implementation notes). Only content published in
peer-reviewed venues, standard textbooks, or otherwise public sources goes into the design.
