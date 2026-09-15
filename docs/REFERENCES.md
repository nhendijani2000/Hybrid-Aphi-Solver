# Reference literature

Public, published sources this project's formulation and numerics are built from. Keep this list
current as new papers are added to the project's literature folder.

## A-Phi formulation, general

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

## Tree-cotree gauge

- S.-C. Lee, J.-F. Lee, R. Lee, "Hierarchical Vector Finite Elements for Analyzing Waveguiding
  Structures," IEEE Trans. Microwave Theory Tech., vol. 51, no. 8, 2003.
- I. Munteanu, "Tree-cotree condensation properties" (technical article) -- projection-based
  interpretation of tree-cotree gauging techniques (Albanese-Rubinacci, Manges, Ticar, Munteanu
  sym/unsym) with condition-number comparison.

## Coulomb gauge

- Y.-L. Li, S. Sun, Q. I. Dai, W. C. Chew, "Vectorial Solution to Double Curl Equation With
  Generalized Coulomb Gauge for Magnetostatic Problems," IEEE Trans. Magn., vol. 51, no. 8, 2015.
- S. Yan, "Continuous Discontinuous Galerkin Method for Electromagnetic Simulations Based on an
  All Frequency Stable Formulation," Progress In Electromagnetics Research M, vol. 106, 2021.
- (folder) AphiMstaticColoumbgauge_yanpu2017, AphiTdomainColoumbgauge_yanpu2017,
  AphiFreqDomainCoulombgauge_yanpu2017, AphifreqdomainGauge_guanghuaJFlee1999 -- add full
  citations here as each is reviewed in detail.

## Matrix conditioning

- J. Jin et al., "Improving matrix conditioning" (tree-cotree context), 2008 (add full citation).

---

**Note:** this project intentionally does not draw on any employer-authored internal material
(slide decks, internal reports, or internal implementation notes). Only content published in
peer-reviewed venues, standard textbooks, or otherwise public sources goes into the design.
