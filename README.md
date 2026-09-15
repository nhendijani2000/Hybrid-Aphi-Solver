# A-Phi Solver

An independent, from-scratch finite-element solver built on the magnetic-vector /
electric-scalar potential (**A-Φ**) formulation, using an all-frequency-stable
variant so it's valid from DC through full-wave, not just low-frequency
eddy-current problems (see `docs/FULL_WAVE_SCOPE.md`). Targets both
signal/power-integrity (EDA) analysis and scattering problems across that
frequency range. An exact open boundary via a FEM-boundary-integral hybrid is
planned as a deferred later phase, not part of the initial (v1) solver.

This project is built from the published literature only (see `docs/REFERENCES.md`). No proprietary
or employer-derived material is used in its design.

## Status

Early development. Implemented so far: a dependency-free dense complex-matrix
utility (`complex_matrix.hpp`), the two matrix-conditioning strategies for the
A-Φ block system (`conditioning.hpp`, see `docs/CONDITIONING.md`), and symmetric
diagonal equilibration (`equilibration.hpp`, see `docs/LINEAR_SOLVER.md`, which
also records the sparse-solver build-vs-buy decision — leaning MUMPS). Mesh,
basis functions, gauge treatment, and real assembly are not implemented yet —
everything above currently operates on hand-built block systems (see
`tests/test_conditioning.cpp`, `tests/test_equilibration.cpp`), not a real solve.

## Planned capabilities

- A-Φ finite-element formulation, all-frequency-stable (DC through full-wave); low-frequency-reduced form available as a cheaper option
- Gauge treatment: tree-cotree splitting first; generalized/implicit Coulomb gauge as a second track
- Whitney edge elements (A) + nodal elements (Φ)
- Unstructured tetrahedral meshing support
- Iterative and direct linear solvers with preconditioning tuned for low-frequency conditioning
- ABC/PML open boundary initially; exact FEM-boundary-integral hybrid as a deferred later phase (`docs/FULL_WAVE_SCOPE.md`)
- Verification against analytical benchmarks and cross-checks against commercial tools

## Build

Requires CMake 3.20+ and a C++17-capable compiler (MSVC / Visual Studio 2022 recommended on Windows).

```
cmake -S . -B build
cmake --build build
```

## License

TBD — decide before any external release or collaborator access (see roadmap, Phase 0).

## Repository layout

```
include/aphi_solver/   Public headers
src/                   Implementation
tests/                 Unit and verification tests
third_party/           Vendored or fetched dependencies (kept out of git; see .gitignore)
cmake/                 CMake helper modules
docs/                  Roadmap, references, design notes
```
