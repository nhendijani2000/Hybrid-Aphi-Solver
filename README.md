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

**To set up and run a case, start with [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md)** -- the input
file section by section, how to run it, and how to get pictures out.

## Status

**The solver runs end to end** — mesh in, fields out, validated against closed
forms. 1915 checks across 20 test executables, plus five regression cases under
`regression_tests/`, four of which carry a written validation report.

### The pipeline

- **Mesh ingestion** — Gmsh `.msh` 2.x and 4.1, with physical names, region
  and surface tags, derived edge/face topology and face→tet adjacency,
  validated at read time (`gmsh_reader.hpp`, `mesh.hpp`). Checked against
  meshes Gmsh itself wrote, in both formats.
- **Basis functions** — first-order Whitney/Nédélec **A** with global edge
  orientation, second-order (P2) nodal Φ (`basis_functions.hpp`; the
  mixed-order pairing is locked in, `docs/FORMULATION.md` §5.1).
- **Gauge** — boundary-first tree-cotree over `n×A = 0` surfaces, plus both
  Albanese-Rubinacci and Munteanu reductions (`tree_cotree.hpp`,
  `gauge_variants.hpp`, `docs/TREE_COTREE_GAUGE.md`).
- **Input file** — a sectioned text format, parsed and validated with errors
  that name the file and line, before the mesh is opened
  (`input_file.hpp`, `docs/USER_GUIDE.md`, `examples/`).
- **Binding and DOF map** — names resolved against the mesh, ports checked for
  planarity and placement, unknowns numbered with the gauge applied
  (`problem_binding.hpp`, `dof_map.hpp`).
- **Assembly** — element matrices over a quadrature rule, global sparsity built
  symbolically, then filled per frequency (`element_matrix.hpp`,
  `sparsity.hpp`, `assembly.hpp`).
- **Linear solve** — an in-house direct `LDLᵀ` with AMD/RCM ordering, symmetric
  equilibration and iterative refinement; **MUMPS** available as an optional,
  never-required backend (`factorization.hpp`, `docs/SOLVER_PLAN.md`,
  `docs/MUMPS_SETUP.md`).
- **Post-processing** — Φ, **A**, **B**, **H**, **E**, **J** as text and VTK,
  plus a `[postprocess]` section that selects pictures and a renderer that
  draws them (`postprocess.hpp`, `tools/postprocess.py`).

### Validated against

| case | what it tests | agreement |
|---|---|---|
| `02_Ansys_Cylinder_50Hz` | a voltage-driven wire, no skin effect | `R` to **1.3e-05** relative; Φ linear in `z` to 0.08 %; phase lag −0.628° against −0.639° |
| `04_Cylinder_SkinDepth` | skin effect at `a/δ = 3` | `R_ac/R_dc` **1.7719** against Kelvin's **1.7680**; core/surface phase lag −134.24° against −133.85° |
| `05_Loop_1A_50Hz` | a ring driven through an internal cut — the multiply-connected gauge | `E` **0.016 %**, `R` **0.13 %**, `L` **0.09 %**, `B` on the axis **0.37 %**, `B` in the air within **±0.75 %** over a decade of distance |

Case 05 is the one that exercises the gauge hardest: the current circulates with
no terminal anywhere on the boundary, and 1 A injected across the cut is
measured crossing a plane where nothing is prescribed.

### Not implemented

- **Open boundary.** Only `flux_tangential` (`n × A = 0`). `pec` and `abc` parse
  and are refused. The ABC/PML work and the FEM-BI hybrid are both still ahead
  (`docs/OPEN_BOUNDARY_ABC.md`, `docs/FULL_WAVE_SCOPE.md`). Until then the air
  box has to be large enough that the wall does not interfere — case 05's report
  measures what happens when it is not.
- **Iterative solvers and preconditioning.** Direct only.
- **Unsymmetric factorization.** The matrix is treated as symmetric throughout.
- **Supernodal factorization** — the largest known performance win, deliberately
  deferred (`docs/SOLVER_PLAN.md` §10).
- **Materials** are linear and isotropic: scalar `sigma`, `eps_r`, `mu_r` per
  body. No nonlinearity, no anisotropy.
- **The scattering track** (`docs/ROADMAP.md` 07).

See `docs/ROADMAP.md` for where each of these sits.

## Scope and direction

What is built today is in [Status](#status) above; this is the shape of the whole
thing, including the parts still ahead.

- A-Φ finite-element formulation, all-frequency-stable (DC through full-wave); low-frequency-reduced form available as a cheaper option
- Gauge treatment: tree-cotree splitting first; generalized/implicit Coulomb gauge as a second track
- Whitney edge elements (A) + nodal elements (Φ)
- Unstructured tetrahedral meshing support
- Iterative and direct linear solvers with preconditioning tuned for low-frequency conditioning
- ABC/PML open boundary initially; exact FEM-boundary-integral hybrid as a deferred later phase (`docs/FULL_WAVE_SCOPE.md`)
- Verification against analytical benchmarks and cross-checks against commercial tools

## Build and run

Requires a C++17 compiler and CMake 3.20+. On Windows both ship with Visual
Studio's "Desktop development with C++" workload — but neither is on `PATH`
in an ordinary terminal, which is what `build.bat` exists to sort out.

```
build.bat          REM configure + build into build\
run-tests.bat      REM run all ten test executables
```

Then run the solver on one of the worked examples:

```
build\aphi_solver.exe examples\cylinder_box.aphi
```

If you are already inside a *x64 Native Tools Command Prompt for VS*, or on
Linux/macOS, the scripts are unnecessary:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

Visual Studio can also open the project directly: **File → Open → CMake…**
and point it at `CMakeLists.txt`.

### What running it does today

`aphi_solver` **reads and checks, but does not solve.** It parses the input
file, opens the mesh it names, resolves every body and port against it, and
prints what it found — tet counts per body, faces and vertices per port, the
resolved port directions, and where Φ will live. What it does not have is a
DOF map, assembly or a linear solver, and its closing lines say so.

Five examples, each documenting the physics it encodes:

| file | what it is |
|---|---|
| `examples/cylinder_box.aphi` | a wire in a square box at DC. Targets `R = 0.1388 mΩ` (exact against the meshed cross-section) and `L = 0.3870 nH` |
| `examples/cylinder_box_sweep.aphi` | the same geometry swept 1 kHz → 10 MHz, which spans the whole skin-effect transition |
| `examples/loop_internal_port.aphi` | a ring driven through an internal cut, showing `current_direction` as a hint. Built by `tools/loop_cut.geo` |
| `examples/cylinder_ac.aphi` | the same wire at 100 MHz — it produces the `potential.out` quoted in `docs/POSTPROCESSING_PLAN.md` |
| `examples/loop_sweep.aphi` | the ring swept 10 kHz → 100 MHz, and where `conditioning = row_scaled` or `scaled_phi` can be tried: both are refused at DC |

Worth trying deliberately: misspell a key, delete `current_direction` from
the loop example, or turn the cylinder's 0 V port into a second current
source. Each produces a specific message rather than a silent wrong answer —
that behaviour is the point of the stage, and every case has a test.

## License

TBD — decide before any external release or collaborator access (see roadmap, Phase 0).

## Repository layout

```
include/aphi_solver/   Public headers
src/                   Implementation
tests/                 Unit and verification tests (ten executables, ~1000 checks)
examples/              Worked input files, each documenting its own physics
tools/                 .geo mesh sources, and the gauge-comparison CLI
meshes/                Test meshes -- small fixtures are tracked, see .gitignore
third_party/           Vendored or fetched dependencies (kept out of git; see .gitignore)
docs/                  Roadmap, references, design notes and the design record
build.bat              Configure + build (finds Visual Studio for you)
run-tests.bat          Run every test executable
```
