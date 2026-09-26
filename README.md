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

Early development — **no field has been solved yet.** Everything below is
built and tested; the assembly that would join it together is the next
phase. Roughly 1000 checks across ten test executables.

Implemented:

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
- **Sparse matrices** — a complex-capable CSR type with Gustavson multiply
  (`sparse_matrix.hpp`), and the two A-Φ conditioning transforms plus
  symmetric equilibration (`conditioning.hpp`, `equilibration.hpp`).
- **Problem description and input file** — a sectioned text format, parsed
  and validated with errors that name the file and line
  (`problem.hpp`, `input_file.hpp`, `examples/`).

Not yet implemented: binding a parsed problem to a mesh, the DOF map,
element matrices, global assembly, and the linear solve. See
`docs/ROADMAP.md` Phase 04 and `docs/INPUT_FILE_PLAN.md`.

## Planned capabilities

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

Three examples, each documenting the physics it encodes:

| file | what it is |
|---|---|
| `examples/cylinder_box.aphi` | a wire in a square box at DC. Targets `R = 0.1388 mΩ` (exact against the meshed cross-section) and `L = 0.3870 nH` |
| `examples/cylinder_box_sweep.aphi` | the same geometry swept 1 kHz → 10 MHz, which spans the whole skin-effect transition |
| `examples/loop_internal_port.aphi` | a ring driven through an internal cut, showing `current_direction` as a hint. Built by `tools/loop_cut.geo` |
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
