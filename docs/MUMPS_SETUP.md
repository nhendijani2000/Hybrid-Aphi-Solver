# Setting up the MUMPS backend on Windows

What has to exist on the machine before `-DAPHI_WITH_MUMPS=ON` can work, and in
what order. Nothing here is required to build or run this project — the
internal solver is the default and needs none of it. See `SOLVER_PLAN.md` §12
for why the option exists at all.

Checked against Intel's download page on 2026-09-29.

---

## Step 1 — Intel oneAPI Toolkit

**The HPC Toolkit no longer exists as a separate download.** As of the 2026.0
release Intel merged the Base Toolkit and the HPC Toolkit into one *Intel®
oneAPI Toolkit*. Looking for "HPC Toolkit" will send you to an archive page for
an old version; don't.

One download covers both things MUMPS needs:

| component | what it is for |
|---|---|
| Intel® Fortran Compiler (`ifx`) | MUMPS is Fortran; MSVC cannot build it |
| Intel® oneAPI Math Kernel Library (oneMKL) | the BLAS and LAPACK MUMPS links against |

```
https://www.intel.com/content/www/us/en/developer/tools/oneapi/oneapi-toolkit-download.html
```

- Version **2026.1.1**, **1824 MB**, released 2026-09-08.
- Choose **Windows**, then the **offline installer**.
- Use **"Continue as a Guest (download starts immediately)"** — no Intel account
  is needed.

**Prerequisite:** Visual Studio must already be installed; the oneAPI installer
checks for it and integrates with it. This machine has VS 18 Community, which
is what `build.bat` already uses, so that is satisfied.

**Take the Recommended installation.** An earlier draft of this file advised a
custom install to trim the toolkit; that advice was written before seeing the
installer's own numbers and is withdrawn.

    installation size   8.5 GB
    location            C:\Program Files (x86)\Intel\oneAPI

Three reasons Recommended is the better choice:

- 8.5 GB against ~800 GB free is not a saving worth optimising for.
- **VTune Profiler is bundled**, and it is the right tool for the open question
  at the bottom of this file -- where the 550 s of factorization actually goes.
  A trimmed install would have dropped it.
- The risk is asymmetric. Deselecting something the MUMPS build turns out to
  need -- a runtime, a threading layer, an OpenMP library -- surfaces as a
  confusing link error an hour into step 2. oneMKL in particular has
  threading-layer dependencies that are not obvious from the component names.

The two components that matter are both in the Recommended set: **Intel Fortran
Compiler** and **Intel oneAPI Math Kernel Library**.

**Silent install**, if you prefer the command line — Intel's parameter list is
at `Command Line Installation Parameters` on the download page:

```
bootstrapper.exe -s -a --silent --eula accept
```

Admin rights are not strictly required: the installer supports a user-scope
install. A user-scope install does sometimes leave `setvars.bat` out of the
places tooling expects, so note where it lands.

### Verifying step 1

Open a **new** terminal (the installer changes the environment; an existing
shell will not see it) and run Intel's environment script, then check both
pieces:

```
"C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
ifx --version
echo %MKLROOT%
```

`ifx` should print a version, and `MKLROOT` should be a path, not empty. If
`setvars.bat` is somewhere else, a user-scope install usually puts it under
`%USERPROFILE%\intel\oneapi\`.

---

## Step 2 — MUMPS itself

Not yet done, and deliberately not attempted before step 1 exists.

MUMPS does not ship a Windows build. The least painful route is the CMake
wrapper at `github.com/scivision/mumps`, which fetches the MUMPS sources and
builds them against a Fortran compiler and a BLAS/LAPACK of your choosing —
`ifx` and oneMKL from step 1.

What this project needs from it:

- the **complex double** interface, `zmumps` — the `RowScaled` and `ScaledPhi`
  conditionings produce complex symmetric indefinite systems, which is exactly
  what MUMPS handles natively and why it was chosen over SuperLU
- `include/zmumps_c.h` and a `zmumps` import library in `lib/`

METIS is optional. Without it MUMPS orders with AMD, which is what the internal
solver already uses — so a first comparison without METIS is the fairer one
anyway, since it isolates the factorization from the ordering.

Then:

```
cmake -S . -B build -DAPHI_WITH_MUMPS=ON -DMUMPS_ROOT=<prefix>
```

The CMake option errors out if `zmumps` is not found rather than guessing, so a
misconfigured prefix fails at configure time, not at link time.

---

## Step 3 — the driver

`src/mumps_backend.cpp` is referenced by the build when the option is ON and
**has not been written**. It was left unwritten on purpose: with no Fortran
compiler, no BLAS and no MUMPS on the machine, none of it could be compiled or
run, and untested Fortran-interop code is worse than none.

Once step 2 produces a linkable `zmumps`, the driver is the remaining work.

---

## Is this worth doing?

Measured on `02_Ansys_Cylinder_50Hz`, 240990 unknowns:

    factorization    550.29 s      96.75 % of the run
    nnz(L)           211483290
    per nonzero      ~2.6 us

2.6 µs per factor nonzero is slow enough to suggest the gap is at least partly
**our own factorization being scalar rather than blocked** — the supernodal /
BLAS3 work is Stage 2 in `SOLVER_PLAN.md` §10 and is unstarted. A blocked
factorization should be roughly an order of magnitude faster on the same
hardware.

So MUMPS may turn out to be a comparison that tells us our inner loop is the
problem, rather than a solution in itself. Profiling where those 550 seconds go
costs nothing and needs no installs, and would make this decision an informed
one. Recorded here because it is the first thing to check if MUMPS disappoints.

---

## Result: 51x on case 02

Measured 2026-09-29, both solvers single-threaded (MUMPS built without OpenMP
precisely so this comparison is like for like).

`02_Ansys_Cylinder_50Hz`, 240990 unknowns:

|  | internal | MUMPS | |
|---|---|---|---|
| factorization | 550.29 s | **10.67 s** | **51.6x** |
| total run | 568.79 s | **22.95 s** | 24.8x |
| backward error | 4.80e-21 | **1.19e-21** | |
| residual | 6.01e-13 | **1.49e-13** | |

`examples/cylinder_50hz.aphi`, 37368 unknowns: 55.25 s -> 1.60 s, **34x**.

**The answers are identical.** `tools/pv_extract_rl.py` on both output
directories gives `R = 9.796995e-05 ohm` and `L = 22.5956 nH` to every digit
printed. MUMPS is not trading accuracy for speed -- its backward error is four
times better.

**What this settles.** The question behind the whole option was whether the gap
to commercial tools is the formulation or the linear algebra. It is the linear
algebra: same mesh, same matrix, same ordering family, same answer, 51x. The
internal solver's scalar LDL^T is the bottleneck, and `SOLVER_PLAN.md` §10's
Stage 2 -- supernodal, BLAS3 -- is what would close it without a dependency.

**What has NOT been measured**: MUMPS with OpenMP, or with METIS ordering.
Both were deliberately left off for this first comparison. Either could make
the gap wider still; neither changes the conclusion.

### Running it

The built executable needs the Intel Fortran runtime (`libifcoremd.dll`), so it
must run with the compiler environment loaded:

```
call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat"
call "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\env\vars.bat"
build-mumps\solve_mesh.exe <case>.aphi
```

and the case must ask for it:

```
[solver]
conditioning = row_scaled
backend      = mumps
```

---

## OpenMP and METIS: measured 2026-09-29

A second MUMPS was built with both enabled (`third_party/build_mumps_omp.bat`,
installing to `mumps-install-omp`) so the sequential build stays available for
comparison. Case 02, 240990 unknowns, **factorization time only**:

| build | threads | factorization | vs internal |
|---|---|---|---|
| internal LDL^T | 1 | 550.29 s | — |
| MUMPS, AMD, no OpenMP | 1 | 10.67 s | 51.6x |
| MUMPS + METIS + OpenMP | 1 | 12.67 s | 43.4x |
| MUMPS + METIS + OpenMP | 8 | **6.05 s** | **91.0x** |
| MUMPS + METIS + OpenMP | 24 | 8.65 s | 63.6x |

**24 threads is worse than 8.** The machine has 24 logical cores and the best
result is at a third of them. A factorization of this size does not have enough
parallel work to feed 24 threads, and past the knee the synchronisation and
memory traffic cost more than the extra cores return. Anyone quoting a thread
count should measure it rather than assume more is better.

**METIS did not help here.** At one thread the METIS build is *slower* than the
AMD one, 12.67 s against 10.67 s. Part of that is OpenMP runtime overhead
present even at one thread, so the two are not a clean isolation of the
ordering -- but there is certainly no METIS win to collect on this problem.
Nested dissection pays off on larger, less regular systems than this.

**The answers are unchanged**: `R = 9.796995e-05 ohm`, `L = 22.5956 nH`, same
as both the internal solver and the sequential MUMPS build.

### A real inefficiency this exposed

The MUMPS path still runs **our** symbolic analysis first:

    ordering      amd   nnz(L) = 211483290   4033.72 MB   [8982.16 ms]

That is 9 seconds of AMD ordering whose result MUMPS never uses -- it does its
own analysis. It is now LARGER than the factorization it precedes. Skipping
`analyze()` when the backend is MUMPS would take case 02 from ~22 s to ~13 s.
Not yet done.
