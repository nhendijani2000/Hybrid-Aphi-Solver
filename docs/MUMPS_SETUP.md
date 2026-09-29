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

**A custom install is worth the extra minute.** The full toolkit is far more
than this project needs — it carries SYCL, VTune, Advisor, the DPC++ compiler
and more. In the installer choose **Custom** and select only:

- Intel Fortran Compiler
- Intel oneAPI Math Kernel Library

That keeps the install down substantially and avoids the Visual Studio
integration for components you will never use.

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
