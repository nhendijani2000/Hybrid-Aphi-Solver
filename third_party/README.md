# third_party

The scripts that build the MUMPS backend. **They are ours and are tracked; what
they build is not.**

| | tracked | why |
|---|---|---|
| `build_mumps.bat` | yes | sequential MUMPS, no OpenMP |
| `build_mumps_omp.bat` | yes | + MUMPS OpenMP and METIS |
| `build_mumps_omp_mkl.bat` | yes | **+ threaded MKL — use this one** |
| `mumps/` | no | upstream checkout, 572 KB |
| `mumps-build*/`, `mumps-install*/` | no | generated, ~190 MB |

`.gitignore` ignores `third_party/*` and then re-admits `*.bat` and this file.
**A new script here needs its own negation or it is silently untracked** — which
is how the MUMPS build came to be recorded nowhere for as long as it was.

## The source

```
git clone https://github.com/scivision/mumps.git third_party/mumps
git -C third_party/mumps checkout v5.9.1.2
```

`v5.9.1.2` is the tag these measurements were made against, and it is pinned on
purpose: `MUMPS_SETUP.md` quotes backend timings and backward errors to four
significant figures, and a floating `master` makes those unreproducible. The
repository is scivision's CMake wrapper, which fetches the MUMPS sources itself
— the upstream MUMPS release is a tarball with a hand-edited Makefile.inc, which
is not something to build on Windows.

## Building

Each script locates the source and its build tree through **`MUMPS_TP`**, which
defaults to this directory. Nothing has to move: point it at an existing tree.

```bat
set MUMPS_TP=C:\Research\APhi_Solver_Project_LowFrequency_EDA\third_party
third_party\build_mumps_omp_mkl.bat
```

Then configure the solver against what it installed:

```bat
cmake -S . -B build-mumps-mkl -G Ninja -DCMAKE_BUILD_TYPE=Release ^
  -DAPHI_WITH_MUMPS=ON ^
  -DMUMPS_ROOT=%MUMPS_TP%/mumps-install-omp-mkl ^
  -DMUMPS_DIR=%MUMPS_TP%/mumps-install-omp-mkl/cmake ^
  -DMETIS_LIBRARY=%MUMPS_TP%/mumps-install-omp-mkl/lib/metis.lib ^
  -DMETIS_INCLUDE_DIR=%MUMPS_TP%/mumps-install-omp-mkl/include ^
  -DMKL_THREADING=intel_thread
```

`run_case.bat` then picks `build-mumps-mkl` up automatically for any input with
`backend = mumps`, preferring it over `build-mumps-omp` and `build-mumps`.

## `-DMKL_THREADING=intel_thread` is not optional, and it goes HERE

It is worth **10x** — 266.78 s to 25.94 s to factor 418318 unknowns, same answer
to every printed digit. Three things make it easy to get wrong, and all three
were gotten wrong before this was written down:

1. **`MUMPS_openmp=ON` does not do it.** That is MUMPS's own tree parallelism.
   The BLAS is chosen separately, and MUMPS's `FindLAPACK.cmake` defaults
   `MKL_THREADING` to `sequential` unless `OpenMP` is in
   `LAPACK_FIND_COMPONENTS`, which `MUMPS_openmp=ON` does not put there. No
   warning is printed.
2. **Setting it on the MUMPS build alone is not enough.** `MUMPSConfig.cmake`
   re-runs `find_dependency(LAPACK COMPONENTS MKL)` at the *consumer's* configure
   time, where `MKL_THREADING` is undefined again, so the link line takes
   `mkl_sequential` however MUMPS was built. It must be passed to the APhi
   configure, as above.
3. **The CMake cache lies.** The build that still linked sequential had
   `MKL_THREADING:UNINITIALIZED=intel_thread` sitting in `CMakeCache.txt`.

So verify the binary, not the cache:

```bat
dumpbin /dependents build-mumps-mkl\solve_mesh.exe | findstr /i mkl
```

must print `mkl_intel_thread.3.dll`. If it prints `mkl_sequential.3.dll`, MUMPS
is running its dense frontal kernels on one core — about 4 of 24 busy — and
nearly all the flops of a large factorization are in those kernels.

## A note on oneAPI

oneAPI's own `setvars.bat` is **broken** on this install: it reports
`'vars.bat' is not recognized` for every component, then claims success, leaving
`MKLROOT` empty and `ifx` off `PATH`. The per-component scripts work, so the
build scripts call them directly:

```bat
call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat"
call "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\env\vars.bat"
```

The same applies when running a MUMPS-linked binary by hand — it needs
`libifcoremd.dll` and will not start without that environment. `run_case.bat`
loads it for you.

See `docs/MUMPS_SETUP.md` for the full setup and the measured timings.
