@echo off
setlocal enabledelayedexpansion
rem ---------------------------------------------------------------------------
rem Build MUMPS (zmumps only) against Intel ifx + oneMKL, for use as the
rem optional backend of the A-Phi solver. See APhi_Solver/docs/MUMPS_SETUP.md.
rem
rem     build_mumps.bat
rem
rem Installs to third_party\mumps-install-omp, which is what MUMPS_ROOT should
rem point at when configuring the solver with -DAPHI_WITH_MUMPS=ON.
rem
rem WHAT IS TURNED OFF, AND WHY
rem   MUMPS_parallel   OFF - sequential. No MPI to install, and the comparison
rem                          against our own single-process solver is the
rem                          meaningful one.
rem   MUMPS_scalapack  OFF - requires MPI.
rem   MUMPS_openmp     OFF - deliberately, for the FIRST build. Both solvers
rem                          then run single-threaded and the factorization
rem                          times compare directly. Threading is a separate
rem                          measurement, made later and on purpose.
rem   all arithmetics  OFF except COMPLEX16 - the RowScaled and ScaledPhi
rem                          conditionings produce complex symmetric indefinite
rem                          systems, so zmumps is the only one wanted.
rem ---------------------------------------------------------------------------

rem --- where the MUMPS source and build trees live --------------------------
rem THIS SCRIPT IS TRACKED; what it builds is not. The MUMPS source is an
rem upstream checkout and the build tree it produces is ~60 MB, so neither
rem belongs in the repository -- see third_party/README.md.
rem
rem MUMPS_TP says where they are. It defaults to this script's own directory, so
rem a fresh checkout can clone MUMPS here and build in place; .gitignore already
rem excludes third_party/mumps*. Point it elsewhere when the trees already exist
rem outside the repository:
rem
rem     set MUMPS_TP=C:\Research\APhi_Solver_Project_LowFrequency_EDA\third_party
rem
rem Do NOT reintroduce %~dp0 below. These scripts used to live beside the source
rem and derived every path from their own location; moving them into the
rem repository silently retargeted cmake -S at a path that did not exist and
rem would have dropped the install prefix inside the working tree.
set "TP=%MUMPS_TP%"
if not defined TP set "TP=%~dp0"
if not "%TP:~-1%"=="\" set "TP=%TP%\"
if not exist "%TP%mumps\CMakeLists.txt" (
  echo ERROR: no MUMPS source at "%TP%mumps"
  echo        Clone it at the tag this project builds against:
  echo            git clone https://github.com/scivision/mumps.git "%TP%mumps"
  echo            git -C "%TP%mumps" checkout v5.9.1.2
  echo        or set MUMPS_TP to a directory that already has it.
  exit /b 1
)
set "PREFIX=%TP%mumps-install-omp"

rem --- MSVC first (cl, cmake, ninja), then Intel (ifx, icx, MKL) -------------
rem VS detection copied from APhi_Solver\build.bat, which already works. The
rem `(x86)` in the Program Files path breaks a naive `for /f` because the
rem closing paren terminates the FOR block, so the path is assigned first.
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSPATH="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -property installationPath`) do set "VSPATH=%%i"
)
if not defined VSPATH set "VSPATH=C:\Program Files\Microsoft Visual Studio\18\Community"
if not exist "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" (
  echo ERROR: vcvars64.bat not found under "!VSPATH!"
  exit /b 1
)
call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" >nul
if not "%ERRORLEVEL%"=="0" ( echo ERROR: vcvars64 failed & exit /b 1 )

rem oneAPI's own setvars.bat is BROKEN on this install: it reports
rem "'vars.bat' is not recognized" for every component and then claims the
rem environment was initialised, leaving MKLROOT empty and ifx off PATH. The
rem per-component scripts it is supposed to call are present and work, so they
rem are called directly. Verified 2026-09-29.
set "ONEAPI=C:\Program Files (x86)\Intel\oneAPI"
call "%ONEAPI%\compiler\latest\env\vars.bat" >nul 2>&1
call "%ONEAPI%\mkl\latest\env\vars.bat" >nul 2>&1

where ifx >nul 2>&1 || ( echo ERROR: ifx not on PATH after setvars & exit /b 1 )
if not defined MKLROOT ( echo ERROR: MKLROOT not set after setvars & exit /b 1 )
echo === toolchain
ifx --version 2>&1 | findstr /i "Fortran"
echo   MKLROOT = %MKLROOT%
echo   prefix  = %PREFIX%
echo.

rem --- NOTE ON ERROR CHECKING ------------------------------------------------
rem `if errorlevel 1` is WRONG here and was silently wrong for a long time. It
rem tests "exit code >= 1", and `cmake --build` returns -1 when Ninja fails, so
rem a broken build sailed past the check and the script printed "done" over a
rem failed compile. Compare against 0 instead. Verified 2026-09-30: a genuinely
rem failing build gives %ERRORLEVEL% = -1.
rem
rem MUMPS_BUILD_TESTING is OFF because MUMPS's own test_mumps_openmp example does not
rem link on this toolchain (unresolved externals, LNK exit 1120) while every
rem library we actually consume builds and works. With the error check fixed,
rem that example would fail the whole script for no reason. The variable is
rem MUMPS_BUILD_TESTING, guarded at CMakeLists.txt:340 as
rem ${PROJECT_NAME}_BUILD_TESTING -- the generic BUILD_TESTING does NOTHING
rem here, which cost a rebuild to discover.
rem ---------------------------------------------------------------------------
echo === configuring
cmake -S "%TP%mumps" -B "%TP%mumps-build-omp" -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_Fortran_COMPILER=ifx ^
  -DCMAKE_C_COMPILER=icx ^
  -DCMAKE_INSTALL_PREFIX="%PREFIX%" ^
  -DMUMPS_parallel=OFF ^
  -DMUMPS_scalapack=OFF ^
  -DMUMPS_openmp=ON ^
  -DMUMPS_metis=ON ^
  -DBUILD_SHARED_LIBS=OFF ^
  -DBUILD_SINGLE=OFF ^
  -DBUILD_DOUBLE=OFF ^
  -DBUILD_COMPLEX=OFF ^
  -DMUMPS_BUILD_TESTING=OFF ^
  -DBUILD_COMPLEX16=ON
if not "%ERRORLEVEL%"=="0" ( echo === CONFIGURE FAILED & exit /b 1 )

echo.
echo === building
cmake --build "%TP%mumps-build-omp" --parallel
if not "%ERRORLEVEL%"=="0" ( echo === BUILD FAILED & exit /b 1 )

echo.
echo === installing
cmake --install "%TP%mumps-build-omp"
if not "%ERRORLEVEL%"=="0" ( echo === INSTALL FAILED & exit /b 1 )

echo.
echo === done. Check for zmumps:
dir /b "%PREFIX%\lib" 2>nul
dir /b "%PREFIX%\include\zmumps_c.h" 2>nul
exit /b 0
