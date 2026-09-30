@echo off
setlocal enabledelayedexpansion
rem ---------------------------------------------------------------------------
rem Run one regression case end to end, from that case's own folder:
rem
rem     cd regression_tests\01_OneCylinder
rem     ..\run_case.bat cylinder_50hz.aphi
rem
rem It builds the solver if needed, meshes the .geo if the .msh is missing or
rem older than it, runs the case, and leaves everything in the case's output\.
rem
rem Nothing here writes outside the case folder and the build directory.
rem ---------------------------------------------------------------------------

if "%~1"=="" (
  echo usage: run_case.bat ^<case.aphi^>
  echo.
  echo   run from inside a case folder, e.g.
  echo     cd 01_OneCylinder ^&^& ..\run_case.bat cylinder_50hz.aphi
  exit /b 2
)

set "CASE=%~1"
if not exist "%CASE%" (
  echo error: no such input file: %CASE%
  echo        run this from inside the case folder.
  exit /b 2
)

rem --- locate the repository, which is two levels above this script ----------
set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "REPO=%%~fI"
set "BUILD=%REPO%\build"

rem --- build, unless SKIP_BUILD is set ---------------------------------------
rem build.bat does `cd /d %~dp0` and does not come back, so the case folder has
rem to be restored afterwards or every relative path below resolves against the
rem repository root instead.
rem The build output goes to a log rather than the console -- a full compile is
rem dozens of lines and would bury the solve. But it is SUMMARISED afterwards:
rem silence here used to mean either "recompiled everything" or "did nothing",
rem with no way to tell them apart, which reads as a broken build step.
if defined SKIP_BUILD (
  echo === build      skipped ^(SKIP_BUILD is set^)
  goto :built
)
set "BUILD_LOG=%TEMP%\aphi_build.log"
pushd "%CD%"
call "%REPO%\build.bat" > "%BUILD_LOG%" 2>&1
set "BUILD_RC=%ERRORLEVEL%"
popd
if not "%BUILD_RC%"=="0" (
  echo === build      FAILED. First errors:
  echo.
  findstr /i /c:"error" "%BUILD_LOG%"
  echo.
  echo   full output: %BUILD_LOG%
  exit /b 1
)
findstr /c:"ninja: no work to do" "%BUILD_LOG%" >nul 2>&1
if not errorlevel 1 (
  echo === build      already up to date
) else (
  set "NBUILT=0"
  for /f %%N in ('findstr /c:"Building CXX" "%BUILD_LOG%" 2^>nul ^| find /c /v ""') do set "NBUILT=%%N"
  echo === build      recompiled !NBUILT! file^(s^)
)
:built

rem --- pick the binary: MUMPS build if the case asks for it -----------------
rem `backend = mumps` needs a build configured with -DAPHI_WITH_MUMPS=ON, which
rem the default build is NOT -- that is deliberate, so the ordinary workflow
rem stays free of the oneAPI dependency. A case that asks for MUMPS gets the
rem build-mumps-omp binary instead, with the Intel environment loaded, because
rem that executable links the Fortran runtime (libifcoremd.dll) and will not
rem start without it.
rem
rem 04_Cylinder_SkinDepth is why this exists: 418318 unknowns is not something
rem the internal scalar LDL^T can factor in a sensible time, so the case would
rem otherwise be impossible to put in check.bat at all.
set "SOLVER=%BUILD%\solve_mesh.exe"
set "WANTS_MUMPS="
findstr /r /i /c:"^ *backend *= *mumps" "%CASE%" >nul 2>&1 && set "WANTS_MUMPS=1"
if defined WANTS_MUMPS (
  rem PREFER build-mumps-mkl. It is the same MUMPS on a THREADED MKL, and that
  rem is worth 10x: 266.78 s -> 25.94 s to factor 418318 unknowns, same answer
  rem to every printed digit. The other two link mkl_sequential and leave 20 of
  rem 24 cores idle in the dense frontal kernels where nearly all the flops are.
  rem See docs/MUMPS_SETUP.md.
  set "MBUILD=%REPO%\build-mumps-mkl"
  if not exist "!MBUILD!\solve_mesh.exe" set "MBUILD=%REPO%\build-mumps-omp"
  if not exist "!MBUILD!\solve_mesh.exe" set "MBUILD=%REPO%\build-mumps"
  if not exist "!MBUILD!\solve_mesh.exe" (
    echo error: %CASE% asks for 'backend = mumps' but none of
    echo        %REPO%\build-mumps-mkl, build-mumps-omp or build-mumps
    echo        has a solve_mesh.exe. See docs/MUMPS_SETUP.md to build one, or
    echo        remove the backend line to use the internal solver.
    exit /b 1
  )
  set "SOLVER=!MBUILD!\solve_mesh.exe"
  echo === backend    mumps ^(!MBUILD!^)
  rem setvars.bat is broken on this oneAPI install -- it reports "'vars.bat' is
  rem not recognized" and leaves MKLROOT empty -- so call the component scripts
  rem directly. See docs/MUMPS_SETUP.md.
  set "ONEAPI=%ProgramFiles(x86)%\Intel\oneAPI"
  if exist "!ONEAPI!\compiler\latest\env\vars.bat" (
    call "!ONEAPI!\compiler\latest\env\vars.bat" >nul 2>&1
    call "!ONEAPI!\mkl\latest\env\vars.bat" >nul 2>&1
  ) else (
    echo warning: oneAPI environment not found at !ONEAPI!
    echo          the MUMPS binary needs libifcoremd.dll and may fail to start.
  )
)
if not exist "%SOLVER%" (
  echo error: solver not found at %SOLVER%
  echo        run build.bat in %REPO% first.
  exit /b 1
)

rem --- mesh, if the .msh is missing or older than the .geo -------------------
rem Read the mesh name out of the input file rather than guessing it.
set "MSH="
for /f "tokens=2 delims==" %%A in ('findstr /r /c:"^ *file *=" "%CASE%"') do (
  for /f "tokens=1 delims=#" %%B in ("%%A") do (
    set "MSH=%%B"
  )
)
for /f "tokens=* delims= " %%A in ("!MSH!") do set "MSH=%%A"
rem strip trailing blanks
:trim
if defined MSH if "!MSH:~-1!"==" " set "MSH=!MSH:~0,-1!" & goto :trim

set "GEO=!MSH:.msh=.geo!"
set "NEEDMESH="
if not exist "!MSH!" set "NEEDMESH=1"
if exist "!GEO!" if exist "!MSH!" (
  for /f %%I in ('dir /b /o-d "!GEO!" "!MSH!" 2^>nul') do (
    if not defined NEWEST set "NEWEST=%%I"
  )
  if /i "!NEWEST!"=="!GEO!" set "NEEDMESH=1"
)

if defined NEEDMESH (
  if not exist "!GEO!" (
    echo error: !MSH! is missing and there is no !GEO! to build it from.
    exit /b 1
  )
  echo === meshing !GEO!
  call :find_gmsh
  if not defined GMSH (
    echo error: gmsh not found. Install it, or put gmsh.exe on PATH.
    echo        winget install gmsh.gmsh
    exit /b 1
  )
  "!GMSH!" -3 "!GEO!" -o "!MSH!" -format msh4
  if errorlevel 1 (
    echo error: meshing failed.
    exit /b 1
  )
) else (
  echo === mesh !MSH! is up to date
)

rem --- solve -----------------------------------------------------------------
echo === solving %CASE%
"%SOLVER%" "%CASE%"
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" (
  echo === FAILED with exit code %RC%
  exit /b %RC%
)

echo.
echo === done. outputs:
if exist output dir /b output
exit /b 0

rem ---------------------------------------------------------------------------
:find_gmsh
set "GMSH="
where gmsh >nul 2>&1 && set "GMSH=gmsh" && goto :eof
set "WG=%LOCALAPPDATA%\Microsoft\WinGet\Packages"
for /d %%D in ("%WG%\gmsh.gmsh*") do (
  for /d %%E in ("%%D\gmsh-*") do (
    if exist "%%E\gmsh.exe" set "GMSH=%%E\gmsh.exe"
  )
)
goto :eof
