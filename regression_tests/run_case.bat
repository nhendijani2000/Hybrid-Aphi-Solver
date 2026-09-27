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
if defined SKIP_BUILD goto :built
echo === building
pushd "%CD%"
call "%REPO%\build.bat" >nul 2>&1
set "BUILD_RC=%ERRORLEVEL%"
popd
if not "%BUILD_RC%"=="0" (
  echo error: build failed. Run build.bat directly to see why.
  exit /b 1
)
:built

set "SOLVER=%BUILD%\solve_mesh.exe"
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
