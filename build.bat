@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

REM Builds the solver and its tests.
REM
REM Neither cmake nor the MSVC compiler is on PATH in an ordinary terminal --
REM both ship inside Visual Studio -- so this finds them first. If you are
REM already in a "x64 Native Tools Command Prompt for VS", plain
REM `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release` does the same.

echo === Locating Visual Studio ===
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSPATH="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -property installationPath`) do set "VSPATH=%%i"
)
if not defined VSPATH set "VSPATH=C:\Program Files\Microsoft Visual Studio\18\Community"

if not exist "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" (
  echo.
  echo ERROR: could not find vcvars64.bat under "!VSPATH!"
  echo.
  echo Install the "Desktop development with C++" workload, or open a
  echo "x64 Native Tools Command Prompt for VS" and run cmake there by hand.
  exit /b 1
)
echo     !VSPATH!
call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" >nul

set "CMAKE=cmake"
where cmake >nul 2>&1
if errorlevel 1 set "CMAKE=!VSPATH!\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"

echo === Configuring ===
"!CMAKE!" -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
if errorlevel 1 exit /b 1

echo === Building ===
"!CMAKE!" --build build
if errorlevel 1 exit /b 1

echo.
echo Build complete.
echo.
echo   Run the solver:   build\aphi_solver.exe examples\cylinder_box.aphi
echo   Run the tests:    run-tests.bat
echo.
