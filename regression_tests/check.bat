@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"
rem ---------------------------------------------------------------------------
rem Solve a regression case and ASSERT ITS PHYSICS. Run before committing a
rem change that touches the solver, the assembly or the post-processing.
rem
rem     check.bat                          02_Ansys_Cylinder_50Hz, the default
rem     check.bat 01_OneCylinder
rem     check.bat 02_Ansys_Cylinder_50Hz verify-only
rem
rem `verify-only` skips the solve and re-checks the output already on disk,
rem which is what you want while iterating on post-processing.
rem
rem WHY THIS IS SEPARATE FROM run-tests.bat. run-tests.bat is unit tests and
rem finishes in seconds, so it can run on every build. This solves a real case:
rem about 6 minutes with the internal solver, or 20 seconds with
rem `backend = mumps`. Different cadence, different command.
rem ---------------------------------------------------------------------------

set "CASE=%~1"
if "%CASE%"=="" set "CASE=02_Ansys_Cylinder_50Hz"
set "MODE=%~2"

if not exist "%CASE%\expected.txt" (
  echo error: %CASE%\expected.txt not found.
  echo        A case without expected values cannot be checked -- see
  echo        02_Ansys_Cylinder_50Hz\expected.txt for the format.
  exit /b 2
)

rem --- the .aphi to run: the one the case names, preferring a non-mumps one ---
set "APHI="
for %%F in ("%CASE%\*.aphi") do (
  echo %%~nF | findstr /i "mumps" >nul || set "APHI=%%~nxF"
)
if not defined APHI (
  echo error: no .aphi found in %CASE%
  exit /b 2
)

if /i "%MODE%"=="verify-only" goto :verify

echo === solving %CASE%\%APHI%
pushd "%CASE%"
call "%~dp0run_case.bat" "%APHI%"
set "RC=!ERRORLEVEL!"
popd
if not "!RC!"=="0" (
  echo === SOLVE FAILED ^(exit !RC!^); not verifying
  exit /b 1
)
echo.

:verify
call :find_pvbatch
if not defined PVBATCH (
  echo error: pvbatch not found. Install ParaView, or put pvbatch on PATH.
  exit /b 2
)
"%PVBATCH%" verify.py "%CASE%"
if errorlevel 1 (
  echo.
  echo === PHYSICS CHECK FAILED for %CASE%
  echo     Do not commit this without understanding which number moved.
  exit /b 1
)
echo.
echo === %CASE% OK
exit /b 0

rem ---------------------------------------------------------------------------
:find_pvbatch
set "PVBATCH="
where pvbatch >nul 2>&1 && set "PVBATCH=pvbatch" && goto :eof
for /d %%D in ("%ProgramFiles%\ParaView*") do (
  if exist "%%D\bin\pvbatch.exe" set "PVBATCH=%%D\bin\pvbatch.exe"
)
goto :eof
