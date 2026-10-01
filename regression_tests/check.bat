@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"
rem ---------------------------------------------------------------------------
rem Solve every regression case and ASSERT ITS PHYSICS. Run this before making
rem a change to the solver, the assembly or the post-processing, and again
rem after, so the two can be compared.
rem
rem     check.bat                       every case with an expected.txt
rem     check.bat all verify-only       re-check existing output, no solving
rem     check.bat 03_Cylinder_1A_50Hz   one case
rem     check.bat 02_Ansys_Cylinder_50Hz verify-only
rem
rem `verify-only` skips the solve and re-checks what is already on disk, which
rem is what you want while iterating on post-processing -- seconds instead of
rem minutes.
rem
rem WHY THIS IS SEPARATE FROM run-tests.bat. That is unit tests and finishes in
rem seconds, so it runs on every build. This solves real cases: about 6 minutes
rem each with the internal solver, or 20 seconds with `backend = mumps`.
rem Different cadence, different command.
rem
rem WHAT THE CASES COVER, and why running both matters:
rem
rem   01_OneCylinder           10 mm radius at 50 Hz, so a/delta = 1.07 -- REAL
rem                            SKIN EFFECT, and the only case that checks it
rem                            against theory (Kelvin's R_ac/R_dc). It is also
rem                            inductive, omega*L/R = 2.89, where Phi is gauge
rem                            dominated; its expected.txt omits the two Phi
rem                            checks for that reason and says so.
rem   02_Ansys_Cylinder_50Hz   1 V drive. Reads the current out.
rem   03_Cylinder_1A_50Hz      1 A drive, the SAME mesh. Reads the voltage out,
rem                            and checks that exactly 1 A is collected at the
rem                            far terminal.
rem
rem They must report the same R and L -- those belong to the geometry and the
rem material, not to how the thing is driven. A change that breaks the duality
rem shows up as the two cases disagreeing, which neither case alone can see.
rem ---------------------------------------------------------------------------

set "CASE=%~1"
set "MODE=%~2"
if "%CASE%"=="" set "CASE=all"

call :find_pvbatch
if not defined PVBATCH (
  echo error: pvbatch not found. Install ParaView, or put pvbatch on PATH.
  exit /b 2
)

set "FAILED="
set "RAN=0"

if /i not "%CASE%"=="all" (
  if not exist "%CASE%\expected.txt" (
    echo error: %CASE%\expected.txt not found.
    echo        A case without expected values cannot be checked -- see
    echo        02_Ansys_Cylinder_50Hz\expected.txt for the format.
    exit /b 2
  )
  call :one "%CASE%"
  goto :summary
)

for /d %%D in (*) do (
  if exist "%%D\expected.txt" call :one "%%D"
)

:summary
echo.
echo ===========================================================
if "%RAN%"=="0" (
  echo no cases found with an expected.txt
  exit /b 2
)
if defined FAILED (
  echo PHYSICS CHECK FAILED:!FAILED!
  echo.
  echo Do not commit this without understanding which number moved.
  exit /b 1
)
echo all %RAN% case^(s^) OK
exit /b 0

rem ---------------------------------------------------------------------------
:one
set "C=%~1"
set /a RAN+=1
echo.
echo === %C%
if /i "%MODE%"=="verify-only" goto :one_verify

rem The .aphi to run. A case whose expected.txt names one with an `input` line
rem gets that one -- necessary whenever a folder holds several inputs, where
rem guessing picks whichever comes last alphabetically. 01_OneCylinder holds a
rem 1 Hz case, a 50 Hz case and a 4-frequency sweep, and the sweep would write
rem a different frequency into the same output directory that verify.py reads.
rem Otherwise: the case's own, preferring one that is not a backend variant, so
rem `check.bat` exercises the default solver.
set "APHI="
for /f "tokens=2" %%V in ('findstr /r /c:"^input[ 	]" "%C%\expected.txt"') do set "APHI=%%V"
if defined APHI (
  if not exist "%C%\!APHI!" (
    echo     expected.txt names input "!APHI!" but %C%\!APHI! does not exist
    set "FAILED=!FAILED! %C%(bad-input)"
    goto :eof
  )
  echo     input !APHI! ^(from expected.txt^)
)
if not defined APHI (
  for %%F in ("%C%\*.aphi") do (
    echo %%~nF | findstr /i "mumps" >nul || set "APHI=%%~nxF"
  )
)
if not defined APHI (
  echo     no .aphi found in %C%
  set "FAILED=!FAILED! %C%(no-input)"
  goto :eof
)
pushd "%C%"
call "%~dp0run_case.bat" "!APHI!"
set "RC=!ERRORLEVEL!"
popd
if not "!RC!"=="0" (
  echo     SOLVE FAILED ^(exit !RC!^)
  set "FAILED=!FAILED! %C%(solve)"
  goto :eof
)

:one_verify
rem A case may bring its OWN verifier. verify.py is built around a straight
rem rod: it differences Phi between two END CAPS to find the drive, and bins J
rem by distance from the axis. 05_Loop_1A_50Hz is a closed ring with neither --
rem the current circulates and the terminal voltage is a jump ACROSS A CUT --
rem so it ships verify_loop.py instead. Contorting one script to cover both
rem would make neither readable, so a case-local verify_*.py wins where it
rem NOTE: this keeps the LAST match. One verify_*.py per case -- a second one
rem silently takes over as the guard, decided by filename order. Analysis
rem scripts that only print go under another name (see 05's analytic_comparison.py).

rem exists.
set "VERIFIER="
set "VNAME="
for %%F in ("%C%\verify_*.py") do ( set "VERIFIER=%%~fF" & set "VNAME=%%~nxF" )
if defined VERIFIER (
  echo     using the case's own verifier, !VNAME!
  "%PVBATCH%" "!VERIFIER!"
) else (
  "%PVBATCH%" verify.py "%C%"
)
if errorlevel 1 set "FAILED=!FAILED! %C%"
goto :eof

rem ---------------------------------------------------------------------------
:find_pvbatch
set "PVBATCH="
where pvbatch >nul 2>&1 && set "PVBATCH=pvbatch" && goto :eof
for /d %%D in ("%ProgramFiles%\ParaView*") do (
  if exist "%%D\bin\pvbatch.exe" set "PVBATCH=%%D\bin\pvbatch.exe"
)
goto :eof
