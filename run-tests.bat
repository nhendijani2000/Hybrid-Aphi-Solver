@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

REM Runs every test executable and reports a single pass/fail at the end.
REM Each one is a self-contained program that prints "N/M checks passed" and
REM exits non-zero on any failure -- see tests/CMakeLists.txt for why there is
REM no external test framework.

if not exist build\tests\aphi_tests.exe (
  echo Not built yet. Run build.bat first.
  exit /b 1
)

set "TESTS=aphi_tests aphi_tests_equilibration aphi_tests_incidence aphi_tests_gmsh_reader aphi_tests_basis_functions aphi_tests_tree_cotree aphi_tests_gauge_variants aphi_tests_sparse_matrix aphi_tests_problem aphi_tests_input_file aphi_tests_problem_binding aphi_tests_dof_map aphi_tests_quadrature aphi_tests_element_matrix aphi_tests_sparsity aphi_tests_assembly aphi_tests_ordering"

set "FAILED="
for %%T in (%TESTS%) do (
  echo --- %%T
  build\tests\%%T.exe
  if errorlevel 1 set "FAILED=!FAILED! %%T"
)

echo.
if defined FAILED (
  echo FAILED:!FAILED!
  exit /b 1
)
echo All tests passed.
