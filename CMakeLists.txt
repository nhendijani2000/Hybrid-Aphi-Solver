cmake_minimum_required(VERSION 3.20)
project(aphi_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# MSVC does not define M_PI (and other <cmath> math constants) unless this is
# set before <cmath>/<math.h> is first included -- GCC/Clang define them
# unconditionally, which is why this only shows up when building with MSVC.
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

option(APHI_BUILD_TESTS "Build unit/verification tests" ON)

# Core library (headers + implementation)
add_library(aphi_core
    src/complex_matrix.cpp
    src/conditioning.cpp
    src/equilibration.cpp
    src/mesh.cpp
    src/problem.cpp
    src/input_file.cpp
    src/problem_binding.cpp
    src/dof_map.cpp
    src/quadrature.cpp
    src/element_matrix.cpp
    src/sparsity.cpp
    src/assembly.cpp
    src/ordering.cpp
    src/symbolic.cpp
    src/factorization.cpp
    src/postprocess.cpp
    src/gmsh_reader.cpp
    src/incidence.cpp
    src/basis_functions.cpp
    src/tree_cotree.cpp
    src/gauge_variants.cpp
)
target_include_directories(aphi_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)


# ---------------------------------------------------------------------------
# Optional MUMPS backend. OFF by default, and that default is the point: a
# clean checkout builds and runs the whole suite with nothing installed.
# SOLVER_PLAN.md Sec. 12. With it OFF, `backend = mumps` in an input file is a
# clear parse error rather than a silent fall back to the internal solver.
#
#     cmake -S . -B build -DAPHI_WITH_MUMPS=ON -DMUMPS_ROOT=<prefix>
#
# MUMPS is Fortran and needs BLAS/LAPACK, so this wants a Fortran toolchain
# present. See docs/LINEAR_SOLVER.md for what that costs on Windows.
option(APHI_WITH_MUMPS "Link the MUMPS sparse direct solver as an optional backend" OFF)
if(APHI_WITH_MUMPS)
    find_package(MUMPS COMPONENTS zmumps)
    if(NOT MUMPS_FOUND)
        find_library(MUMPS_ZMUMPS_LIB NAMES zmumps HINTS ${MUMPS_ROOT}/lib)
        find_path(MUMPS_INCLUDE_DIR NAMES zmumps_c.h HINTS ${MUMPS_ROOT}/include)
        if(NOT MUMPS_ZMUMPS_LIB OR NOT MUMPS_INCLUDE_DIR)
            message(FATAL_ERROR
                "APHI_WITH_MUMPS=ON but zmumps was not found. Set MUMPS_ROOT to the "
                "install prefix that holds include/zmumps_c.h and lib/zmumps.*, or turn "
                "the option back OFF to build with the internal solver only.")
        endif()
    endif()
    target_compile_definitions(aphi_core PUBLIC APHI_WITH_MUMPS)
    target_sources(aphi_core PRIVATE src/mumps_backend.cpp)
    if(MUMPS_FOUND)
        target_link_libraries(aphi_core PUBLIC MUMPS::zmumps)
    else()
        target_include_directories(aphi_core PUBLIC ${MUMPS_INCLUDE_DIR})
        target_link_libraries(aphi_core PUBLIC ${MUMPS_ZMUMPS_LIB})
    endif()
    message(STATUS "MUMPS backend: ENABLED")
else()
    message(STATUS "MUMPS backend: off (internal solver only)")
endif()

# Command-line driver
add_executable(aphi_solver src/main.cpp)
target_link_libraries(aphi_solver PRIVATE aphi_core)

# Developer tool: compares the two tree-cotree gauge variants on a real mesh
# (docs/TREE_COTREE_GAUGE.md, tools/compare_gauges.cpp).
add_executable(compare_gauges tools/compare_gauges.cpp)
target_link_libraries(compare_gauges PRIVATE aphi_core)

add_executable(solve_mesh tools/solve_mesh.cpp)
target_link_libraries(solve_mesh PRIVATE aphi_core)

if(APHI_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
