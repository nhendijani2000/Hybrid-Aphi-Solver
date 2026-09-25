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
    src/gmsh_reader.cpp
    src/incidence.cpp
    src/basis_functions.cpp
    src/tree_cotree.cpp
    src/gauge_variants.cpp
)
target_include_directories(aphi_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Command-line driver
add_executable(aphi_solver src/main.cpp)
target_link_libraries(aphi_solver PRIVATE aphi_core)

# Developer tool: compares the two tree-cotree gauge variants on a real mesh
# (docs/TREE_COTREE_GAUGE.md, tools/compare_gauges.cpp).
add_executable(compare_gauges tools/compare_gauges.cpp)
target_link_libraries(compare_gauges PRIVATE aphi_core)

if(APHI_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
