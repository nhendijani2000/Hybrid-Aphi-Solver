cmake_minimum_required(VERSION 3.20)
project(aphi_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(APHI_BUILD_TESTS "Build unit/verification tests" ON)

# Core library (headers + implementation)
add_library(aphi_core
    src/complex_matrix.cpp
    src/conditioning.cpp
    src/equilibration.cpp
)
target_include_directories(aphi_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Command-line driver
add_executable(aphi_solver src/main.cpp)
target_link_libraries(aphi_solver PRIVATE aphi_core)

if(APHI_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
