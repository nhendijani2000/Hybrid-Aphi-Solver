# Verification / unit tests.
#
# No external test framework (Catch2/GoogleTest) is fetched here yet: pulling one
# in via CMake FetchContent needs network access at configure time, and this
# project doesn't yet need more than pass/fail checks with a message, so
# test_conditioning.cpp uses a small hand-written runner instead. Revisit this
# once the test suite grows past what that comfortably supports (parameterized
# cases, fixtures, etc. -- Catch2 is a reasonable choice then).

add_executable(aphi_tests test_conditioning.cpp)
target_link_libraries(aphi_tests PRIVATE aphi_core)
add_test(NAME aphi_tests COMMAND aphi_tests)

add_executable(aphi_tests_equilibration test_equilibration.cpp)
target_link_libraries(aphi_tests_equilibration PRIVATE aphi_core)
add_test(NAME aphi_tests_equilibration COMMAND aphi_tests_equilibration)

add_executable(aphi_tests_incidence test_incidence.cpp)
target_link_libraries(aphi_tests_incidence PRIVATE aphi_core)
add_test(NAME aphi_tests_incidence COMMAND aphi_tests_incidence)

add_executable(aphi_tests_gmsh_reader test_gmsh_reader.cpp)
target_link_libraries(aphi_tests_gmsh_reader PRIVATE aphi_core)
# The reader test also reads the two meshes written by Gmsh itself
# (meshes/validation_box_v*.msh), so it needs their location regardless of
# the working directory it is run from. Every other fixture is a string
# literal written to a temp file; these two cannot be, because the whole
# point of them is that this project did not author them.
target_compile_definitions(aphi_tests_gmsh_reader PRIVATE
    APHI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../meshes")
add_test(NAME aphi_tests_gmsh_reader COMMAND aphi_tests_gmsh_reader)

add_executable(aphi_tests_basis_functions test_basis_functions.cpp)
target_link_libraries(aphi_tests_basis_functions PRIVATE aphi_core)
add_test(NAME aphi_tests_basis_functions COMMAND aphi_tests_basis_functions)

add_executable(aphi_tests_tree_cotree test_tree_cotree.cpp)
target_link_libraries(aphi_tests_tree_cotree PRIVATE aphi_core)
# Reads meshes/coax_via.msh to check the tree's exact counts on the real
# target geometry (27,732 free A unknowns).
target_compile_definitions(aphi_tests_tree_cotree PRIVATE
    APHI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../meshes")
add_test(NAME aphi_tests_tree_cotree COMMAND aphi_tests_tree_cotree)

add_executable(aphi_tests_gauge_variants test_gauge_variants.cpp)
target_link_libraries(aphi_tests_gauge_variants PRIVATE aphi_core)
add_test(NAME aphi_tests_gauge_variants COMMAND aphi_tests_gauge_variants)

add_executable(aphi_tests_sparse_matrix test_sparse_matrix.cpp)
target_link_libraries(aphi_tests_sparse_matrix PRIVATE aphi_core)
add_test(NAME aphi_tests_sparse_matrix COMMAND aphi_tests_sparse_matrix)

add_executable(aphi_tests_problem test_problem.cpp)
target_link_libraries(aphi_tests_problem PRIVATE aphi_core)
add_test(NAME aphi_tests_problem COMMAND aphi_tests_problem)

add_executable(aphi_tests_input_file test_input_file.cpp)
target_link_libraries(aphi_tests_input_file PRIVATE aphi_core)
add_test(NAME aphi_tests_input_file COMMAND aphi_tests_input_file)

add_executable(aphi_tests_problem_binding test_problem_binding.cpp)
target_link_libraries(aphi_tests_problem_binding PRIVATE aphi_core)
target_compile_definitions(aphi_tests_problem_binding PRIVATE
    APHI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../meshes")
add_test(NAME aphi_tests_problem_binding COMMAND aphi_tests_problem_binding)

add_executable(aphi_tests_dof_map test_dof_map.cpp)
target_link_libraries(aphi_tests_dof_map PRIVATE aphi_core)
target_compile_definitions(aphi_tests_dof_map PRIVATE
    APHI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../meshes")
add_test(NAME aphi_tests_dof_map COMMAND aphi_tests_dof_map)

add_executable(aphi_tests_quadrature test_quadrature.cpp)
target_link_libraries(aphi_tests_quadrature PRIVATE aphi_core)
add_test(NAME aphi_tests_quadrature COMMAND aphi_tests_quadrature)
