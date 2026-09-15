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
