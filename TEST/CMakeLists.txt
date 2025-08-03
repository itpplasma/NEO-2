project(NEO-2-TESTS)
enable_testing()

# Test executable for nrutil module
add_executable(test_nrutil test_nrutil.f90)
target_link_libraries(test_nrutil common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Add the tests
add_test(NAME nrutil_test 
         COMMAND test_nrutil
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Add Python test for variable renaming
add_test(NAME variable_renaming_test
         COMMAND python3 test_variable_renaming.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Set test properties
set_tests_properties(nrutil_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

set_tests_properties(variable_renaming_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All variable renaming tests passed!"
    FAIL_REGULAR_EXPRESSION "AssertionError"
    TIMEOUT 30
)
