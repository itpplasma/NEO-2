project(NEO-2-TESTS)
enable_testing()

find_program(PYTHON_EXECUTABLE python3 REQUIRED)

# Test executable for nrutil module
add_executable(test_nrutil test_nrutil.f90)
target_link_libraries(test_nrutil common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Add the test
add_test(NAME nrutil_test
         COMMAND test_nrutil
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(nrutil_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

# Test executable for er_rotation module
add_executable(test_er_rotation test_er_rotation.f90)
target_link_libraries(test_er_rotation neo2_ql)

add_test(NAME er_rotation_test
         COMMAND test_er_rotation
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(er_rotation_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_test(
    NAME python_compute_omte_test
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../python/test/run_compute_omte_checks.py
)
set_tests_properties(python_compute_omte_test PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
    TIMEOUT 60
)
