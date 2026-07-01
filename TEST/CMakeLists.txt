project(NEO-2-TESTS)
enable_testing()

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

# Test executable for the fortsparse sparse-solve wrapper
add_executable(test_sparse_solve_fortsparse test_sparse_solve_fortsparse.f90)
target_link_libraries(test_sparse_solve_fortsparse common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(NAME sparse_solve_fortsparse_test
         COMMAND test_sparse_solve_fortsparse
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_solve_fortsparse_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>"
    TIMEOUT 120
)

# Boundary before/after equivalence vs dense LAPACK reference
add_executable(test_sparse_boundary_equivalence test_sparse_boundary_equivalence.f90)
target_link_libraries(test_sparse_boundary_equivalence common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(NAME sparse_boundary_equivalence_test
         COMMAND test_sparse_boundary_equivalence
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_boundary_equivalence_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>"
    TIMEOUT 120
)

add_subdirectory(mympilib)

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

add_executable(test_spline_sparse_assembly test_spline_sparse_assembly.f90)
target_link_libraries(test_spline_sparse_assembly common)

add_test(NAME spline_sparse_assembly_test
         COMMAND test_spline_sparse_assembly
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(spline_sparse_assembly_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    ENVIRONMENT "FORTSPARSE_UMFPACK_HELPER=$<TARGET_FILE:fortsparse_umfpack_helper>"
    TIMEOUT 120
)
