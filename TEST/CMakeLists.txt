project(NEO-2-TESTS)
enable_testing()

# Coverage flags are now handled globally via CMAKE_BUILD_TYPE=Coverage

# Test executable
# QODO NOTE: spline_cof_original_dense.f90 is intentionally kept as a reference
# implementation for mathematical validation. This ensures continued correctness
# verification and serves as a golden standard for regression testing.
add_executable(test_spline_comparison
    test_spline_comparison.f90
    spline_cof_original_dense.f90
)

# Set compiler flags
target_compile_options(test_spline_comparison PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_comparison
    common
)

# Include directories
target_include_directories(test_spline_comparison PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_comparison_test 
         COMMAND test_spline_comparison
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_comparison_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Unit test executable
add_executable(test_spline_unit
    test_spline_unit.f90
)

# Set compiler flags
target_compile_options(test_spline_unit PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_unit
    common
)

# Include directories
target_include_directories(test_spline_unit PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_unit_test 
         COMMAND test_spline_unit
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_unit_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All large spline tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some large spline tests FAILED!"
)

# Analytical test executable
add_executable(test_spline_analytical
    test_spline_analytical.f90
    spline_cof_original_dense.f90
)

# Set compiler flags
target_compile_options(test_spline_analytical PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_analytical
    common
)

# Include directories
target_include_directories(test_spline_analytical PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_analytical_test 
         COMMAND test_spline_analytical
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_analytical_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All analytical tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some analytical tests FAILED!"
)

# Simple sparse vs dense comparison test executable
add_executable(test_spline_simple
    test_spline_simple.f90
    spline_cof_original_dense.f90
)

# Set compiler flags
target_compile_options(test_spline_simple PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_simple
    common
)

# Include directories
target_include_directories(test_spline_simple PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_simple_test 
         COMMAND test_spline_simple
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_simple_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)


# Error path coverage test executable
add_executable(test_spline_error_paths
    test_spline_error_paths.f90
)

# Set compiler flags
target_compile_options(test_spline_error_paths PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_error_paths
    common
)

# Include directories
target_include_directories(test_spline_error_paths PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_error_paths_test 
         COMMAND test_spline_error_paths
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_error_paths_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All error path tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some error path tests FAILED!"
)


# Sparse module legacy test executable
add_executable(test_sparse_legacy
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_legacy.f90
)

# Set compiler flags
target_compile_options(test_sparse_legacy PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_legacy
    common
)

# Include directories
target_include_directories(test_sparse_legacy PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_legacy_test 
         COMMAND test_sparse_legacy
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_legacy_test PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "Tests passed: *21"
    FAIL_REGULAR_EXPRESSION "STOP 1"
)

# Sparse types test executable
add_executable(test_sparse_types
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_types.f90
)

# Set compiler flags
target_compile_options(test_sparse_types PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_types
    common
)

# Include directories
target_include_directories(test_sparse_types PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_types_test 
         COMMAND test_sparse_types
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_types_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Sparse conversion test executable
add_executable(test_sparse_conversion
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_conversion.f90
)

# Set compiler flags
target_compile_options(test_sparse_conversion PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_conversion
    common
)

# Include directories
target_include_directories(test_sparse_conversion PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_conversion_test 
         COMMAND test_sparse_conversion
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_conversion_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Sparse I/O test executable
add_executable(test_sparse_io
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_io.f90
)

# Set compiler flags
target_compile_options(test_sparse_io PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_io
    common
)

# Include directories
target_include_directories(test_sparse_io PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_io_test 
         COMMAND test_sparse_io
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_io_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Sparse arithmetic test executable
add_executable(test_sparse_arithmetic
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_arithmetic.f90
)

# Set compiler flags
target_compile_options(test_sparse_arithmetic PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_arithmetic
    common
)

# Include directories
target_include_directories(test_sparse_arithmetic PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_arithmetic_test 
         COMMAND test_sparse_arithmetic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_arithmetic_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Sparse solvers test executable
add_executable(test_sparse_solvers
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_solvers.f90
)

# Set compiler flags
target_compile_options(test_sparse_solvers PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_solvers
    common
)

# Include directories
target_include_directories(test_sparse_solvers PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_solvers_test 
         COMMAND test_sparse_solvers
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_solvers_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!"
)

# Sparse utils test executable
add_executable(test_sparse_utils
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sparse_utils.f90
)

# Set compiler flags
target_compile_options(test_sparse_utils PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_sparse_utils
    common
)

# Include directories
target_include_directories(test_sparse_utils PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME sparse_utils_test 
         COMMAND test_sparse_utils
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(sparse_utils_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All sparse utils tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some sparse utils tests FAILED!"
)

# ILU preconditioner test executable
add_executable(test_ilu_precond
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_ilu_precond.f90
)

# Set compiler flags
target_compile_options(test_ilu_precond PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_ilu_precond
    common
)

# Include directories
target_include_directories(test_ilu_precond PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME ilu_precond_test 
         COMMAND test_ilu_precond
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(ilu_precond_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All ILU preconditioner tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some ILU preconditioner tests FAILED!"
)

# BiCGSTAB solver test executable
add_executable(test_bicgstab
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_bicgstab.f90
)

# Set compiler flags
target_compile_options(test_bicgstab PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_bicgstab
    common
)

# Include directories
target_include_directories(test_bicgstab PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME bicgstab_test 
         COMMAND test_bicgstab
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(bicgstab_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All BiCGSTAB tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some BiCGSTAB tests FAILED!"
)

# Solver integration test executable
add_executable(test_solver_integration
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_solver_integration.f90
)

# Set compiler flags
target_compile_options(test_solver_integration PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_solver_integration
    common
)

# Include directories
target_include_directories(test_solver_integration PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME solver_integration_test 
         COMMAND test_solver_integration
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(solver_integration_test PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "All solver integration tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some solver integration tests FAILED!"
)

# BiCGSTAB accuracy test (diagnostic tool, not a regular test)
add_executable(test_spline_bicgstab_accuracy
    test_spline_bicgstab_accuracy.f90
)

# Simple BiCGSTAB test for debugging
add_executable(test_bicgstab_simple
    test_bicgstab_simple.f90
)

# Ill-conditioned BiCGSTAB test
add_executable(test_bicgstab_illcond
    test_bicgstab_illcond.f90
)

# Set compiler flags
target_compile_options(test_spline_bicgstab_accuracy PRIVATE
    -g -fbacktrace
)

target_compile_options(test_bicgstab_simple PRIVATE
    -g -fbacktrace
)

target_compile_options(test_bicgstab_illcond PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_bicgstab_accuracy
    common
)

target_link_libraries(test_bicgstab_simple
    common
)

target_link_libraries(test_bicgstab_illcond
    common
)

# Include directories
target_include_directories(test_spline_bicgstab_accuracy PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

target_include_directories(test_bicgstab_simple PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

target_include_directories(test_bicgstab_illcond PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Adaptive tolerance test executable
add_executable(test_adaptive_tolerance
    test_adaptive_tolerance.f90
)

# Set compiler flags
target_compile_options(test_adaptive_tolerance PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_adaptive_tolerance
    common
)

# Include directories
target_include_directories(test_adaptive_tolerance PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test (expecting failure in RED phase)
add_test(NAME adaptive_tolerance_test 
         COMMAND test_adaptive_tolerance
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties (expecting failure until implementation is complete)
set_tests_properties(adaptive_tolerance_test PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "All adaptive tolerance tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some adaptive tolerance tests FAILED!"
)

# ILU fill level progression test executable  
add_executable(test_ilu_fill_levels
    test_ilu_fill_levels.f90
)

# Set compiler flags
target_compile_options(test_ilu_fill_levels PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_ilu_fill_levels
    common
)

# Include directories
target_include_directories(test_ilu_fill_levels PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME ilu_fill_levels_test 
         COMMAND test_ilu_fill_levels
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(ilu_fill_levels_test PROPERTIES
    TIMEOUT 120
    PASS_REGULAR_EXPRESSION "SUCCESS: BiCGSTAB converged"
    FAIL_REGULAR_EXPRESSION "FAILURE: BiCGSTAB did not converge"
)
