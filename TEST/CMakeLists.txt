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

# Comprehensive coverage test executable
add_executable(test_spline_coverage
    test_spline_coverage.f90
)

# Set compiler flags
target_compile_options(test_spline_coverage PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_spline_coverage
    common
)

# Include directories
target_include_directories(test_spline_coverage PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME spline_coverage_test 
         COMMAND test_spline_coverage
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(spline_coverage_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All coverage tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some coverage tests FAILED!"
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

# Patch coverage test executable  
add_executable(test_patch_coverage
    test_patch_coverage.f90
)

# Set compiler flags
target_compile_options(test_patch_coverage PRIVATE
    -g -fbacktrace
)

# Link to the common library which contains all our modules
target_link_libraries(test_patch_coverage
    common
)

# Include directories
target_include_directories(test_patch_coverage PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON
)

# Add the test
add_test(NAME patch_coverage_test 
         COMMAND test_patch_coverage
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set test properties
set_tests_properties(patch_coverage_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All patch coverage tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some patch coverage tests FAILED!"
)


