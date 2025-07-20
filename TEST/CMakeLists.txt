project(NEO-2-TESTS)
enable_testing()

# Test executable
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