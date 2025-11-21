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

# Sparse module tests

# Test executable for sparse types module
add_executable(test_sparse_types test_sparse_types.f90)
target_compile_options(test_sparse_types PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_types common)
target_include_directories(test_sparse_types PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_types_test 
         COMMAND test_sparse_types
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_types_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!")

# Test executable for sparse conversion module
add_executable(test_sparse_conversion test_sparse_conversion.f90)
target_compile_options(test_sparse_conversion PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_conversion common)
target_include_directories(test_sparse_conversion PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_conversion_test 
         COMMAND test_sparse_conversion
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_conversion_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!")

# Test executable for sparse I/O module
add_executable(test_sparse_io test_sparse_io.f90)
target_compile_options(test_sparse_io PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_io common)
target_include_directories(test_sparse_io PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_io_test 
         COMMAND test_sparse_io
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_io_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!")

# Test executable for sparse arithmetic module
add_executable(test_sparse_arithmetic test_sparse_arithmetic.f90)
target_compile_options(test_sparse_arithmetic PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_arithmetic common)
target_include_directories(test_sparse_arithmetic PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_arithmetic_test 
         COMMAND test_sparse_arithmetic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_arithmetic_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!")

# Test executable for sparse solvers module
add_executable(test_sparse_solvers test_sparse_solvers.f90)
target_compile_options(test_sparse_solvers PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_solvers common)
target_include_directories(test_sparse_solvers PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_solvers_test 
         COMMAND test_sparse_solvers
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_solvers_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some tests FAILED!")

# Test executable for sparse utils module
add_executable(test_sparse_utils test_sparse_utils.f90)
target_compile_options(test_sparse_utils PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_utils common)
target_include_directories(test_sparse_utils PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_utils_test 
         COMMAND test_sparse_utils
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_utils_test PROPERTIES
    TIMEOUT 30
    PASS_REGULAR_EXPRESSION "All sparse utils tests PASSED!"
    FAIL_REGULAR_EXPRESSION "Some sparse utils tests FAILED!")

# Test executable for sparse legacy module
add_executable(test_sparse_legacy test_sparse_legacy.f90)
target_compile_options(test_sparse_legacy PRIVATE -g -fbacktrace)
target_link_libraries(test_sparse_legacy common)
target_include_directories(test_sparse_legacy PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../COMMON
    ${CMAKE_BINARY_DIR}/COMMON)

add_test(NAME sparse_legacy_test 
         COMMAND test_sparse_legacy
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(sparse_legacy_test PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "Tests passed: *21"
    FAIL_REGULAR_EXPRESSION "STOP 1")
