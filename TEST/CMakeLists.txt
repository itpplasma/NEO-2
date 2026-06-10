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

# Value-equivalence oracle test for the former FGSL/GSL wrappers (issue #88).
# GSL is linked here only; no shipped target depends on it.
find_package(GSL QUIET)
if(GSL_FOUND)
  add_executable(test_gsl_replacements test_gsl_replacements.f90)
  target_link_libraries(test_gsl_replacements common GSL::gsl)

  add_test(NAME gsl_replacements_test
           COMMAND test_gsl_replacements
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  set_tests_properties(gsl_replacements_test PROPERTIES
      PASS_REGULAR_EXPRESSION "All tests passed!"
      FAIL_REGULAR_EXPRESSION "FAIL"
      TIMEOUT 60
  )
else()
  message(STATUS "GSL not found: skipping gsl_replacements oracle test")
endif()

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
