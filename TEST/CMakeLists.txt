project(NEO-2-TESTS)
enable_testing()

# Test executable for boozmn reading (INP_SWI_BOOZMN path)
add_executable(test_boozmn_read test_boozmn_read.f90)
target_link_libraries(test_boozmn_read common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# The test needs a boozmn_test.nc in the working directory.
# If BOOZMN_TEST_FILE is set (absolute path), use it; otherwise skip.
if(DEFINED BOOZMN_TEST_FILE AND EXISTS "${BOOZMN_TEST_FILE}")
  add_test(NAME boozmn_read_test
           COMMAND ${CMAKE_COMMAND}
             -E create_symlink "${BOOZMN_TEST_FILE}"
             "${CMAKE_CURRENT_BINARY_DIR}/boozmn_test.nc"
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME boozmn_read_run
           COMMAND test_boozmn_read
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(boozmn_read_run PROPERTIES
      DEPENDS boozmn_read_test
      PASS_REGULAR_EXPRESSION "All tests passed!"
      FAIL_REGULAR_EXPRESSION "FAIL"
      TIMEOUT 30
  )
endif()

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
