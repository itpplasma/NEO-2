project(NEO-2-TESTS)
enable_testing()

# Test executable for boozmn reading (INP_SWI_BOOZMN path).
add_executable(test_boozmn_read test_boozmn_read.f90)
target_link_libraries(test_boozmn_read common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# The boozmn fixture is committed into the source tree.
# The test reads it via the BOOZMN_TEST_FILE environment variable so no
# configure-time generation, no /tmp references, and no libneo path at runtime.
set(BOOZMN_FIXTURE "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/boozmn_test.nc")

add_test(NAME boozmn_read_test
         COMMAND test_boozmn_read
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(boozmn_read_test PROPERTIES
    ENVIRONMENT "BOOZMN_TEST_FILE=${BOOZMN_FIXTURE}"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

# Test executable for nrutil module
add_executable(test_nrutil test_nrutil.f90)
target_link_libraries(test_nrutil common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(NAME nrutil_test
         COMMAND test_nrutil
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(nrutil_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
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
    TIMEOUT 30
)

add_executable(test_arnoldi_breakdown test_arnoldi_breakdown.f90)
target_link_libraries(test_arnoldi_breakdown neo2_ql)

add_test(NAME arnoldi_breakdown_test
         COMMAND test_arnoldi_breakdown
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(arnoldi_breakdown_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

# eval_bslimit must fail (nonzero exit) when the eta-level index exceeds
# bsfunc_max_index instead of exiting with status 0.
add_executable(test_binarysplit_max_index test_binarysplit_max_index.f90)
target_link_libraries(test_binarysplit_max_index common)

add_test(NAME binarysplit_max_index_test
         COMMAND test_binarysplit_max_index
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(binarysplit_max_index_test PROPERTIES
    WILL_FAIL TRUE
    TIMEOUT 30
)

add_executable(test_fax_basis test_fax_basis.f90)
target_link_libraries(test_fax_basis neo2_ql)
add_test(NAME fax_basis_test COMMAND test_fax_basis)
set_tests_properties(fax_basis_test PROPERTIES TIMEOUT 30)

# Optional native-solve fixtures for the distribution/transport reconstruction oracle.
# Each file is an unmodified lsw_export_fax output, not a fitted reference table.
set(NEO2_FAX_ORACLE_FILES "" CACHE STRING "Native f_AX export files for reconstruction oracle")
if(NEO2_FAX_ORACLE_FILES)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME fax_native_distribution_transport
    COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/python/src/neo2_ql/fax.py ${NEO2_FAX_ORACLE_FILES}
      --out ${CMAKE_CURRENT_BINARY_DIR}/fax_oracle.json)
  set_tests_properties(fax_native_distribution_transport PROPERTIES TIMEOUT 180)
endif()
