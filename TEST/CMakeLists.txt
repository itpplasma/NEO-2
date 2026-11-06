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

# theta_B from (R,Z) on a closed curve (used by calc_thetaB_RZloc for
# isw_Vphi_loc=1); analytic-curve oracle.
add_executable(test_theta_rz_solver test_theta_rz_solver.f90)
target_link_libraries(test_theta_rz_solver common)

add_test(NAME theta_rz_solver_test
         COMMAND test_theta_rz_solver
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(theta_rz_solver_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
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

# Reading a stored propagator must allocate its fixed-shape qflux.
add_executable(test_propagator_io test_propagator_io.f90)
target_link_libraries(test_propagator_io neo2_ql)

add_test(NAME propagator_io_test
         COMMAND test_propagator_io
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(propagator_io_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

# isw_calc_Er=2 (prescribed Om_tE) must prepare a multi-species scan from an
# input without V_phi. Drives the real neo_2_ql.x on synthetic HDF5 input.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, h5py"
                  RESULT_VARIABLE _neo2_py_deps OUTPUT_QUIET ERROR_QUIET)
endif()
if(Python3_Interpreter_FOUND AND _neo2_py_deps EQUAL 0)
  add_test(NAME er_direct_prepare_scan_test
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_er_direct_prepare_scan.py
                   $<TARGET_FILE:neo_2_ql.x>
                   ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/neo_ql_prepare_scan.in
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(er_direct_prepare_scan_test PROPERTIES
      PASS_REGULAR_EXPRESSION "All tests passed!"
      FAIL_REGULAR_EXPRESSION "FAIL"
      TIMEOUT 120
  )
else()
  message(STATUS "numpy/h5py not found: skipping er_direct_prepare_scan_test")
endif()

# Restarted complex GMRES of the QL fallback solve.
add_executable(test_ql_gmres test_ql_gmres.f90)
target_link_libraries(test_ql_gmres neo2_ql)
add_test(NAME ql_gmres_test COMMAND test_ql_gmres)
set_tests_properties(ql_gmres_test PROPERTIES
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
    PASS_REGULAR_EXPRESSION "number of eta levels exceeds bsfunc_max_index"
    FAIL_REGULAR_EXPRESSION "FAIL|unexpected"
    TIMEOUT 30
)

# Dry run of the golden-change CI path; skipped without jq.
add_test(NAME golden_change_ci_test
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test_golden_change.sh)
set_tests_properties(golden_change_ci_test PROPERTIES SKIP_RETURN_CODE 77)
