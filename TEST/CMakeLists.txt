project(NEO-2-TESTS)
enable_testing()

# Test executable for boozmn reading (INP_SWI_BOOZMN path)
add_executable(test_boozmn_read test_boozmn_read.f90)
target_link_libraries(test_boozmn_read common ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Generate boozmn_test.nc from circ.bc via bc_to_booz_xform.py.
# circ.bc: 63 surfaces spread 0..1 in s, tokamak, nfp=1.
if(DEFINED ENV{CODE})
  set(CIRC_BC "$ENV{CODE}/NEO-RT/examples/circ.bc")
else()
  set(CIRC_BC "/home/ert/code/NEO-RT/examples/circ.bc")
endif()
# Locate bc_to_booz_xform.py: try the booz-wt eqdsk-booz worktree (task environment)
# then fall back to $CODE/libneo (if it has the eqdsk-boozer branch).
foreach(_candidate
    "/tmp/booz-wt/eqdsk-booz/python/libneo/bc_to_booz_xform.py"
    "$ENV{CODE}/libneo/python/libneo/bc_to_booz_xform.py"
)
  if(EXISTS "${_candidate}")
    set(BC_TO_BOOZMN_SCRIPT "${_candidate}")
    break()
  endif()
endforeach()
set(BOOZMN_TEST_NC "${CMAKE_CURRENT_BINARY_DIR}/boozmn_test.nc")

find_program(PYTHON3_EXE python3)

if(PYTHON3_EXE AND EXISTS "${CIRC_BC}" AND EXISTS "${BC_TO_BOOZMN_SCRIPT}")
  # Generate the reference boozmn at configure time so it is ready before ctest.
  execute_process(
    COMMAND ${PYTHON3_EXE}
            "${BC_TO_BOOZMN_SCRIPT}"
            "${CIRC_BC}"
            "${BOOZMN_TEST_NC}"
    RESULT_VARIABLE GEN_RESULT
    OUTPUT_VARIABLE GEN_OUT
    ERROR_VARIABLE  GEN_ERR
  )
  if(GEN_RESULT EQUAL 0)
    message(STATUS "Generated boozmn_test.nc from circ.bc")
    add_test(NAME boozmn_read_test
             COMMAND test_boozmn_read
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(boozmn_read_test PROPERTIES
        ENVIRONMENT "BOOZMN_TEST_FILE=${BOOZMN_TEST_NC}"
        PASS_REGULAR_EXPRESSION "All tests passed!"
        FAIL_REGULAR_EXPRESSION "FAIL"
        TIMEOUT 30
    )
  else()
    message(WARNING "bc_to_booz_xform.py failed (${GEN_RESULT}): ${GEN_ERR}; boozmn_read_test skipped")
  endif()
else()
  message(WARNING "boozmn_read_test skipped: python3=${PYTHON3_EXE} circ.bc=${CIRC_BC} script=${BC_TO_BOOZMN_SCRIPT}")
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
