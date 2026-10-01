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

# isw_calc_Er=2 (prescribed Om_tE) must prepare a multi-species scan from an
# input without V_phi. Drives the real neo_2_ql.x on synthetic HDF5 input.
find_package(Python3 COMPONENTS Interpreter)
set(NEO2_PY_TEST_DEPS FALSE)
if(Python3_Interpreter_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, h5py"
                  RESULT_VARIABLE _neo2_py_deps OUTPUT_QUIET ERROR_QUIET)
  if(_neo2_py_deps EQUAL 0)
    set(NEO2_PY_TEST_DEPS TRUE)
  else()
    message(STATUS "numpy/h5py not found: skipping er_direct tests")
  endif()
endif()
if(NEO2_PY_TEST_DEPS)
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
endif()

# Per-surface isw_calc_Er=2 solve without V_phi. Needs the private golden
# record data; enabled only when NEO2_QL_GOLDEN_DIR points to
# <data>/TESTS/NEO-2/golden_record/ql at configure time.
if(NEO2_PY_TEST_DEPS AND DEFINED ENV{NEO2_QL_GOLDEN_DIR})
  add_test(NAME er_direct_solve_test
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/test_er_direct_solve.py
                   $<TARGET_FILE:neo_2_ql.x> $ENV{NEO2_QL_GOLDEN_DIR}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(er_direct_solve_test PROPERTIES
      PASS_REGULAR_EXPRESSION "All tests passed!"
      FAIL_REGULAR_EXPRESSION "FAIL"
      TIMEOUT 3600
  )
endif()

# Python force-balance checks (issue #75). Registered only when the neo2_ql
# package and its dependencies import at configure time.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    execute_process(COMMAND ${CMAKE_COMMAND} -E env
                            "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
                            ${Python3_EXECUTABLE} -c "import h5py, neo2_ql.force_balance"
                    RESULT_VARIABLE NEO2_PYTHON_DEPS_MISSING
                    OUTPUT_QUIET ERROR_QUIET)
endif()
if(Python3_Interpreter_FOUND AND NOT NEO2_PYTHON_DEPS_MISSING)
    add_test(NAME python_force_balance_test
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/../python/test/test_force_balance.py
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(python_force_balance_test PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
        PASS_REGULAR_EXPRESSION "All tests passed!"
        TIMEOUT 60
    )
    add_test(NAME python_force_balance_aug30835_test
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/../python/test/test_force_balance_aug30835.py
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(python_force_balance_aug30835_test PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
        PASS_REGULAR_EXPRESSION "All tests passed!"
        TIMEOUT 60
    )
else()
    message(STATUS "neo2_ql Python package not importable: skipping python force-balance tests")
endif()
