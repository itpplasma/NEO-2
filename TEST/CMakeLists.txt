project(NEO-2-TESTS)
enable_testing()

add_executable(test_sparse_memory_config test_sparse_memory_config.f90)
target_link_libraries(test_sparse_memory_config PRIVATE common)
add_test(NAME sparse_memory_config COMMAND test_sparse_memory_config)
set_tests_properties(sparse_memory_config PROPERTIES TIMEOUT 1 LABELS "sparse;fast"
    PASS_REGULAR_EXPRESSION "All tests passed!")
add_test(NAME sparse_memory_negative
    COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:test_sparse_memory_config>
        -DARGUMENT=--negative "-DEXPECTED=sparse_host_memory_limit_bytes must be nonnegative"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_sparse_config_failure.cmake)
set_tests_properties(sparse_memory_negative PROPERTIES TIMEOUT 1 LABELS "sparse;fast")
if(NOT NEO2_NATIVE_SPARSE)
    foreach(option cap low)
        add_test(NAME sparse_memory_legacy_${option}
            COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:test_sparse_memory_config>
                -DARGUMENT=--legacy-${option} "-DEXPECTED=memory controls require NEO2_NATIVE_SPARSE=ON"
                -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_sparse_config_failure.cmake)
        set_tests_properties(sparse_memory_legacy_${option} PROPERTIES
            TIMEOUT 1 LABELS "sparse;fast")
    endforeach()
endif()

if(NEO2_NATIVE_SPARSE)
  add_executable(test_native_sparse_memory test_native_sparse_memory.f90)
  target_link_libraries(test_native_sparse_memory PRIVATE common)
  add_test(NAME native_sparse_memory COMMAND test_native_sparse_memory)
  set_tests_properties(native_sparse_memory PROPERTIES
      PASS_REGULAR_EXPRESSION "All tests passed!" TIMEOUT 1 LABELS "sparse;fast")
  foreach(case cap-real cap-complex retained-cap)
    add_test(NAME native_sparse_memory_${case}
        COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:test_native_sparse_memory>
            -DARGUMENT=--${case} -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_native_failure.cmake)
    set_tests_properties(native_sparse_memory_${case} PROPERTIES
        TIMEOUT 1 LABELS "sparse;fast")
  endforeach()
  add_test(NAME native_sparse_memory_mpi
      COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:test_native_sparse_memory>
          -DARGUMENT=--mpi-cap -DMPIEXEC=${MPIEXEC_EXECUTABLE}
          -DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_native_failure.cmake)
  set_tests_properties(native_sparse_memory_mpi PROPERTIES TIMEOUT 8 PROCESSORS 2
      LABELS "sparse;mpi;fast" ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  add_executable(test_native_sparse_adapter test_native_sparse_adapter.f90)
  target_link_libraries(test_native_sparse_adapter PRIVATE common)
  add_test(NAME native_sparse_adapter COMMAND test_native_sparse_adapter)
  set_tests_properties(native_sparse_adapter PROPERTIES
      PASS_REGULAR_EXPRESSION "All tests passed!" TIMEOUT 1 LABELS "sparse;fast")
  add_test(NAME native_sparse_singular
      COMMAND ${CMAKE_COMMAND} -DEXECUTABLE=$<TARGET_FILE:test_native_sparse_adapter>
          -P ${CMAKE_CURRENT_SOURCE_DIR}/expect_native_failure.cmake)
  set_tests_properties(native_sparse_singular PROPERTIES
      TIMEOUT 1 LABELS "sparse;fast")
endif()

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

add_test(NAME arnoldi_projection_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_arnoldi_breakdown> ${MPIEXEC_POSTFLAGS} --distributed)
set_tests_properties(arnoldi_projection_mpi PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!" FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30 PROCESSORS 2 LABELS "arnoldi;mpi;fast"
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

set_tests_properties(arnoldi_breakdown_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

# Generated point kernel checked against independent covariant-field derivatives.
add_executable(test_covariant_curl_leaf test_covariant_curl_leaf.f90)
target_link_libraries(test_covariant_curl_leaf PRIVATE neo_covariant_curl)
add_test(NAME covariant_curl_leaf COMMAND test_covariant_curl_leaf)
set_tests_properties(covariant_curl_leaf PROPERTIES TIMEOUT 1 LABELS "field;fast")

add_executable(test_spline_value_derivatives test_spline_value_derivatives.f90)
target_link_libraries(test_spline_value_derivatives PRIVATE common)
add_test(NAME spline_value_derivatives COMMAND test_spline_value_derivatives)
set_tests_properties(spline_value_derivatives PROPERTIES TIMEOUT 1 LABELS "field;fast")

add_executable(test_periodic_spline_underflow test_periodic_spline_underflow.f90)
target_link_libraries(test_periodic_spline_underflow PRIVATE common)
add_test(NAME periodic_spline_underflow COMMAND test_periodic_spline_underflow)
set_tests_properties(periodic_spline_underflow PROPERTIES TIMEOUT 1 LABELS "field;fast")
