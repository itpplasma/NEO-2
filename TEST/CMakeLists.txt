project(NEO-2-TESTS)
enable_testing()

add_executable(test_collision_krook test_collision_krook.f90)
add_executable(test_linear_system_diagnostic test_linear_system_diagnostic.f90)
target_link_libraries(test_linear_system_diagnostic common)
add_test(NAME linear_system_diagnostic_disabled COMMAND test_linear_system_diagnostic)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME linear_system_diagnostic_dense_oracle
    COMMAND ${CMAKE_COMMAND} -E env
        "NEO2_DIAGNOSTIC_TEST_EXE=$<TARGET_FILE:test_linear_system_diagnostic>"
        "OMP_NUM_THREADS=1" "OPENBLAS_NUM_THREADS=1"
        ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
        --basetemp=${CMAKE_CURRENT_BINARY_DIR}/diagnostic-pytest
        ${CMAKE_CURRENT_SOURCE_DIR}/test_dense_linear_system_oracle.py)
set_tests_properties(linear_system_diagnostic_dense_oracle PROPERTIES TIMEOUT 30)
target_link_libraries(test_collision_krook common)
add_test(NAME collision_krook_test COMMAND test_collision_krook)
set_tests_properties(collision_krook_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!" TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
foreach(rejection unknown-model negative-rate missing-rate nan-rate integral grid cache
                  zero-trapped negative-passing nan-trapped nan-passing
                  infinite-trapped infinite-passing)
    add_test(NAME collision_krook_rejects_${rejection}
             COMMAND test_collision_krook ${rejection})
    set_tests_properties(collision_krook_rejects_${rejection} PROPERTIES
        WILL_FAIL TRUE TIMEOUT 10
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

add_executable(test_electromagnetic_work test_electromagnetic_work.f90)
target_link_libraries(test_electromagnetic_work neo2_ql)
add_test(NAME electromagnetic_work_test COMMAND test_electromagnetic_work)
set_tests_properties(electromagnetic_work_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL" TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_executable(test_electromagnetic_work_geometry test_electromagnetic_work_geometry.f90)
target_link_libraries(test_electromagnetic_work_geometry neo2_ql)
foreach(EM_LAB 8 10)
    set(EM_GEOMETRY_DIR "${CMAKE_CURRENT_BINARY_DIR}/em_work_lab${EM_LAB}")
    file(MAKE_DIRECTORY "${EM_GEOMETRY_DIR}")
    configure_file(fixtures/em_work_geometry/neo.in.in
                   "${EM_GEOMETRY_DIR}/neo.in" @ONLY)
    add_test(NAME electromagnetic_work_lab${EM_LAB}_geometry_test
             COMMAND test_electromagnetic_work_geometry ${EM_LAB}
             WORKING_DIRECTORY "${EM_GEOMETRY_DIR}")
    set_tests_properties(electromagnetic_work_lab${EM_LAB}_geometry_test PROPERTIES
        PASS_REGULAR_EXPRESSION "All tests passed!"
        FAIL_REGULAR_EXPRESSION "FAIL" TIMEOUT 30
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()

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

add_executable(test_gmres_iterator test_gmres_iterator.f90)
target_link_libraries(test_gmres_iterator neo2_ql)
add_test(NAME gmres_iterator_test COMMAND test_gmres_iterator)
set_tests_properties(gmres_iterator_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!" FAIL_REGULAR_EXPRESSION "FAIL" TIMEOUT 30)
if(MPI_SUPPORT)
    add_test(NAME gmres_iterator_mpi_test
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_gmres_iterator> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(gmres_iterator_mpi_test PROPERTIES
        PASS_REGULAR_EXPRESSION "All tests passed!" FAIL_REGULAR_EXPRESSION "FAIL"
        TIMEOUT 30 PROCESSORS 2)
endif()

add_executable(test_iterator_convergence test_iterator_convergence.f90)
target_link_libraries(test_iterator_convergence neo2_ql)
add_test(NAME iterator_convergence_test COMMAND test_iterator_convergence)
set_tests_properties(iterator_convergence_test PROPERTIES TIMEOUT 30)
if(MPI_SUPPORT)
    add_test(NAME iterator_convergence_mpi_test
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_iterator_convergence> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(iterator_convergence_mpi_test PROPERTIES TIMEOUT 30 PROCESSORS 2)
endif()
