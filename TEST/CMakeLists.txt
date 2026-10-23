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

# Existing Coulomb switch contract; native preserves non-default controls.
add_executable(test_coulomb_model_controls test_coulomb_model_controls.f90)
target_link_libraries(test_coulomb_model_controls common)
foreach(model native pitch_angle test_particle full_fp monoenergetic_lorentz)
    add_test(NAME coulomb_controls_${model}
        COMMAND test_coulomb_model_controls ${model})
    set_tests_properties(coulomb_controls_${model} PROPERTIES TIMEOUT 30)
endforeach()
add_test(NAME coulomb_controls_invalid COMMAND test_coulomb_model_controls invalid)
set_tests_properties(coulomb_controls_invalid PROPERTIES WILL_FAIL TRUE TIMEOUT 30)

# Existing UMFPACK refinement selector: preserve default and reject invalid input.
add_executable(test_sparse_method_controls test_sparse_method_controls.f90)
target_link_libraries(test_sparse_method_controls common)
add_test(NAME sparse_method_default
    COMMAND ${CMAKE_COMMAND} -E env --unset=NEO2_SPARSE_SOLVE_METHOD
        $<TARGET_FILE:test_sparse_method_controls> 3)
add_test(NAME sparse_method_refined COMMAND test_sparse_method_controls 2)
set_tests_properties(sparse_method_refined PROPERTIES
    ENVIRONMENT "NEO2_SPARSE_SOLVE_METHOD=2")
add_test(NAME sparse_method_explicit_default COMMAND test_sparse_method_controls 3)
set_tests_properties(sparse_method_explicit_default PROPERTIES
    ENVIRONMENT "NEO2_SPARSE_SOLVE_METHOD=3")
add_test(NAME sparse_method_invalid COMMAND test_sparse_method_controls invalid)
set_tests_properties(sparse_method_invalid PROPERTIES
    ENVIRONMENT "NEO2_SPARSE_SOLVE_METHOD=invalid" WILL_FAIL TRUE)
set_tests_properties(sparse_method_default sparse_method_refined
    sparse_method_explicit_default sparse_method_invalid PROPERTIES TIMEOUT 30)

add_executable(test_operator_ket_actions test_operator_ket_actions.f90)
target_link_libraries(test_operator_ket_actions common)
add_test(NAME operator_ket_actions COMMAND test_operator_ket_actions)
set_tests_properties(operator_ket_actions PROPERTIES TIMEOUT 30)
