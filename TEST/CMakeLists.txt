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

add_executable(test_return_distance test_return_distance.f90)
target_link_libraries(test_return_distance common)

add_test(NAME return_distance_test
         COMMAND test_return_distance
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(return_distance_test PROPERTIES
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

add_executable(test_ntv_perturbation_input test_ntv_perturbation_input.f90)
target_link_libraries(test_ntv_perturbation_input neo2_ql)

add_test(NAME ntv_perturbation_input_test
         COMMAND test_ntv_perturbation_input
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ntv_perturbation_input_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_ntv_hel_drive_input test_ntv_hel_drive_input.f90)
target_link_libraries(test_ntv_hel_drive_input neo2_ql)

add_test(NAME ntv_hel_drive_input_test
         COMMAND test_ntv_hel_drive_input
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ntv_hel_drive_input_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_ntv_source_guard test_ntv_source_guard.f90)
target_link_libraries(test_ntv_source_guard neo2_ql)

add_test(NAME ntv_zero_source_guard_test
         COMMAND test_ntv_source_guard zero
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ntv_zero_source_guard_test PROPERTIES
    PASS_REGULAR_EXPRESSION "ERROR: non-axisymmetric solve has no source."
    FAIL_REGULAR_EXPRESSION "FAIL;NaN;zhseqr;returned from rejected source"
    TIMEOUT 5
)

add_test(NAME ntv_driven_source_guard_test
         COMMAND test_ntv_source_guard driven
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ntv_driven_source_guard_test PROPERTIES
    PASS_REGULAR_EXPRESSION "Driven source accepted!"
    FAIL_REGULAR_EXPRESSION "FAIL;NaN;zhseqr;non-axisymmetric solve has no source"
    TIMEOUT 5
)

add_executable(test_phi_divider_parity test_phi_divider_parity.f90)
target_link_libraries(test_phi_divider_parity common)

add_test(NAME phi_divider_parity_test
         COMMAND test_phi_divider_parity
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(phi_divider_parity_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_phi_placement_diagnostic test_phi_placement_diagnostic.f90)
target_link_libraries(test_phi_placement_diagnostic common)

add_test(NAME phi_placement_diagnostic_test
         COMMAND test_phi_placement_diagnostic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(phi_placement_diagnostic_test PROPERTIES
    ENVIRONMENT "NEO2_PHI_PLACEMENT_FILE=${CMAKE_CURRENT_BINARY_DIR}/phi_placement.csv"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_phi_crossing_alignment test_phi_crossing_alignment.f90)
target_link_libraries(test_phi_crossing_alignment neo2_ql)

add_test(NAME phi_crossing_alignment_test
         COMMAND test_phi_crossing_alignment)

set_tests_properties(phi_crossing_alignment_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_crossing_geometry_padding test_crossing_geometry_padding.f90)
target_link_libraries(test_crossing_geometry_padding neo2_ql)

add_test(NAME crossing_geometry_padding_test
         COMMAND test_crossing_geometry_padding
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(crossing_geometry_padding_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_join_failure_diagnostic test_join_failure_diagnostic.f90)
target_link_libraries(test_join_failure_diagnostic neo2_ql)

add_test(NAME join_failure_diagnostic_test
         COMMAND test_join_failure_diagnostic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(join_failure_diagnostic_test PROPERTIES
    ENVIRONMENT "NEO2_JOIN_NORMALIZATION_FILE=${CMAKE_CURRENT_BINARY_DIR}/join_normalization.csv;NEO2_JOIN_END_TRACE_FILE=${CMAKE_CURRENT_BINARY_DIR}/join_end_trace.csv"
    PASS_REGULAR_EXPRESSION "join_ripples: DGBSV info=2"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_lorentz_projection_diagnostic
    test_lorentz_projection_diagnostic.f90)
target_link_libraries(test_lorentz_projection_diagnostic neo2_ql)

add_test(NAME lorentz_projection_diagnostic_test
         COMMAND test_lorentz_projection_diagnostic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(lorentz_projection_diagnostic_test PROPERTIES
    ENVIRONMENT "NEO2_LOCAL_PROJECTION_TRACE_FILE=${CMAKE_CURRENT_BINARY_DIR}/local_projection_trace.csv"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_stage2_attribution test_stage2_attribution.f90)
target_link_libraries(test_stage2_attribution neo2_ql)

add_test(NAME stage2_attribution_test
         COMMAND test_stage2_attribution
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(stage2_attribution_test PROPERTIES
    ENVIRONMENT "NEO2_STAGE2_ATTRIBUTION_FILE=${CMAKE_CURRENT_BINARY_DIR}/stage2_attribution_footprint.csv"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_qflux_interface_diagnostic test_qflux_interface_diagnostic.f90)
target_link_libraries(test_qflux_interface_diagnostic neo2_ql)

add_test(NAME qflux_interface_diagnostic_test
         COMMAND test_qflux_interface_diagnostic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(qflux_interface_diagnostic_test PROPERTIES
    ENVIRONMENT "NEO2_INTERFACE_TRACE_FILE=${CMAKE_CURRENT_BINARY_DIR}/qflux_interface_trace.csv"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_ripple_solver_route test_ripple_solver_route.f90)
target_link_libraries(test_ripple_solver_route neo2_ql)

add_test(NAME ripple_solver_route_test
         COMMAND test_ripple_solver_route
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ripple_solver_route_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_helical_source test_helical_source.f90)
target_link_libraries(test_helical_source neo2_ql)

add_test(NAME helical_source_test
         COMMAND test_helical_source
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(helical_source_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_helical_response test_helical_response.f90)
target_link_libraries(test_helical_response neo2_ql)

add_test(NAME helical_response_test
         COMMAND test_helical_response
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(helical_response_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

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

add_executable(test_qflux_profile test_qflux_profile.f90)
target_link_libraries(test_qflux_profile neo2_ql)

add_test(NAME qflux_profile
         COMMAND test_qflux_profile
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(qflux_profile PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_peak_rss test_peak_rss.f90)
target_link_libraries(test_peak_rss neo2_ql)

add_test(NAME peak_rss_test
         COMMAND test_peak_rss
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(peak_rss_test PROPERTIES
    ENVIRONMENT "NEO2_PEAK_RSS_FILE=${CMAKE_CURRENT_BINARY_DIR}/peak_rss.csv"
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_ntv_output_geometry test_ntv_output_geometry.f90)
target_link_libraries(test_ntv_output_geometry neo2_ql)

add_test(NAME ntv_output_geometry_test
         COMMAND test_ntv_output_geometry
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(ntv_output_geometry_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_collop_zero_endpoint test_collop_zero_endpoint.f90)
target_link_libraries(test_collop_zero_endpoint common)

add_test(NAME collop_zero_endpoint_test
         COMMAND test_collop_zero_endpoint
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(collop_zero_endpoint_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_collop_precompute_mass test_collop_precompute_mass.f90)
target_link_libraries(test_collop_precompute_mass common)

add_test(NAME collop_precompute_mass_test
         COMMAND test_collop_precompute_mass
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(collop_precompute_mass_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 30
)

add_executable(test_collop_bspline_precompute test_collop_bspline_precompute.f90)
target_link_libraries(test_collop_bspline_precompute common)

add_test(NAME collop_bspline_precompute_test
         COMMAND test_collop_bspline_precompute
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(collop_bspline_precompute_test PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL;Floating point exception"
    TIMEOUT 30
)

add_executable(test_cquad_ieee_guard test_cquad_ieee_guard.f90)
target_link_libraries(test_cquad_ieee_guard common)

add_test(NAME cquad_ieee_guard_test
         COMMAND test_cquad_ieee_guard
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(cquad_ieee_guard_test PROPERTIES
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
