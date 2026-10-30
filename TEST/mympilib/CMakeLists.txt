add_executable(test_mympilib
  test_mympilib.f90
  derived_scheduler_module.f90
)
target_link_libraries(test_mympilib PRIVATE NEO2MyMPILib)

add_test(NAME test_mympilib COMMAND test_mympilib)
set_tests_properties(test_mympilib PROPERTIES
  PASS_REGULAR_EXPRESSION "Derived initMaster"
)

add_executable(test_packbuffer_capacity test_packbuffer_capacity.f90)
target_link_libraries(test_packbuffer_capacity PRIVATE NEO2MyMPILib)
add_test(NAME packbuffer_capacity_test
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_packbuffer_capacity> ${MPIEXEC_POSTFLAGS}
)
set_tests_properties(packbuffer_capacity_test PROPERTIES
  PASS_REGULAR_EXPRESSION "All tests passed!"
  PROCESSORS 2
  TIMEOUT 30
)
