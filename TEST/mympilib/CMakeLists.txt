add_executable(test_mympilib
  test_mympilib.f90
  derived_scheduler_module.f90
)
target_link_libraries(test_mympilib PRIVATE NEO2MyMPILib)

add_test(NAME test_mympilib COMMAND test_mympilib)
set_tests_properties(test_mympilib PROPERTIES
  PASS_REGULAR_EXPRESSION "Derived initMaster"
)
