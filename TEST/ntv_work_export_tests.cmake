find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_executable(test_ntv_work_export
  ${CMAKE_CURRENT_LIST_DIR}/test_ntv_work_export.f90)
target_link_libraries(test_ntv_work_export PRIVATE ntv_work_export)
set(WORK_EXPORT_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/work_export_fixture)
file(MAKE_DIRECTORY ${WORK_EXPORT_TEST_DIR})
add_test(NAME ntv_work_export_writer COMMAND test_ntv_work_export)
set_tests_properties(ntv_work_export_writer PROPERTIES
  WORKING_DIRECTORY ${WORK_EXPORT_TEST_DIR}
  FIXTURES_SETUP ntv_work_export_fixture TIMEOUT 30)
add_test(NAME ntv_work_export_reader COMMAND ${Python3_EXECUTABLE}
  ${CMAKE_CURRENT_LIST_DIR}/test_ntv_work_reader.py ${WORK_EXPORT_TEST_DIR})
set_tests_properties(ntv_work_export_reader PROPERTIES
  FIXTURES_REQUIRED ntv_work_export_fixture TIMEOUT 30)
