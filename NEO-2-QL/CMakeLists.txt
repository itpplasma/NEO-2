cmake_minimum_required (VERSION 3.18)

project(NEO-2-QL)
enable_testing()

set(EXE_NAME "neo_2_ql.x")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

include(CMakeSources.in)

add_executable(${EXE_NAME}
  ${NEO2_QL_SRC_FILES}
  ${COLLOPLIB_FILES}
)

target_link_libraries(${EXE_NAME} common)
target_include_directories(${EXE_NAME} PUBLIC "$<TARGET_PROPERTY:common,Fortran_MODULE_DIRECTORY>")

add_custom_target(copy
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME}
  ${PROJECT_SOURCE_DIR}/../../TEST_RUNS/${EXE_NAME}
  )

ADD_TEST(Test_Build_Code "make")
ADD_TEST(Extended_QL_Multispecies  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies" "2" "QL")
ADD_TEST(Extended_QL_Multispecies_Read  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies_Read" "2" "QL")
#ADD_TEST(Extended_QL_Multispecies  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_scan.sh "Extended_QL_Scan_Kink" "QL")
