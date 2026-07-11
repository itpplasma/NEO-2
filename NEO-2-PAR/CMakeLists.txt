cmake_minimum_required(VERSION 3.24)

project(NEO-2-PAR)
enable_language(C Fortran)
enable_testing()

set(EXE_NAME "neo_2_par.x")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

include (CMakeSources.in)

add_executable(${EXE_NAME}
  ${NEO2_PAR_SRC_FILES}
  ${MPI_SRC_FILES}
)
set_target_properties(${EXE_NAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

target_link_libraries(${EXE_NAME} common)
target_include_directories(${EXE_NAME} PUBLIC "${CMAKE_BINARY_DIR}/COMMON")

install (TARGETS ${EXE_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
set_property(TARGET ${EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

if (COPY_TO_ARCHIVE)
  add_custom_command(
    TARGET ${EXE_NAME}
    POST_BUILD
    COMMAND ../Scripts/cp_to_archive.sh
    ARGS ${EXE_NAME} ${EXE_ARCHIVE_PATH} ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_BUILD_TYPE}
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()

ADD_TEST(NAME Test_Build_Code COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})
set_tests_properties(Test_Build_Code PROPERTIES RESOURCE_LOCK neo2_build_tree)
ADD_TEST(Test_Testcase1-Seq       ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase1" "0")
ADD_TEST(Test_Testcase1-Par       ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase1" "6")
ADD_TEST(Test_Relativistic_BeliaevBudker-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Relativistic_BeliaevBudker" "0")
ADD_TEST(Test_Relativistic_BeliaevBudker-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Relativistic_BeliaevBudker" "6")
ADD_TEST(Test_Relativistic_BraamsKarney-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Relativistic_BraamsKarney" "0")
ADD_TEST(Test_Relativistic_BraamsKarney-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Relativistic_BraamsKarney" "6")

ADD_TEST(Test_par_standard_polynomials-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_standard_polynomials" "0")
ADD_TEST(Test_par_standard_polynomials-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_standard_polynomials" "6")
ADD_TEST(Test_par_quadratic_polynomials-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_quadratic_polynomials" "0")
ADD_TEST(Test_par_quadratic_polynomials-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_quadratic_polynomials" "6")
ADD_TEST(Test_par_cubic_splines-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_cubic_splines" "0")
ADD_TEST(Test_par_cubic_splines-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_cubic_splines" "6")
ADD_TEST(Test_par_general_bsplines-Seq ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_general_bsplines" "0")
ADD_TEST(Test_par_general_bsplines-Par ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "par_general_bsplines" "6")

set(EXTERNAL_DATA_TESTS
  Test_Testcase1-Seq
  Test_Testcase1-Par
  Test_Relativistic_BeliaevBudker-Seq
  Test_Relativistic_BeliaevBudker-Par
  Test_Relativistic_BraamsKarney-Seq
  Test_Relativistic_BraamsKarney-Par
  Test_par_standard_polynomials-Seq
  Test_par_standard_polynomials-Par
  Test_par_quadratic_polynomials-Seq
  Test_par_quadratic_polynomials-Par
  Test_par_cubic_splines-Seq
  Test_par_cubic_splines-Par
  Test_par_general_bsplines-Seq
  Test_par_general_bsplines-Par
)
set_tests_properties(${EXTERNAL_DATA_TESTS} PROPERTIES
  LABELS "regression;external-data"
)
