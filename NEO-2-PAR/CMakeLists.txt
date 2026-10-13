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
target_sources(${EXE_NAME} PRIVATE k_sparse_backend_mod.F90)

option(NEO2_USE_FORTSPARSE "Enable experimental K-only fortsparse backends" OFF)
if(NEO2_USE_FORTSPARSE)
  include(FetchContent)
  FetchContent_Declare(fortsparse
    GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
    GIT_TAG ee14df08b19883f4d91eba4aaddd5c71dbc4eceb)
  FetchContent_MakeAvailable(fortsparse)
  target_compile_definitions(${EXE_NAME} PRIVATE NEO2_HAVE_FORTSPARSE)
  target_link_libraries(${EXE_NAME} fortsparse)
  fortsparse_colocate_helper(${EXE_NAME})
endif()

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
