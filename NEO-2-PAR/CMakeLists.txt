cmake_minimum_required (VERSION 3.18)

project (NEO-2-PAR)
enable_testing()

set(EXE_NAME "neo_2_par.x")

# Generic compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "-O2")

# Compiler dependent settings
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-cpp -ffpe-summary=invalid,zero,overflow")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -Wall -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label -Wshadow")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays")
endif ()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/OBJS)

# Coverage and unit tests
set (CMAKE_Fortran_FLAGS_COVERAGE "--coverage")
set (CMAKE_Fortran_FLAGS_TEST "-O2 --coverage")

### Source files
include (CMakeSources.in)

### Define executable
add_executable(${EXE_NAME}
  ${NEO2_SRC_FILES}
  ${MPI_SRC_FILES}
)
set_target_properties(${EXE_NAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

### Add libraries linked with executable
target_link_libraries(${EXE_NAME} common)
target_include_directories(${EXE_NAME} PUBLIC "$<TARGET_PROPERTY:common,Fortran_MODULE_DIRECTORY>")

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

### Testing section
ADD_TEST(Test_Build_Code      "make")
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
