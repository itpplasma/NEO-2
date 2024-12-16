cmake_minimum_required (VERSION 3.18)

project (NEO-2-PAR)
enable_testing()

set(EXE_NAME "neo_2_par.x")


# Include BuildConfig, where paths are defined
if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include(${PROJECT_SOURCE_DIR}/common/ProjectConfig.cmake.in)

# MPI support
option(MPI_SUPPORT "Should the MPI-library be built with MPE-Support?" ON)

# MPE allows graphical work load evaluation of parallel processes
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" OFF)


# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${NEO2_Inc})
include_directories(${NEO2_Inc}/src/)
include_directories(${NEO2_Libs}/include)

# Generic compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "-O2")

# Compiler dependent settings
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # Determine compiler version
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
  string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})

  set (CMAKE_Fortran_FLAGS "-cpp -ffpe-summary=invalid,zero,overflow")
  set (CMAKE_Fortran_FLAGS_DEBUG "-openmp -g -fcheck=all -Wall -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label -Wshadow")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays")
endif ()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/OBJS)

# Create Fortran file with build information
configure_file (
  "${PROJECT_SOURCE_DIR}/common/cmake_version.f90.in"
  "${PROJECT_SOURCE_DIR}/cmake_version.f90"
  )

# Coverage and unit tests
set (CMAKE_Fortran_FLAGS_COVERAGE "--coverage")
set (CMAKE_Fortran_FLAGS_TEST "-O2 --coverage")

if (MPI_SUPPORT)
  message(STATUS "MPI-Support activated")
  if (DEFINED MPI_PATH)
    include_directories(${MPI_PATH}/include)
  endif ()
  add_definitions(-DMPI_SUPPORT)
endif ()

if (MPE_SUPPORT)
  message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
    if (DEFINED MPE_PATH)
      include_directories(${MPE_PATH}/include)
    endif ()
    add_definitions(-DMPE_SUPPORT)
  find_library(MPE_lib mpe ${MPE_PATH}/lib)
endif ()

### Find libraries
find_library(LIBNEO neo ${NEO2_Libs} NO_DEFAULT_PATH)
find_library(LIBMAGFIE magfie ${NEO2_Libs} NO_DEFAULT_PATH)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# FGSL
find_library(LIBFGSL fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${NEO2_Libs}/src/MyMPILib/ NO_DEFAULT_PATH)
# HDF5 interface
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs}/src/hdf5_tools/ NO_DEFAULT_PATH)

### Source files
include (CMakeSources.in)

### Define executable
add_executable(${EXE_NAME}
  ${NEO2_SRC_FILES}
  ${SUITESPARSE_SRC_FILES}
  ${MPI_SRC_FILES}
)
set_target_properties(${EXE_NAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

### Add libraries linked with executable
target_link_libraries(${EXE_NAME}
  common
  ${LIBFGSL} ${GSL_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${LIBNEO}
  ${LIBMAGFIE}
  ${MyMPILib_lib}
  ${MPE_lib}
  ${MPI_Fortran_LIBRARIES}
  ${HDF5_tools_lib}
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${NETCDF_LIBS}
  -lz -lpthread
)
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
