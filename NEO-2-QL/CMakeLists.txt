### Initialize CMake (some internal stuff)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../COMMON/CMAKE_TOOLS")

# Enable project
project(NEO-2-NTV)
enable_language(Fortran)
enable_testing()

### Check Build-Type and set Build-Name extension

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  # Release mode as default
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif ()

set(EXE_NAME "neo_2.x")

### Macro for comparison of version strings
include(../COMMON/CMAKE_TOOLS/MacroCompareVersionStrings.cmake)

### Define paths to external libraries (load external file)
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  include(${CMAKE_SOURCE_DIR}/common/ProjectConfig.cmake.in)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  include(${CMAKE_SOURCE_DIR}/common/ProjectConfig.cmake.in.intel)
endif ()

### MPI is required by this NEO-2 version
find_package(MPI REQUIRED)

# MPE allows graphical work load evaluation of parallel processes
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" OFF)

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${NEO2_Inc})
include_directories(${NEO2_Inc}/src/)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${NEO2_Libs}/include)

### Define flags for Fortran/C compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl -heap-arrays -stand f08")
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(MY_DEBUG_FLAG_LIST
  " -O0 -cpp -g"
 " -ffpe-summary=invalid,zero,overflow "
 " -fcheck=all "
 " -fbacktrace "
 " -Wall "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
 " -Wshadow "
  )
  string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
  set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -mkl -g -warn all -stand f08 ")
endif ()
set(CFLAGS "-O2 -DDOUBLE_APPEND_FORTRAN")

### Define container for *.o and *.mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)


set (MPI_SUPPORT 1)
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

### HDF5
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

### SuiteSparse and METIS
find_package(SuiteSparse REQUIRED)

### SuperLU
find_package(SuperLU REQUIRED)

### Find libraries
find_library(LIBNEO neo ${NEO2_Libs} NO_DEFAULT_PATH)
find_library(LIBMAGFIE magfie ${NEO2_Libs} NO_DEFAULT_PATH)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# GSL
find_library(LIBGSL gsl)
# FGSL
find_library(LIBFGSL fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
include_directories(${FGSL_INC})

# LAPACK and BLAS
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  find_library(lapack_lib lapack ${LAPACK_LIB} NO_DEFAULT_PATH)
  find_library(f77blas_lib f77blas ${F77BLAS_LIB} NO_DEFAULT_PATH)
  find_library(atlas_lib atlas ${ATLASBLAS_LIB} NO_DEFAULT_PATH)
  find_library(open_blas_lapack_lib lapack ${OPEN_BLAS_LAPACK_LIB} NO_DEFAULT_PATH)
  find_library(open_blas_lib blas ${OPEN_BLAS_LIB} NO_DEFAULT_PATH)
endif ()

# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${NEO2_Libs}/src/MyMPILib/ NO_DEFAULT_PATH)

# HDF5 interface
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs}/src/hdf5_tools/ NO_DEFAULT_PATH)

### Source files
include(CMakeSources.in)

### Select C/Fortran-Interface for SuiteSparse
COMPARE_VERSION_STRINGS(${SUITESPARSE_VERSION} "4.4" ISW_SUITESPARSE_VERSION)
if(${ISW_SUITESPARSE_VERSION} LESS 0)
  message(STATUS "Select C/Fortran-Interface for SuiteSparse version <  4.4 ...")
  set(SUITESPARSE_SRC_FILES "${SUITESPARSE_SRC_FILES_VER_3_6}")
else(${ISW_SUITESPARSE_VERSION} LESS 0)
  message(STATUS "Select C/Fortran-Interface for SuiteSparse version >= 4.4 ...")
  set(SUITESPARSE_SRC_FILES "${SUITESPARSE_SRC_FILES_VER_4_5}")
endif(${ISW_SUITESPARSE_VERSION} LESS 0)

### Select C/Fortran-Interface for SuperLU
COMPARE_VERSION_STRINGS(${SUPERLU_VERSION_VAR} "5.0" ISW_SUPERLU_VERSION)
if(${ISW_SUPERLU_VERSION} LESS 0)
  message(STATUS "Select C/Fortran-Interface for SuperLU version <  5.0 ...")
  set(SUPERLU_SRC_FILES "${SUPERLU_SRC_FILES_VER_4_1}")
else(${ISW_SUPERLU_VERSION} LESS 0)
  message(STATUS "Select C/Fortran-Interface for SuperLU version >= 5.0 ...")
  set(SUPERLU_SRC_FILES "${SUPERLU_SRC_FILES_VER_5_0}")
endif(${ISW_SUPERLU_VERSION} LESS 0)

configure_file (
  "${PROJECT_SOURCE_DIR}/common/cmake_version.f90.in"
  "${PROJECT_SOURCE_DIR}/cmake_version.f90"
  )

### Define executable
add_executable(${EXE_NAME}
        ${NEO2_SRC_FILES}
        ${SUITESPARSE_SRC_FILES}
        ${SUPERLU_SRC_FILES}
        ${COLLOPLIB_FILES}
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(LINALG_LIBRARIES ${open_blas_lib} ${open_blas_lapack_lib} pthread)
endif ()

add_subdirectory(common/)

### Add libraries linked with executable
target_link_libraries(${EXE_NAME}
        common
        ${LIBFGSL} ${LIBGSL}
        ${SUPERLU_LIBRARIES}
        ${SUITESPARSE_LIBRARIES}
        ${MyMPILib_lib}
        ${LIBNEO}
        ${LIBMAGFIE}
        ${MPE_lib}
        ${MPI_Fortran_LIBRARIES}
        ${HDF5_tools_lib}
        ${HDF5_Fortran_HL_LIBRARIES}
        ${HDF5_Fortran_LIBRARIES}
        ${NETCDF_LIBS}
        ${LINALG_LIBRARIES}
        )
target_include_directories(${EXE_NAME} PUBLIC "$<TARGET_PROPERTY:common,Fortran_MODULE_DIRECTORY>")

add_custom_target(
    NEO-2-Versioning
    COMMAND ../ShellScripts/tag_version.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_dependencies(${EXE_NAME} NEO-2-Versioning)


### Specify command "make copy"
add_custom_target(copy
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME}
  ${PROJECT_SOURCE_DIR}/../../TEST_RUNS/${EXE_NAME}
  )

### Testing section
ADD_TEST(Test_Build_Code "make")
ADD_TEST(Extended_QL_Multispecies  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies" "2" "QL")
ADD_TEST(Extended_QL_Multispecies_Read  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies_Read" "2" "QL")
#ADD_TEST(Extended_QL_Multispecies  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_scan.sh "Extended_QL_Scan_Kink" "QL")
