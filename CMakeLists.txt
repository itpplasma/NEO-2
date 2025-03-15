cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(NEO-2)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
enable_language(C Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

if(EXISTS ${PROJECT_BINARY_DIR}/BuildConfig.cmake.in)
  include (${PROJECT_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include(${PROJECT_SOURCE_DIR}/COMMON/ProjectConfig.cmake.in)

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "Release")
endif()

option(MPI_SUPPORT "Enable MPI support?" ON)
option(MPE_SUPPORT "Enable MPE for graphical MPI work load evaluation?" OFF)

### Compiler flags
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

add_compile_options(-g -fPIC)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -fbacktrace -fopenmp -ffpe-summary=invalid,zero,overflow")
  if(APPLE)
    add_compile_options(-march=native)
  else()
    add_compile_options(-march=x86-64-v2 -mtune=generic)
  endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -qopenmp -qmkl -heap-arrays -stand f08")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(MY_DEBUG_FLAG_LIST
  " -O0"
  " -fbacktrace "
 " -ffpe-summary=invalid,zero,overflow "
 " -fcheck=all "
 " -Wall "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
 " -Wshadow "
  )
  string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
  set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -mkl -g -warn all -stand f08 ")
endif()

set (CMAKE_Fortran_FLAGS_COVERAGE "--coverage")
set (CMAKE_Fortran_FLAGS_TEST "-O2 --coverage")

### Upstream dependencies
find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

if(MPI_SUPPORT)
  message(STATUS "MPI-Support activated")
  if (DEFINED MPI_PATH)
    include_directories(${MPI_PATH}/include)
  endif ()
  add_definitions(-DMPI_SUPPORT)
endif()

if(MPE_SUPPORT)
  message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
    if (DEFINED MPE_PATH)
      include_directories(${MPE_PATH}/include)
    endif ()
    add_definitions(-DMPE_SUPPORT)
  find_library(MPE_lib mpe ${MPE_PATH}/lib)
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(SuiteSparse)
if(SUITESPARSE_FOUND)
  # In Debian bookworm being standard in 2025, SuiteSparse v5.12 packages don't have
  # their own CMake files for find_package. As soon as this changes,
  # find_package(SuiteSparse) should be replaced by find_package(UMFPACK) and
  # the legacy switch removed.
  set(LEGACY_SUITESPARSE ON)
  message(STATUS "SuiteSparse found.")
else()
  # Use more recent SuiteSparse v7.6.1 from GitHub.
  message(STATUS "SuiteSparse not found. Building.")
  include(suitesparse)
endif()

find_package(HDF5 COMPONENTS C HL Fortran Fortran_HL REQUIRED)

find_package(GSL REQUIRED)

include_directories(${PROJECT_BINARY_DIR}/COMMON)
include_directories(${FGSL_INC})

### Projects to build
add_subdirectory(COMMON)
add_subdirectory(NEO-2-PAR)
add_subdirectory(NEO-2-QL)
add_subdirectory(MULTI-SPEC-TOOLS)
add_subdirectory(tools/create_surfaces)

### Optional target for docs
find_package(LATEX)
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if(LATEX_FOUND AND DOXYGEN_FOUND AND NOT DOXYGEN_DOT_EXECUTABLE STREQUAL "")
  message(STATUS "Documentation will be built.")
  add_subdirectory(DOC EXCLUDE_FROM_ALL)
else()
  message(STATUS "Documentation will not be built.")
endif()
