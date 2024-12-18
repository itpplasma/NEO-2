cmake_minimum_required(VERSION 3.20)

project(NEO-2)
enable_language(C Fortran)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/COMMON/CMAKE_TOOLS")

if(EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include(${CMAKE_SOURCE_DIR}/COMMON/ProjectConfig.cmake.in)

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "Release")
endif()

option(MPI_SUPPORT "Enable MPI support?" ON)
option(MPE_SUPPORT "Enable MPE for graphical MPI work load evaluation?" OFF)

### Compiler flags
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl -heap-arrays -stand f08")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
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

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

find_package(GSL REQUIRED)

include_directories(${NEO2_Inc})
include_directories(${NEO2_Inc}/src)
include_directories(${NEO2_Libs}/include)
include_directories(${CMAKE_BINARY_DIR}/COMMON)

### Projects to build
add_subdirectory(COMMON)
add_subdirectory(NEO-2-PAR)
add_subdirectory(NEO-2-QL)
add_subdirectory(MULTI-SPEC-TOOLS)
add_subdirectory(tools/create_surfaces)

### Optional target for docs
find_package(LATEX)
find_package(Doxygen)
if(LATEX_FOUND AND DOXYGEN_FOUND)
  add_subdirectory(DOC EXCLUDE_FROM_ALL)
endif()
