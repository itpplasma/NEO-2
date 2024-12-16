cmake_minimum_required(VERSION 3.18)

project(NEO-2)
enable_language(Fortran C)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/COMMON/CMAKE_TOOLS")

if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include(${CMAKE_SOURCE_DIR}/COMMON/ProjectConfig.cmake.in)

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "Release")
endif ()

option(MPI_SUPPORT "Enable MPI support?" ON)
option(MPE_SUPPORT "Enable MPE for graphical MPI work load evaluation?" OFF)

message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

### Upstream dependencies
find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

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

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

find_package(GSL REQUIRED)
find_package(SuiteSparse REQUIRED)

### Own libraries and includes
find_library(LIBNEO neo ${NEO2_Libs} NO_DEFAULT_PATH)
find_library(LIBMAGFIE magfie ${NEO2_Libs} NO_DEFAULT_PATH)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(LIBFGSL fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
find_library(MyMPILib_lib MyMPILib ${NEO2_Libs}/src/MyMPILib/ NO_DEFAULT_PATH)
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs}/src/hdf5_tools/ NO_DEFAULT_PATH)

include_directories(${NEO2_Inc})
include_directories(${NEO2_Inc}/src)
include_directories(${NEO2_Libs}/include)
include_directories(${CMAKE_BINARY_DIR}/COMMON)

### Projects to build
add_subdirectory(COMMON)
add_subdirectory(NEO-2-PAR)
add_subdirectory(NEO-2-QL)
