cmake_minimum_required (VERSION 3.18)

project (NEO-2)
enable_language(Fortran C)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/COMMON/CMAKE_TOOLS")

# Upstream dependencies
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

find_package(GSL REQUIRED)
find_package(SuiteSparse REQUIRED)

add_subdirectory(COMMON)
add_subdirectory(NEO-2-PAR)
add_subdirectory(NEO-2-QL)
