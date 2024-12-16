### Basic settings
cmake_minimum_required (VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../COMMON/CMAKE_TOOLS")

### Define the project
project(HDF5-MERGE-MULTISPEC-OUTPUT)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

### Define the name of the executable
set(PROJECT_EXE_NAME "h5merge_multispec.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -O2 -Wall")

### HDF5
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)

# HDF5 interface
include_directories(${NEO2_Inc})
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs} NO_DEFAULT_PATH)

### Define executable
add_executable(${PROJECT_EXE_NAME}
	./h5merge_multispec.f90
)

set(Interface_Libs
	${HDF5_tools_lib}
	${HDF5_Fortran_HL_LIBRARIES}
	${HDF5_HL_LIBRARIES}
	${HDF5_Fortran_LIBRARIES}
	${HDF5_LIBRARIES}
        lapack
        blas
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
	${Interface_Libs}
)
