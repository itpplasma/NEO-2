### Basic settings
cmake_minimum_required (VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

### Define the project
project(HDF5-MERGE-MULTISPEC-OUTPUT)

### Define the name of the executable
set(PROJECT_EXE_NAME "h5merge_multispec.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -O2 -Wall")

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
	${BLAS_LIBRARIES}
	${LAPACK_LIBRARIES}
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
	${Interface_Libs}
)
