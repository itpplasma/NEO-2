set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
project(HDF5-MERGE-MULTISPEC-OUTPUT)
set(PROJECT_EXE_NAME "h5merge_multispec.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -O2 -Wall")

### Define executable
add_executable(${PROJECT_EXE_NAME}
	./h5merge_multispec.f90
)

set(Interface_Libs
	hdf5::hdf5
	hdf5::hdf5_hl
	hdf5::hdf5_fortran
	hdf5::hdf5_hl_fortran
	BLAS::BLAS
	LAPACK::LAPACK
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
	${Interface_Libs}
)

if(DEFINED ENV{CODE})
	target_link_libraries(${PROJECT_EXE_NAME}
		${LIBNEOLIB}
		${HDF5_tools_lib}
	)
else()
	target_link_libraries(${PROJECT_EXE_NAME}
		LIBNEO::neo
		LIBNEO::hdf5_tools
	)
endif()
