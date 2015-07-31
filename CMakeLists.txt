cmake_minimum_required (VERSION 2.6)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (NEO-2-Par)
enable_language (Fortran)
enable_testing()

include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
include (${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)
#include (${CMAKE_SOURCE_DIR}/ProjectVersion.cmake.in)

option(MPI_SUPPORT "Should NEO-2 be built with MPI-Support?" ON)
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" ON)

message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

include_directories("${HDF5_DIR}/include")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

        #gfortran -I/usr/lib/openmpi/include -pthread -I/usr/lib/openmpi/lib -L/usr//lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lhwloc

        include_directories(${INCLUDEDIRS_MPI_GNU})
	include_directories(/usr/include)
        link_directories(${LINKDIRS_MPI_GNU})
	set(MPI_lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl)

        set (CMAKE_Fortran_FLAGS "-cpp -pthread -ffpe-summary=invalid,zero,overflow")
        #set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
        #set (CMAKE_Fortran_FLAGS_DEBUG "")

	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)

	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
        list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
        list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

        #ifort -I/usr/local/openmpi-1.6.3-intel/include -I/usr/local/openmpi-1.6.3-intel/lib -L/usr/local/openmpi-1.6.3-intel/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lm -lnuma -Wl,--export-$

        include_directories(${INCLUDEDIRS_MPI_INTEL})
        link_directories(${LINKDIRS_MPI_INTEL})

        set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays")
        #set (CMAKE_Fortran_FLAGS_RELEASE "")
        #set (CMAKE_Fortran_FLAGS_DEBUG "")

endif ()

if (CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
    list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
    list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

    #message(STATUS ${GCC_MAJOR})
    #message(STATUS ${GCC_MINOR})
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/cmake_version.f90.in"
  "${PROJECT_SOURCE_DIR}/cmake_version.f90"
  )

set (CMAKE_Fortran_FLAGS_COVERAGE "--coverage")
set (CMAKE_Fortran_FLAGS_TEST "-O3 --coverage")

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${MYMPILIB_PATH}/OBJS")

set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

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
endif ()

### Debug settings
#if(DEBUG)
#	set(DEBUGFLAG "-g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow")
#	set(FFLAG_DEBUG "-Waliasing -Wampersand  -Wline-truncation  -Wsurprising -Wno-tabs  -Wunderflow -cpp")  #-Wnonstd-intrinsics removed in 4.7, had to add -cpp, had to remove -MJOBS
#	set(CDEBUG_DEBUG "-g -ggdb -C -p -fbacktrace")
#	set(CFLAGS_DEBUG "-O0") 
#	set(DBGLIBS "efence")
#endif(DEBUG)

### Build name extension
set(EXE_NAME "neo_2.x")

### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(SuperLU_lib superlu_4.1 ${SUPERLU_DIR}lib NO_DEFAULT_PATH)
find_library(Umfpack_lib umfpack ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Amd_lib amd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Cholmod_lib cholmod ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Colamd_lib colamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Camd_lib camd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(Metis_lib metis ${SUITESPARSE_DIR}lib ${METIS_PATH} NO_DEFAULT_PATH)
find_library(Ccolamd_lib ccolamd ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(SuiteSparseConfig_lib suitesparseconfig ${SUITESPARSE_DIR}lib NO_DEFAULT_PATH)
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})
find_library(MPE_lib mpe ${MPE_PATH}/lib)
#find_library(NetCDF_lib netcdf)
#find_library(NetCDFf_lib netcdff)
find_library(HDF5_lib hdf5 ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5_fortran_lib hdf5_fortran ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5HL_lib hdf5_hl ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5HL_fortran_lib hdf5hl_fortran ${HDF5_DIR} NO_DEFAULT_PATH)
#set(NetCDF_dyn_lib -lnetcdff -lnetcdf)# -lhdf5)

### Source files
include (CMakeSources.in)

### Define executable
add_executable(${EXE_NAME}
	${NEO2_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${MPI_SRC_FILES}
)

target_link_libraries(${EXE_NAME}
        ${SuperLU_lib}
        ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
	${SuiteSparseConfig_lib}
        ${MyMPILib_lib}
        ${MPE_lib}
	${MPI_lib}
	${HDF5HL_fortran_lib}
	${HDF5HL_lib}
	${HDF5_fortran_lib} 
	${HDF5_lib} -lpthread -lz -lrt -lm 
	#${NetCDF_dyn_lib}
	${BLAS_lib}
	#debug efence
        )

install (TARGETS ${EXE_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
set_property(TARGET ${EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

add_custom_target(
    NEO-2-Versioning
    COMMAND ./Scripts/tag_version.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_dependencies(${EXE_NAME} NEO-2-Versioning)

if (COPY_TO_ARCHIVE)
	add_custom_command(
		TARGET ${EXE_NAME}
		POST_BUILD
		COMMAND ../Scripts/cp_to_archive.sh
		ARGS ${EXE_NAME} ${EXE_ARCHIVE_PATH} ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_BUILD_TYPE}
        	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		)
endif()

### TESTING
ADD_TEST(Test_Build_Code      "make")
ADD_TEST(Test_Testcase1-Seq       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase1" "0")
ADD_TEST(Test_Testcase1-Par       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase1" "6")
ADD_TEST(Test_Testcase2-Seq       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase2" "0")
ADD_TEST(Test_Testcase2-Par       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase2" "6")


#ADD_TEST(Test_Testcase2       ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase2" "0")
#ADD_TEST(Test_Testcase2-par4  ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase2" "4")
#ADD_TEST(Test_Testcase3-par6  ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase3" "6")
