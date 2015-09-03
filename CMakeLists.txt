cmake_minimum_required (VERSION 2.8.9)

# Enable project
project (NEO-2-Par)
enable_language (Fortran)
enable_testing()

### Build name
set(EXE_NAME "neo_2.x")

# MPI is required by this NEO-2 version
find_package(MPI REQUIRED)

# MPE allows graphical work load evaluation of parallel processes 
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" OFF)

# Include BuildConfig, where paths are defined
if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include (${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${HDF5_DIR}/include)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${HDF5_TOOLS}/include)
include_directories(${MYMPILIB_PATH}/OBJS)

# Generic compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "-O2")

# Compiler dependent settings
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
	# Determine compiler version
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
        list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
        list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

	if (GCC_MINOR GREATER 8)
		set (CMAKE_Fortran_FLAGS "-cpp -ffpe-summary=invalid,zero,overflow")
	else()
		set (CMAKE_Fortran_FLAGS "-cpp")
	endif()
	set (CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all")     

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays")
endif ()
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

# Create Fortran file with build information
configure_file (
  "${PROJECT_SOURCE_DIR}/cmake_version.f90.in"
  "${PROJECT_SOURCE_DIR}/cmake_version.f90"
  )

# Coverage and unit tests
set (CMAKE_Fortran_FLAGS_COVERAGE "--coverage")
set (CMAKE_Fortran_FLAGS_TEST "-O2 --coverage")

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

### Debug settings
#if(DEBUG)
#	set(DEBUGFLAG "-g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow")
#	set(FFLAG_DEBUG "-Waliasing -Wampersand  -Wline-truncation  -Wsurprising -Wno-tabs  -Wunderflow -cpp")  #-Wnonstd-intrinsics removed in 4.7, had to add -cpp, had to remove -MJOBS
#	set(CDEBUG_DEBUG "-g -ggdb -C -p -fbacktrace")
#	set(CFLAGS_DEBUG "-O0") 
#	set(DBGLIBS "efence")
#endif(DEBUG)

### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(SuperLU_lib superlu_4.1 ${SUPERLU_DIR}/lib NO_DEFAULT_PATH)
find_library(Umfpack_lib umfpack ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(Amd_lib amd ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(Cholmod_lib cholmod ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(Colamd_lib colamd ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(Camd_lib camd ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(Metis_lib metis ${SUITESPARSE_DIR}/lib ${METIS_PATH} NO_DEFAULT_PATH)
find_library(Ccolamd_lib ccolamd ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(SuiteSparseConfig_lib suitesparseconfig ${SUITESPARSE_DIR}/lib NO_DEFAULT_PATH)
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})
find_library(HDF5_tools_lib hdf5_tools ${HDF5_TOOLS}/lib NO_DEFAULT_PATH)
find_library(gsl_lib gsl ${GSL_LIB} NO_DEFAULT_PATH)
find_library(fgsl_lib fgsl_gfortran ${FGSL_LIB} NO_DEFAULT_PATH)
find_library(HDF5_lib hdf5 ${HDF5_DIR}/lib NO_DEFAULT_PATH)
find_library(HDF5_fortran_lib hdf5_fortran ${HDF5_DIR}/lib NO_DEFAULT_PATH)
find_library(HDF5HL_lib hdf5_hl ${HDF5_DIR}/lib NO_DEFAULT_PATH)
find_library(HDF5HL_fortran_lib hdf5hl_fortran ${HDF5_DIR}/lib NO_DEFAULT_PATH)
set(HDF5_Fortran_LIBRARIES ${HDF5_lib} ${HDF5_fortran_lib})
set(HDF5_Fortran_HL_LIBRARIES ${HDF5HL_lib} ${HDF5HL_fortran_lib})

### Source files
include (CMakeSources.in)

### Define executable
add_executable(${EXE_NAME}
	${NEO2_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${FGSL_USERINTER_SRC_FILES}
	${MPI_SRC_FILES}
)
set_target_properties(${EXE_NAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
                                                                                                 
target_link_libraries(${EXE_NAME}
        ${SuperLU_lib}
        ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
	${SuiteSparseConfig_lib}
        ${MyMPILib_lib}
        ${MPE_lib}
	${MPI_Fortran_LIBRARIES}
	${HDF5_tools_lib}
	${HDF5_fortran_lib} ${HDF5HL_fortran_lib}
	${HDF5_lib} ${HDF5HL_lib}
        ${fgsl_lib} ${gsl_lib}
	${BLAS_lib} -lz -lpthread
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

### Testing section
ADD_TEST(Test_Build_Code      "make")
ADD_TEST(Test_Testcase1-Seq       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase1" "0")
ADD_TEST(Test_Testcase1-Par       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase1" "6")
ADD_TEST(Test_Testcase2-Seq       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase2" "0")
ADD_TEST(Test_Testcase2-Par       ${CMAKE_SOURCE_DIR}/Scripts/test_full_reconstruct.sh "Testcase2" "6")

#ADD_TEST(Test_Testcase2       ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase2" "0")
#ADD_TEST(Test_Testcase2-par4  ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase2" "4")
#ADD_TEST(Test_Testcase3-par6  ${CMAKE_SOURCE_DIR}/Scripts/test1.sh "Testcase3" "6")
