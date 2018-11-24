### Initialize CMake (some internal stuff)
cmake_minimum_required (VERSION 2.8.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../COMMON/CMAKE_TOOLS")

# Enable project
project (NEO-2-Par)
enable_language (Fortran)
enable_testing()

### Build name
set(EXE_NAME "neo_2.x")

### Macro for comparison of version strings
include(../COMMON/CMAKE_TOOLS/MacroCompareVersionStrings.cmake)

# Include BuildConfig, where paths are defined
if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)
#include(${CMAKE_SOURCE_DIR}/ProjectConfig.Debian8.cmake.in)

# MPI is required by this NEO-2 version
find_package(MPI REQUIRED)

# MPE allows graphical work load evaluation of parallel processes 
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" OFF)

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${NEO2_Inc})
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
	set (CMAKE_Fortran_FLAGS_DEBUG "-openmp -g -fcheck=all -Wall -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label -Wshadow")

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

### HDF5
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
# message(STATUS "${HDF5_INCLUDE_DIRS}")
# message(STATUS "${HDF5_LIBRARIES}")
# message(STATUS "${HDF5_Fortran_LIBRARIES}")
# message(STATUS "${HDF5_HL_LIBRARIES}")
# message(STATUS "${HDF5_Fortran_HL_LIBRARIES}")
# message(STATUS "${HDF5_VERSION}")

### SuiteSparse and METIS
find_package(SuiteSparse REQUIRED)
# message(STATUS "${SUITESPARSE_INCLUDE_DIRS}")
# message(STATUS "${SUITESPARSE_LIBRARIES}")
# message(STATUS "${SUITESPARSE_VERSION}")
# message(STATUS "${METIS_LIBRARY}")

### SuperLU
find_package(SuperLU REQUIRED)
# message(STATUS "${SUPERLU_INCLUDES}")
# message(STATUS "${SUPERLU_LIBRARIES}")
# message(STATUS "${SUPERLU_VERSION_VAR}")

### Find libraries
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# GSL
find_library(gsl_lib gsl ${GSL_LIB})
# FGSL
#find_library(fgsl_lib fgsl_gfortran ${FGSL_LIB} NO_DEFAULT_PATH)
find_library(fgsl_lib fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})
# HDF5 interface
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs} NO_DEFAULT_PATH)

### Source files
include (CMakeSources.in)

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

### Define executable
add_executable(${EXE_NAME}
        ${NEO2_SRC_FILES}
	${NEO2_COMMON_FILES}
	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${FGSL_USERINTER_SRC_FILES}
	${MPI_SRC_FILES}
)
set_target_properties(${EXE_NAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

### Add libraries linked with executable
target_link_libraries(${EXE_NAME}
        ${fgsl_lib} ${gsl_lib}
        ${SUPERLU_LIBRARIES}
	${SUITESPARSE_LIBRARIES}
        ${MyMPILib_lib}
        ${MPE_lib}
        ${MPI_Fortran_LIBRARIES}
	${HDF5_tools_lib}
	${HDF5_Fortran_HL_LIBRARIES}
	${HDF5_HL_LIBRARIES}
	${HDF5_Fortran_LIBRARIES}
	${HDF5_LIBRARIES}
	${BLAS_lib} -lz -lpthread
	#debug efence
)

install (TARGETS ${EXE_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
set_property(TARGET ${EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

add_custom_target(
    NEO-2-Versioning
    COMMAND ../ShellScripts/tag_version.sh
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
ADD_TEST(Test_Testcase1-Seq       ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase1" "0")
ADD_TEST(Test_Testcase1-Par       ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase1" "6")
#ADD_TEST(Test_Testcase2-Seq       ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase2" "0")
#ADD_TEST(Test_Testcase2-Par       ${CMAKE_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Testcase2" "6")

#ADD_TEST(Test_Testcase2       ${CMAKE_SOURCE_DIR}/../ShellScripts/test1.sh "Testcase2" "0")
#ADD_TEST(Test_Testcase2-par4  ${CMAKE_SOURCE_DIR}/../ShellScripts/test1.sh "Testcase2" "4")
#ADD_TEST(Test_Testcase3-par6  ${CMAKE_SOURCE_DIR}/../ShellScripts/test1.sh "Testcase3" "6")
