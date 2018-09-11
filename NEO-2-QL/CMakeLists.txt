### Initialize CMake (some internal stuff)
cmake_minimum_required(VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../COMMON/CMAKE_TOOLS")

# Enable project
project(NEO-2-NTV)
enable_language(Fortran)

### Check Build-Type and set Build-Name extension

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  # Release mode as default
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif ()

set(EXE_NAME "neo_2.x")

### Macro for comparison of version strings
include(../COMMON/CMAKE_TOOLS/MacroCompareVersionStrings.cmake)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)
#include(${CMAKE_SOURCE_DIR}/ProjectConfig.Debian8.cmake.in)

### MPI is required by this NEO-2 version
find_package(MPI REQUIRED)

# MPE allows graphical work load evaluation of parallel processes 
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" OFF)

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${NEO2_Inc})
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${MYMPILIB_PATH}/OBJS)

### Define flags for Fortran/C compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl")
endif ()

set(MY_DEBUG_FLAG_LIST
  " -O0 -cpp -g"   
# " -O0 -g -ggdb -C -p -pg -fbacktrace " 
# " -ffpe-trap=invalid,zero,overflow,underflow "
 " -ffpe-trap=invalid,zero,overflow "
 " -ffpe-summary=invalid,zero,overflow "
 " -fcheck=all "
# " -fcheck=bounds,mem,array-temps,pointer "
 " -fbacktrace "
# " -Werror " 
# " -fmem-report " 
# " -Waliasing -Wampersand  -Wline-truncation " 
# " -Wnonstd-intrinsics  -Wsurprising -Wno-tabs -Wunderflow " 
 " -Wall "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
 " -Wshadow "
# " -Wno-underflow "
# " -Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation " 
  )
string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}") 
set(CFLAGS "-O2 -DDOUBLE_APPEND_FORTRAN")

### Define container for *.o and *.mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)


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

### HDF5
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

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
include_directories(${FGSL_INC})

# LAPACK and BLAS
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  find_library(lapack_lib lapack ${LAPACK_LIB} NO_DEFAULT_PATH)
  find_library(f77blas_lib f77blas ${F77BLAS_LIB} NO_DEFAULT_PATH)
  find_library(atlas_lib atlas ${ATLASBLAS_LIB} NO_DEFAULT_PATH)
  find_library(open_blas_lapack_lib lapack ${OPEN_BLAS_LAPACK_LIB} NO_DEFAULT_PATH)
  find_library(open_blas_lib blas ${OPEN_BLAS_LIB} NO_DEFAULT_PATH)
endif ()

# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})

# HDF5 interface
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs} NO_DEFAULT_PATH)

### Source files
include(CMakeSources.in)

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
	${COLLOPLIB_FILES}
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(LINALG_LIBRARIES ${lapack_lib} ${f77blas_lib} ${atlas_lib})
endif ()

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
	${LINALG_LIBRARIES}
        )

add_custom_target(
    NEO-2-Versioning
    COMMAND ./Scripts/tag_version.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_dependencies(${EXE_NAME} NEO-2-Versioning)
        
        
        
        
### Get location of the executable (depricated)
#get_target_property(EXE_LOC ${EXE_NAME} LOCATION)
#message(${EXE_LOC})

### Specify command "make copy"
#message(${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(copy
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME}
  ${PROJECT_SOURCE_DIR}/../../TEST_RUNS/${EXE_NAME}
  ) 
      
### Define target directory for "make install" 
#install (TARGETS ${EXE_NAME} DESTINATION /temp/andfmar/)
#set_property(TARGET ${EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
