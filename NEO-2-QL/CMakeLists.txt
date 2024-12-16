cmake_minimum_required (VERSION 3.18)

# Enable project
project(NEO-2-NTV)
enable_testing()

### Check Build-Type and set Build-Name extension

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  # Release mode as default
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif ()

set(EXE_NAME "neo_2_ql.x")

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

# Include directories for NEO-2
include_directories(${NEO2_Inc})
include_directories(${NEO2_Inc}/src/)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${NEO2_Libs}/include)

### Define flags for Fortran/C compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl -heap-arrays -stand f08")
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
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
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -mkl -g -warn all -stand f08 ")
endif ()
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

### Find libraries
find_library(LIBNEO neo ${NEO2_Libs} NO_DEFAULT_PATH)
find_library(LIBMAGFIE magfie ${NEO2_Libs} NO_DEFAULT_PATH)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# GSL
find_library(LIBGSL gsl)
# FGSL
find_library(LIBFGSL fgsl ${FGSL_LIB} NO_DEFAULT_PATH)
include_directories(${FGSL_INC})

# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${NEO2_Libs}/src/MyMPILib/ NO_DEFAULT_PATH)

# HDF5 interface
find_library(HDF5_tools_lib hdf5_tools ${NEO2_Libs}/src/hdf5_tools/ NO_DEFAULT_PATH)

### Source files
include(CMakeSources.in)

configure_file (
  "${PROJECT_SOURCE_DIR}/common/cmake_version.f90.in"
  "${PROJECT_SOURCE_DIR}/cmake_version.f90"
  )

### Define executable
add_executable(${EXE_NAME}
        ${NEO2_SRC_FILES}
        ${SUITESPARSE_SRC_FILES}
        ${COLLOPLIB_FILES}
)

### Add libraries linked with executable
target_link_libraries(${EXE_NAME}
  common
  ${LIBFGSL} ${GSL_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${LIBNEO}
  ${LIBMAGFIE}
  ${MyMPILib_lib}
  ${MPE_lib}
  ${MPI_Fortran_LIBRARIES}
  ${HDF5_tools_lib}
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${NETCDF_LIBS}
  -lz -lpthread
)
target_include_directories(${EXE_NAME} PUBLIC "$<TARGET_PROPERTY:common,Fortran_MODULE_DIRECTORY>")

### Specify command "make copy"
add_custom_target(copy
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME}
  ${PROJECT_SOURCE_DIR}/../../TEST_RUNS/${EXE_NAME}
  )

### Testing section
ADD_TEST(Test_Build_Code "make")
ADD_TEST(Extended_QL_Multispecies  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies" "2" "QL")
ADD_TEST(Extended_QL_Multispecies_Read  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies_Read" "2" "QL")
#ADD_TEST(Extended_QL_Multispecies  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_scan.sh "Extended_QL_Scan_Kink" "QL")
