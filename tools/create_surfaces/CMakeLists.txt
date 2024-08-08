### Basic settings
cmake_minimum_required (VERSION 2.8.12)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

### Define the project
project (create_surfaces_tool)
enable_language(Fortran)

### Define the name of the executable
set(PROJECT_EXE_NAME "create_surfaces.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -Wall")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../COMMON/CMAKE_TOOLS")

include(${CMAKE_SOURCE_DIR}/../../COMMON/ProjectConfig.cmake.in)


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# MPI
find_package(MPI REQUIRED)

### SUPERLU ###
find_package(SuperLU REQUIRED)
set(SUPERLU_SRC_FILES
    ${CMAKE_SOURCE_DIR}/../../COMMON/c_fortran_dgssv_ver_5_0.c
    ${CMAKE_SOURCE_DIR}/../../COMMON/c_fortran_zgssv_ver_5_0.c
)
set(SUPERLU_FLAGS "${CFLAGS_DEBUG} ${CFLAGS} -I${SUPERLU_INCLUDES}")
set_source_files_properties(${SUPERLU_SRC_FILES} PROPERTIES COMPILE_FLAGS ${SUPERLU_FLAGS})

### SUITESPARSE ###
find_package(SuiteSparse REQUIRED)
set(SUITESPARSE_SRC_FILES
        ${CMAKE_SOURCE_DIR}/../../COMMON/umf4_f77wrapper_ver_4_5.c
        ${CMAKE_SOURCE_DIR}/../../COMMON/umf4_f77zwrapper_ver_4_5.c
)
set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
set_source_files_properties(${SUITESPARSE_SRC_FILES}  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")

### net-cdf
execute_process(COMMAND nc-config --includedir
                OUTPUT_STRIP_TRAILING_WHITESPACE
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
include_directories(${NETCDFINCLUDE_DIR})
execute_process(COMMAND nc-config --flibs
                OUTPUT_STRIP_TRAILING_WHITESPACE
                OUTPUT_VARIABLE NETCDF_LIBS)

find_library(LIBNEO libneo ${NEO2_Libs} NO_DEFAULT_PATH)
# MyMPIlib
find_library(MyMPILib_lib MyMPILib ${NEO2_Libs}/src/MyMPILib/ NO_DEFAULT_PATH)
#~ include_directories(${CMAKE_SOURCE_DIR}/../../Build-11-1-0/ ${CMAKE_BINARY_DIR})
include_directories(${NEO2_Inc})
#include_directories(${NEO2_Inc}/src/)

set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/../../COMMON/inter_interfaces.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_actual_fluxs.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_actual_spectra.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_control.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_eval_switch.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_exchange.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_input.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_spline_b00.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_spline_data.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_spline.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_sub.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_support.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_units.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_van.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_work.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/nrtype.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/sizey_bo.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/sizey_cur.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/sizey_pla.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/sparse_mod.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/spline_cof.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/spline_int.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/spline_mod.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/test_function.f90

    ${SUPERLU_SRC_FILES}
    ${SUITESPARSE_SRC_FILES}
)

set(Interface_Libs
  ${LIBNEO}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  ../src/create_surfaces.f90
  ${SOURCE_FILES}
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
  ${Interface_Libs}
  ${MyMPILib_lib}
  ${MPI_Fortran_LIBRARIES}
  ${NETCDF_LIBS}
  ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
)

###
add_executable("nfp.x"
  ../src/nfp.f90
  ${SOURCE_FILES}
)
target_link_libraries("nfp.x"
  ${Interface_Libs}
  ${MyMPILib_lib}
  ${MPI_Fortran_LIBRARIES}
  ${NETCDF_LIBS}
  ${SUPERLU_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
)
### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
