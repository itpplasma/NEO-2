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

include(${CMAKE_SOURCE_DIR}/../../NEO-2-QL/ProjectConfig.cmake.in)

find_library(libneo_lib libneo ${NEO2_Libs} NO_DEFAULT_PATH)
#~ include_directories(${CMAKE_SOURCE_DIR}/../../Build-11-1-0/ ${CMAKE_BINARY_DIR})
include_directories(${NEO2_Inc})

set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/../../COMMON/inter_interfaces.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_actual_fluxs.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_actual_spectra.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_control.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_eval_switch.f90
    ${CMAKE_SOURCE_DIR}/../../COMMON/neo_exchange.f90
    ../src/neo_input.f90
    ../src/neo_spline_b00.f90
    ../src/neo_spline_data.f90
    ../src/neo_spline.f90
    ../src/neo_sub_mod.f90
    ../src/neo_support.f90
    ../src/neo_units.f90
    ../src/neo_van.f90
    ../src/neo_work.f90
    ../src/nrtype.f90
    ../src/sizey_bo.f90
    ../src/sizey_cur.f90
    ../src/sizey_pla.f90
    ../src/spline_cof.f90
    ../src/spline_int.f90
    ../src/spline_mod.f90
    ../src/test_function.f90
)

set(Interface_Libs
  ${libneo_lib}
  lapack
  blas
)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  ../src/create_surfaces.f90
  ${SOURCE_FILES}
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
  ${Interface_Libs}
)

###
add_executable("nfp.x"
  ../src/nfp.f90
  ${SOURCE_FILES}
)
target_link_libraries("nfp.x"
  ${Interface_Libs}
)
### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
