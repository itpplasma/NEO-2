cmake_minimum_required(VERSION 3.20)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (create_surfaces_tool)

set(PROJECT_EXE_NAME "create_surfaces.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -Wall")

set(Interface_Libs
  common
  ${LIBNEO}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  src/create_surfaces.f90
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
  src/nfp.f90
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
