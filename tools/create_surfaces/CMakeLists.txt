cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project (create_surfaces_tool)
enable_language(C Fortran)

set(PROJECT_EXE_NAME "create_surfaces.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -Wall")

### Define executables
add_executable(${PROJECT_EXE_NAME} src/create_surfaces.f90)
target_link_libraries(${PROJECT_EXE_NAME} common)

add_executable("nfp.x" src/nfp.f90)
target_link_libraries("nfp.x" common)
