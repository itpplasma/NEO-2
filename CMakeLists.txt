### Initialize CMake (some internal stuff)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(NEO-2-NTV)
enable_language(Fortran)

### Define paths to external libraries (load external file)
include(${CMAKE_SOURCE_DIR}/ProjectConfig.cmake.in)

### MPI
option(MPI_SUPPORT "Should NEO-2 be built with MPI-Support?" ON)
option(MPE_SUPPORT "Should the MPI-library be built with MPE-Support?" ON)
include_directories("${MYMPILIB_PATH}/OBJS")

if (MPI_SUPPORT)
        message(STATUS "MPI-Support activated")
        if (DEFINED MPI_PATH)
                include_directories(${MPI_PATH}/include)
        endif ()
        add_definitions(-DMPI_SUPPORT)
	find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})
        set(MPI_lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -pthread)

endif ()

if (MPE_SUPPORT)
        message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
        if (DEFINED MPE_PATH)
                include_directories(${MPE_PATH}/include)
        endif ()
        add_definitions(-DMPE_SUPPORT)
	find_library(MPE_lib mpe ${MPE_PATH}/lib)
endif ()

### Check Build-Type and set Build-Name extension
if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  # Release mode as default
  message("Building in release mode (default)...")
  set(CMAKE_BUILD_TYPE "RELEASE")
  set(EXE_NAME "neo_2.x")
else(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    message("debug mode")
    set(EXE_NAME "neo_2.x.DEBUG")  
  else(CMAKE_BUILD_TYPE MATCHES RELEASE)
    message("release mode")
    set(CMAKE_BUILD_TYPE "RELEASE")
    set(EXE_NAME "neo_2.x")
  endif(CMAKE_BUILD_TYPE MATCHES DEBUG)
endif(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")

### Define flags for Fortran/C compiler
set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -ffpe-summary=invalid,zero,overflow")
set(MY_DEBUG_FLAG_LIST
  " -O0 -g"   
# " -O0 -g -ggdb -C -p -pg -fbacktrace " 
# " -ffpe-trap=invalid,zero,overflow,underflow "
 " -ffpe-trap=invalid,zero,overflow "
 " -ffpe-summary=invalid,zero,overflow "
 " -fcheck=bounds,mem,array-temps,pointer "
 " -fbacktrace "
# " -Werror " 
# " -fmem-report " 
# " -Waliasing -Wampersand  -Wline-truncation " 
# " -Wnonstd-intrinsics  -Wsurprising -Wno-tabs -Wunderflow " 
 " -Wall "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
# " -Wno-underflow "
# " -Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation " 
  )
string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}") 
set(CFLAGS "-O2 -DDOUBLE_APPEND_FORTRAN")

### Define container for *.o and *.mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

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
find_library(pspline_lib pspline ${PSPLINE_LIB} NO_DEFAULT_PATH)
find_library(gsl_lib gsl ${GSL_LIB} NO_DEFAULT_PATH)
find_library(fgsl_lib fgsl_gfortran ${FGSL_LIB} NO_DEFAULT_PATH)

### HDF5
set(HDF5_DIR /usr/lib/x86_64-linux-gnu/hdf5/serial/)
include_directories(/usr/include/hdf5/serial/)
find_library(HDF5_lib hdf5 ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5_fortran_lib hdf5_fortran ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5HL_lib hdf5_hl ${HDF5_DIR} NO_DEFAULT_PATH)
find_library(HDF5HL_fortran_lib hdf5hl_fortran ${HDF5_DIR} NO_DEFAULT_PATH)

### Source files (load external file)
include (CMakeSources.txt)

### Define executable
add_executable(${EXE_NAME}
	${NEO2_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${FGSL_USERINTER_SRC_FILES}
#	${PSPLINE_USERINTER_SRC_FILES}
)

### Add libraries linked with executable
target_link_libraries(${EXE_NAME} 
        ${fgsl_lib} ${gsl_lib}
	${pspline_lib}
        ${SuperLU_lib}
        ${Umfpack_lib} ${Amd_lib} ${Cholmod_lib} ${Colamd_lib} ${Camd_lib} ${Metis_lib} ${Ccolamd_lib}
        ${MyMPILib_lib}
        ${MPE_lib}
        ${MPI_lib}
        ${HDF5HL_fortran_lib}
        ${HDF5HL_lib}
        ${HDF5_fortran_lib}
        ${HDF5_lib} -lpthread -lz -lrt -lm 
	blas
	lapack
	#debug efence
        )

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
