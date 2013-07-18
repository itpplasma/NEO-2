### Basic settings
cmake_minimum_required (VERSION 2.6)

### Edit these lines for customization
include (BuildConfig.txt)
include (Version.txt)

set(DEBUG 0) # Enable/Disable debug mode
set(MPIFLAGS "")
if (${MPI_SUPPORT})
	set(MPIFLAGS "${MPIFLAGS} -DMPI_SUPPORT")
endif (${MPI_SUPPORT})

if (${MPE_SUPPORT})
	set(MPIFLAGS "${MPIFLAGS} -DMPE_SUPPORT")
endif (${MPE_SUPPORT})

if (${CLUSTER_ENV} STREQUAL "ITP")
	include (CMakeITP.txt)
endif (${CLUSTER_ENV} STREQUAL "ITP")

if (${CLUSTER_ENV} STREQUAL "VSC2")
	include (CMakeVSC2.txt)
endif (${CLUSTER_ENV} STREQUAL "VSC2")

if (${CLUSTER_ENV} STREQUAL "ZID")
	include (CMakeZID.txt)
endif (${CLUSTER_ENV} STREQUAL "ZID")

### CMake compiler configuration
set(CMAKE_C_COMPILER ${CC})
set(CMAKE_CXX_COMPILER ${CC})
set(CMAKE_Fortran_COMPILER ${F90})

project (Neo2-parallel)
enable_language(Fortran)
enable_language(C)

### MPI Settings
if (${MPI_SUPPORT})
	include_directories(${MYMPILIB_PATH})
	include_directories(${MYMPILIB_PATH}/OBJS)
	#add_subdirectory(${MYMPILIB_PATH})
endif (${MPI_SUPPORT})

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

### Debug settings
if(DEBUG)
	set(DEBUGFLAG "-g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow")
	set(FFLAG_DEBUG "-Waliasing -Wampersand  -Wline-truncation  -Wsurprising -Wno-tabs  -Wunderflow -cpp")  #-Wnonstd-intrinsics removed in 4.7, had to add -cpp, had to remove -MJOBS
	set(CDEBUG_DEBUG "-g -ggdb -C -p -fbacktrace")
	set(CFLAGS_DEBUG "-O0") 
	set(DBGLIBS "efence")
endif(DEBUG)

### Build name extension
exec_program("uname -m" OUTPUT_VARIABLE ARCH)
set(NAME_EXT ".${NAME_EXT}_${ARCH}")

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

### Source files
include (CMakeSources.txt)

### Define executable
add_executable(neo_2.x${NAME_EXT}
	${NEO2_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${MPI_SRC_FILES}
)

### Linking
if (${CLUSTER_ENV} STREQUAL "ITP")
	include (CMakeITP-linker.txt)
endif (${CLUSTER_ENV} STREQUAL "ITP")

if (${CLUSTER_ENV} STREQUAL "VSC2")
	include (CMakeVSC2-linker.txt)
endif (${CLUSTER_ENV} STREQUAL "VSC2")

if (${CLUSTER_ENV} STREQUAL "ZID")
	include (CMakeZID-linker.txt)
endif (${CLUSTER_ENV} STREQUAL "ZID")
