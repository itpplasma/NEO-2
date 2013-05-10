### Basic settings
cmake_minimum_required (VERSION 2.6)
project (Neo2-parallel)

### Edit these lines for customization
include (BuildConfig.txt)

set(MPIFLAGS "")
if (${MPI_SUPPORT})
	set(MPIFLAGS "${MPIFLAGS} -DMPI_SUPPORT")
endif (${MPI_SUPPORT})

if (${MPE_SUPPORT})
	set(MPIFLAGS "${MPIFLAGS} -DMPE_SUPPORT")
endif (${MPE_SUPPORT})

if (${CLUSTER_ENV} STREQUAL "ITP")
	set(GFORTRAN_FULL_VERSION 4.7.2)
	set(OPENMPI_FULL_VERSION 1.4.5)
	set(F90 mpif90)
	set(NAME_EXT "gfortran-${GFORTRAN_FULL_VERSION}")

	if (NOT DEFINED MYMPILIB_PATH)
		set(MYMPILIB_PATH /proj/plasma/Libs/MyMPILib)
	endif (NOT DEFINED MYMPILIB_PATH)

	set(OPENMPI_PATH "/proj/plasma/Libs/openmpi/openmpi-${OPENMPI_FULL_VERSION}-gfortran-${GFORTRAN_FULL_VERSION}")
	set(COMPILER_LIBS_PATH "/afs/itp.tugraz.at/opt/gcc/${GFORTRAN_FULL_VERSION}/arch-amd64/lib64")
	set(MPE_PATH "/afs/itp.tugraz.at/opt/mpe/1.3.0")
	
	set(CFLAGS "-O3 -D DOUBLE_APPEND_FORTRAN")
	set(FFLAGS "-O3 -cpp")	
	set(CLUSTER_FLAG "-DITP")

	### Build library paths
	set(LAPACK_LIB -llapack)
	set(BLAS_LIB -lblas)
	set(PROJLIBS /proj/plasma/Libs/)
	set(SUPERLU_VER 4.1)
	set(SUPERLU_DIR ${PROJLIBS}SuperLU/SuperLU_${SUPERLU_VER}/)
	set(SUPERLU_HDR ${SUPERLU_DIR}SRC/)
	set(SUPERLU_F90 ${SUPERLU_DIR}FORTRAN/)
	set(SUPERLU_LIB "-L${SUPERLU_DIR}lib/ -lsuperlu_4.1")
	set(SUITESPARSE_DIR ${PROJLIBS}SuiteSparse_libandinclude/)
	set(SUITESPARSE_HDR ${SUITESPARSE_DIR}include/)
	set(SUITESPARSE_F90 ${SUITESPARSE_DIR}F90/)
	set(SUITESPARSE_LIB "-L${SUITESPARSE_DIR}lib/ -lumfpack -L${SUITESPARSE_DIR}lib/AMD/Lib/ -lamd -L${SUITESPARSE_DIR}lib/ -lcholmod -L${SUITESPARSE_DIR}lib/ -lcolamd -L${SUITESPARSE_DIR}lib/ -lcamd -L${SUITESPARSE_DIR}lib/ -lmetis -L${SUITESPARSE_DIR}lib/ -lccolamd")

        set(MPIFLAGS "${MPIFLAGS} -L${OPENMPI_PATH}/lib ${CLUSTER_FLAG} -Wl,-rpath,${COMPILER_LIBS_PATH}")

endif (${CLUSTER_ENV} STREQUAL "ITP")

if (${CLUSTER_ENV} STREQUAL "VSC2")
	set(F90 mpiifort)
	set(NAME_EXT "ifort")
	set(MYMPILIB_PATH "/home/lv70337/gernot_k/Src/Libs/MyMPILib")
	
	set(CFLAGS "-D DOUBLE_APPEND_FORTRAN")
	set(FFLAGS "-O2 -msse3 -xhost -ftz -cpp")
	set(CLUSTER_FLAG "-DVSC2")
	
	set(PROJLIBS ${PROJECT_SOURCE_DIR}/../../Libs/)
	set(SUPERLU_VER 4.1)
	set(SUPERLU_DIR ${PROJLIBS}SuperLU_${SUPERLU_VER}/)
	set(SUPERLU_HDR ${SUPERLU_DIR}SRC/)
	set(SUPERLU_F90 ${SUPERLU_DIR}FORTRAN/)
	set(SUPERLU_LIB "-L${SUPERLU_DIR}lib/ -lsuperlu_4.1")
	set(SUITESPARSE_DIR ${PROJLIBS}SuiteSparse_3_6_0/)
	set(SUITESPARSE_HDR ${SUITESPARSE_DIR}include/)
	set(SUITESPARSE_F90 ${SUITESPARSE_DIR}F90/)

        set(MPIFLAGS "${MPIFLAGS} ${CLUSTER_FLAG}")

endif (${CLUSTER_ENV} STREQUAL "VSC2")

if (${CLUSTER_ENV} STREQUAL "ZID")
	if (${COMPILER} STREQUAL "IFORT")
	        set(F90 "/usr/local/openmpi-1.6.3-intel/bin/mpif90")
		set(NAME_EXT "ifort")
	endif (${COMPILER} STREQUAL "IFORT")
	if (${COMPILER} STREQUAL "GFORTRAN")
		set(F90 "mpif90")
		set(NAME_EXT "gfortran-4.7.2")
	endif (${COMPILER} STREQUAL "GFORTRAN")

        set(MYMPILIB_PATH "/home/gernot_k/Src/Libs/MyMPILib")
	set(MPE_PATH "/home/gernot_k/Src/Libs/MPE")

        set(CFLAGS "-D DOUBLE_APPEND_FORTRAN")
        set(FFLAGS "-cpp -O2")
        set(CLUSTER_FLAG "-DZID")

        set(PROJLIBS "/home/gernot_k/Src/Libs/")
        set(SUPERLU_VER 4.1)
        set(SUPERLU_DIR ${PROJLIBS}SuperLU_${SUPERLU_VER}/)
        set(SUPERLU_HDR ${SUPERLU_DIR}SRC/)
        set(SUPERLU_F90 ${SUPERLU_DIR}FORTRAN/)
        set(SUPERLU_LIB "-L${SUPERLU_DIR}lib/ -lsuperlu_4.1")
        set(SUITESPARSE_DIR ${PROJLIBS}SuiteSparse_3_6_0/)
        set(SUITESPARSE_HDR ${SUITESPARSE_DIR}include/)
        set(SUITESPARSE_F90 ${SUITESPARSE_DIR}F90/)

        set(MPIFLAGS "${MPIFLAGS} ${CLUSTER_FLAG}")
endif (${CLUSTER_ENV} STREQUAL "ZID")

set(DEBUG 0) # Enable/Disable debug mode

### CMake compiler configuration
set(CMAKE_Fortran_COMPILER ${F90})
enable_language(Fortran)

### MPI Settings
if (${MPI_SUPPORT})
	include_directories(${MYMPILIB_PATH})
	include_directories(${MYMPILIB_PATH}/OBJS)
	#add_subdirectory(${MYMPILIB_PATH})
endif (${MPI_SUPPORT})

### Set compiler flags

set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
if (DEFINED COMPILER_LIBS_PATH)
	set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,${COMPILER_LIBS_PATH}")
endif (DEFINED COMPILER_LIBS_PATH)


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

### Set linker directories
link_directories(${SUITESPARSE_DIR}lib/)
link_directories(${SUITESPARSE_DIR}lib/AMD/Lib/)
link_directories(${SUPERLU_DIR}lib/)

### SUPERLU ###
set(SUPERLU_SRC_FILES
	${SUPERLU_F90}c_fortran_dgssv.c 
	${SUPERLU_F90}c_fortran_zgssv.c
)
set(SUPERLU_FLAGS "${CFLAGS_DEBUG} ${CLFAGS} -I${SUPERLU_HDR}")
set_source_files_properties(${SUPERLU_SRC_FILES} PROPERTIES COMPILE_FLAGS ${SUPERLU_FLAGS})

### SUITESPARSE###
set(SUITESPARSE_SRC_FILES
	${SUITESPARSE_F90}umf4_f77wrapper.c  
	${SUITESPARSE_F90}umf4_f77zwrapper.c
)
set(SUITESPARSE_FLAGS "${CFLAGS} -I${SUITESPARSE_HDR}")
set_source_files_properties(${SUITESPARSE_F90}umf4_f77wrapper.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(${SUITESPARSE_F90}umf4_f77zwrapper.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")


if (${MPI_SUPPORT})
### MPI FILES ###
set (MPI_SRC_FILES
        ./neo2scheduler_module.f90
        ./parallelStorage_module.f90
        ./wuGenericNeo2Workunit_module.f90
)
set_source_files_properties(${MPI_SRC_FILES} PROPERTIES COMPILE_FLAGS "${DEBUGFLAG} ${FFLAG_DEBUG} ${FFLAGS} ${MPIFLAGS}")
endif (${MPI_SUPPORT})

### NEO2 ###
set(NEO2_SRC_FILES
	./field_divB0.f90 
	./bdivfree.f90 
	./spline5_RZ.f90 
	./compute_aiota.f90 
	./neo_modules.f90 
	./nrutil.f90 
	./inter_interfaces.f90 
	./solve_system.f90 
	./test_function.f90 
	./spline_cof.f90 
	./spline_int.f90 
	./pspline.f90 
	./neo_sub.f90 
	./neo_magfie.f90 
	./magfie.f90 
	./neo_mod.f90 
	./polleg.f90 
	./collision_operator.f90 
	./lapack_band.f90 
	./sparse_mod.f90 
	./binarysplit_mod.f90 
	./magnetics.f90 
	./device_mod.f90 
	./mag_interface.f90 
	./fluxsplitter.f90 
	./propagator.f90 
	./prodband.f90 
	./mag.f90 
	./flint.f90 
	./kin_allocate.f90 
	./rhs_kin.f90 
	./magdata_for_particles.f90 
	./rk4_kin.f90 
	./ripple_solver.f90 
	./join_ripples_int.f90 
	./join_ripples.f90 
	./join_ends.f90 
	./plot_distrf.f90 
	./vvn_tok.f90 
	./vvn_w7as.f 
	./vvn_legendre.f 
	./neo2.f90
)
set_source_files_properties(${NEO2_SRC_FILES} PROPERTIES COMPILE_FLAGS "${DEBUGFLAG} ${FFLAG_DEBUG} ${FFLAGS} ${MPIFLAGS}")

### Define executable
add_executable(neo_2.x${NAME_EXT}
	${NEO2_SRC_FILES} 
  	${SUITESPARSE_SRC_FILES}
	${SUPERLU_SRC_FILES}
	${MPI_SRC_FILES}
)

### Linking
if (${CLUSTER_ENV} STREQUAL "ITP")
	target_link_libraries(neo_2.x${NAME_EXT} 
		lapack
		blas
		superlu_4.1 umfpack amd cholmod colamd camd metis ccolamd
		${MYMPILIB_PATH}/libMyMPILib.a
       		${MPE_PATH}/lib/libmpe.a
		${DBGLIBS}
	)
	install (TARGETS neo_2.x${NAME_EXT} DESTINATION /temp/gernot_k/)
endif (${CLUSTER_ENV} STREQUAL "ITP")

if (${CLUSTER_ENV} STREQUAL "VSC2")
	set(MKL_PATH $env(MKLROOT)/lib/intel64)
	set(BLAS_LIB "-Wl,--start-group /opt/goto/libgoto.a $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group  -liomp5")
	set(LAPACK_LIB "/opt/goto/libgoto.a $(MKL_PATH)/libmkl_intel_lp64.a")

	set(CMAKE_EXE_LINKER_FLAGS "-mkl")

	target_link_libraries(neo_2.x${NAME_EXT} ${MYMPILIB_PATH}/libMyMPILib.a mpe /opt/goto/libgoto.a mkl_intel_lp64 mkl_sequential mkl_core 
	${SUPERLU_DIR}lib/libsuperlu_4.1.a
	${SUITESPARSE_DIR}lib/libumfpack.a
	${SUITESPARSE_DIR}lib/libamd.a ${SUITESPARSE_DIR}lib/libcholmod.a ${SUITESPARSE_DIR}lib/libcolamd.a ${SUITESPARSE_DIR}lib/libcamd.a
	${SUITESPARSE_DIR}/../metis-4.0/libmetis.a ${SUITESPARSE_DIR}lib/libccolamd.a ${DBGLIBS}
	)

	install (TARGETS neo_2.x${NAME_EXT} DESTINATION "/home/lv70337/gernot_k/Neo2")
endif (${CLUSTER_ENV} STREQUAL "VSC2")

if (${CLUSTER_ENV} STREQUAL "ZID")
	if (${COMPILER} STREQUAL "GFORTRAN")

        	target_link_libraries(neo_2.x${NAME_EXT}
        		${MYMPILIB_PATH}/libMyMPILib.a
        		${MPE_PATH}/lib/libmpe.a
		        blas
        		lapack
        		${SUPERLU_DIR}lib/libsuperlu_4.1.a
        		${SUITESPARSE_DIR}lib/libumfpack.a
        		${SUITESPARSE_DIR}lib/libamd.a ${SUITESPARSE_DIR}lib/libcholmod.a ${SUITESPARSE_DIR}lib/libcolamd.a ${SUITESPARSE_DIR}lib/libcamd.a
        		${SUITESPARSE_DIR}/../metis-4.0/libmetis.a ${SUITESPARSE_DIR}lib/libccolamd.a ${DBGLIBS}
        	)

		
	endif (${COMPILER} STREQUAL "GFORTRAN")
	
	if (${COMPILER} STREQUAL "IFORT")

        	set(CMAKE_EXE_LINKER_FLAGS "-mkl")
        	target_link_libraries(neo_2.x${NAME_EXT}
        		${MYMPILIB_PATH}/libMyMPILib.a
        		${MPE_PATH}/lib/libmpe.a
        		#mkl_intel_lp64 mkl_sequential mkl_core
		        ${SUPERLU_DIR}lib/libsuperlu_4.1.a
        		${SUITESPARSE_DIR}lib/libumfpack.a
        		${SUITESPARSE_DIR}lib/libamd.a ${SUITESPARSE_DIR}lib/libcholmod.a ${SUITESPARSE_DIR}lib/libcolamd.a ${SUITESPARSE_DIR}lib/libcamd.a
        		${SUITESPARSE_DIR}/../metis-4.0/libmetis.a ${SUITESPARSE_DIR}lib/libccolamd.a ${DBGLIBS}
        	)

	endif (${COMPILER} STREQUAL "IFORT")

        install (TARGETS neo_2.x${NAME_EXT} DESTINATION /home/gernot_k/Neo2)
endif (${CLUSTER_ENV} STREQUAL "ZID")

### Installing
set_property(TARGET neo_2.x${NAME_EXT} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
