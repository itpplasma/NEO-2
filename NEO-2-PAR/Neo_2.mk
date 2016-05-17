#F90 = f95-lah
VER = 
F90 = gfortran
#VER = -4.7.0
DEBUG =   

OBJS =	\
	OBJS/field_divB0.o \
	OBJS/bdivfree.o \
	OBJS/spline5_RZ.o \
	OBJS/compute_aiota.o \
	OBJS/neo_modules.o \
	OBJS/nrutil.o \
	OBJS/inter_interfaces.o \
	OBJS/solve_system.o \
	OBJS/test_function.o \
	OBJS/spline_cof.o \
	OBJS/spline_int.o \
	OBJS/pspline.o \
	OBJS/neo_sub.o \
	OBJS/neo_magfie.o \
	OBJS/magfie.o \
	OBJS/neo_mod.o \
	OBJS/polleg.o \
	OBJS/collision_operator.o \
	OBJS/lapack_band.o \
	OBJS/sparse_mod.o \
	OBJS/c_fortran_dgssv.o \
	OBJS/c_fortran_zgssv.o \
	OBJS/umf4_f77wrapper64.o \
	OBJS/umf4_f77zwrapper64.o \
	OBJS/binarysplit_mod.o \
	OBJS/magnetics.o \
	OBJS/device_mod.o \
	OBJS/mag_interface.o \
	OBJS/fluxsplitter.o \
	OBJS/propagator.o \
	OBJS/prodband.o \
	OBJS/mag.o \
	OBJS/flint.o \
	OBJS/kin_allocate.o \
	OBJS/rhs_kin.o \
	OBJS/magdata_for_particles.o \
	OBJS/rk4_kin.o \
	OBJS/ripple_solver.o \
	OBJS/join_ripples_int.o \
	OBJS/join_ripples.o \
	OBJS/join_ends.o \
	OBJS/plot_distrf.o \
	OBJS/vvn_tok.o \
	OBJS/vvn_w7as.o \
	OBJS/vvn_legendre.o \
	OBJS/neo2.o


LAPACK_LIB  = -llapack 
BLAS_LIB    = -lblas

# Project Libraries
PROJLIBS = /proj/plasma/Libs/
# SuperLU V. 4.1 and SuiteSparse V. 3.6.0
# Fortran-Interface of SuiteSparse is only tested 
SUPERLU_VER = 4.1
SUPERLU_DIR = $(PROJLIBS)SuperLU/SuperLU_$(SUPERLU_VER)/
SUPERLU_HDR = $(SUPERLU_DIR)SRC/
SUPERLU_F90 = $(SUPERLU_DIR)FORTRAN/
SUPERLU_LIB = -L$(SUPERLU_DIR)lib/ -lsuperlu_4.1
SUITESPARSE_DIR = $(PROJLIBS)SuiteSparse_libandinclude/
SUITESPARSE_HDR= $(SUITESPARSE_DIR)include/
SUITESPARSE_F90 = $(SUITESPARSE_DIR)F90/
SUITESPARSE_LIB = -L$(SUITESPARSE_DIR)lib/ -lumfpack -L$(SUITESPARSE_DIR)lib/AMD/Lib/ -lamd \
-L$(SUITESPARSE_DIR)lib/ -lcholmod -L$(SUITESPARSE_DIR)lib/ -lcolamd \
-L$(SUITESPARSE_DIR)lib/ -lcamd -L$(SUITESPARSE_DIR)lib/ -lmetis \
-L$(SUITESPARSE_DIR)lib/ -lccolamd



OSTYPE  = $(shell uname -s)
# setting architecture
ARCH   := $(shell uname -m)
# setting hostname
HOSTNAME := $(shell uname -n)

# available Fortran90 compilers for Linux:
NAG      = f95-nag
NAG1     = f95-nag1
LAHEY    = f95-lah
INTEL    = f95-intel
INTEL81  = f95-intel8.1
GNU_F95  = g95
GFORTRAN = gfortran
PORTLAND = f90-pgi

# define standard Fortran90 compiler
ifndef F90
  F90 = $(GFORTRAN)
  #F90 = $(LAHEY)
endif


#NAME_EXT := $(ARCH)
NAME_EXT := $(F90)_$(ARCH)
ifdef DEBUG
  NAME_EXT := DEBUG_$(F90)_$(ARCH)
endif
NAME_EXT := .$(NAME_EXT)
# to avoid an extension to the name of the executable uncomment the next line
#NAME_EXT := 

ifeq ($(OSTYPE), Linux)
##### linux #####
  OUTPUT := $(addsuffix -$(ARCH)-linux, $(OUTPUT))
  CC = gcc
  FC = $(F90)$(VER)
  LD = $(F90)$(VER)
  #CFLAGS = -O3
  #CDEBUG =
  #  C preprocessor defs for compilation for the Fortran interface
  #  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
  #
  CDEFS        = -DAdd_

  #ifeq ($(F90), $(NAG))
  ifeq ($(findstring nag, $(F90)), nag)
    ### NAGf95 ###
    ## disable warnings from license manager
    export NAG_LM_OPTS=nowarn
    ifdef DEBUG
      DEBUGFLAG += -g -C=all #-C=undefined -gline #-mtrace
    endif
    FFLAGS    += -strict95 -maxcontin=500 -v -thread_safe -O2 -nan -u
    LDFLAGS   += #-Bstatic -unsharedf95 -thread_safe
    CUT_ASM_WARN =
    LIBS = -llapack -lblas -lg2c
  endif

  ifeq ($(F90), $(LAHEY))
    ### Lahey ###
    ifdef DEBUG
      DEBUGFLAG += -M OBJS --chk aesux --chkglobal -g --trace --trap # -C=all
      FFLAGS    += -M OBJS --wo --warn --f95 -O --ap --tp4 #--parallel
      CDEBUG    += -g -ggdb -C -p -fbacktrace
      CFLAGS    += -O0
    else
      FFLAGS    += -M OBJS --wo --warn --f95 -O --ap --tp4 #--parallel
      CFLAGS    += "-D DOUBLE_APPEND_FORTRAN"
    endif
    CFLAGS    += "-D DOUBLE_APPEND_FORTRAN"
    LDFLAGS   += # -static
    #LIBS = -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -L/proj/plasma/Libs/ -llapack_amd_lah -lblas_amd_lah -lg2c
    LIBS = $(SUPERLU_LIB) $(LAPACK_LIB) $(BLAS_LIB) # -lg2c
    CUT_ASM_WARN = 2>&1 | grep -v "/tmp/asm"
  endif

  ifeq ($(F90), $(ABSOFT))
    ### Absoft ###
    ifdef DEBUG
      DEBUGFLAG += -g -B80
    endif
    FFLAGS    += -en -v -O #-B100
    LDFLAGS   += 
    CUT_ASM_WARN = 
  endif

  ifeq ($(F90), $(PORTLAND))
    ### Portland ###
    ifdef DEBUG
      DEBUGFLAG += -g
    endif
    FFLAGS    += -v -O2 -tp athlon -Mconcur 
    LDFLAGS   += -Mconcur
    CUT_ASM_WARN = 
  endif

#  ifeq ($(F90), $(INTEL), $(INTEL81))
  ifeq ($(findstring intel,$(F90)), intel)
    ### Intel ###
    ifdef DEBUG
      DEBUGFLAG += -g -inline_debug_info -CB -check all -traceback -DD # -pg
    endif
    FFLAGS    += # -inline_debug_info #-warn all
    LDFLAGS   += #-static -Bstatic 
    CUT_ASM_WARN =
    LIBS = -lefence -llapack -lblas -lg2c
  endif

  ifeq ($(F90), $(GNU_F95))
    ### Gnu F95 ###
    ifdef DEBUG
      DEBUGFLAG += -g -ggdb -C -p -ftrace=full -fbacktrace #-pg #-fmem-report
      FFLAGS    += -Wall
    else
      FFLAGS    += -O3 -Wall 
    endif
    CFLAGS    += "-D DOUBLE_APPEND_FORTRAN"
    LDFLAGS   += #-static  # static linkin does not work on AMD64 bit machines
    CUT_ASM_WARN = 
    # LIBS = -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -llapack -lblas -lg2c
    LIBS = -llapack -lblas -lg2c
  endif

  ifeq ($(F90), $(GFORTRAN))
    ### gfortran Gnu F95 ###

    FFLAGS += -std=gnu #-std=f2003
    ifdef DEBUG
      DEBUGFLAG += -g -ggdb -C -p -pg -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow #-Werror #-pg #-fmem-report
      FFLAGS    += -MOBJS -Waliasing -Wampersand  -Wline-truncation  -Wnonstd-intrinsics  -Wsurprising -Wno-tabs  -Wunderflow #-Wall# -Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation
      CDEBUG    += -g -ggdb -C -p -fbacktrace
      CFLAGS    += -O0 
    else	
      CFLAGS    += -O3 #-Wall 
      FFLAGS    += -J OBJS -O3 -cpp #-Wall
    endif
    CFLAGS    += "-D DOUBLE_APPEND_FORTRAN"
    LDFLAGS   += #-static
    CUT_ASM_WARN = 
    LIBS = $(SUPERLU_LIB) $(SUITESPARSE_LIB) $(LAPACK_LIB) $(BLAS_LIB) # -lg2c
  endif

  ifdef DEBUG    # valgrind --leak-check=yes myprog arg1 arg2
    DBGLIBS += -lefence 
  endif
  
  LIBS += $(DBGLIBS)
  
 #### ifneq (,$(findstring t,$(MAKEFLAGS)))

endif  # ifeq ($(OSTYPE), Linux)




#COMP = f95-lah
#COMP = f95-intel
#COMP = g95
#OPTS = -w 
#OPTS = -M OBJS #-O
#OPTS = -M OBJS --chk a,e,s,u,x
#LIB = -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -llapack -lblas -lg2c
#LIB = -L/usr/lib/gcc-lib/i386-linux/2.95.4/ -llapack -lblas -lg2c
#LIB = -llapack -lblas -lg2c

neo_2.x: $(OBJS) Neo_2.mk
	$(LD) $(DEBUGFLAG) $(LDFLAGS) -o neo_2.x$(NAME_EXT) $(OBJS) $(LIBS)
OBJS/neo_modules.o: neo_modules.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c neo_modules.f90
	mv neo_modules.o OBJS
OBJS/nrutil.o: nrutil.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c nrutil.f90
	mv nrutil.o OBJS
OBJS/sparse_mod.o: sparse_mod.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c sparse_mod.f90
	mv sparse_mod.o OBJS
OBJS/inter_interfaces.o: inter_interfaces.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c inter_interfaces.f90
	mv inter_interfaces.o OBJS
OBJS/solve_system.o: solve_system.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c solve_system.f90
	mv solve_system.o OBJS
OBJS/test_function.o: test_function.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c test_function.f90
	mv test_function.o OBJS
OBJS/spline_cof.o: spline_cof.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c spline_cof.f90
	mv spline_cof.o OBJS
OBJS/spline_int.o: spline_int.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c spline_int.f90
	mv spline_int.o OBJS
OBJS/pspline.o: pspline.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c pspline.f90
	mv pspline.o OBJS
OBJS/neo_sub.o: neo_sub.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c neo_sub.f90
	mv neo_sub.o OBJS
OBJS/neo_magfie.o: neo_magfie.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c neo_magfie.f90
	mv neo_magfie.o OBJS
OBJS/field_divB0.o: field_divB0.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c field_divB0.f90
	mv field_divB0.o OBJS
OBJS/bdivfree.o: bdivfree.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c bdivfree.f90
	mv bdivfree.o OBJS
OBJS/spline5_RZ.o: spline5_RZ.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c spline5_RZ.f90
	mv spline5_RZ.o OBJS
OBJS/compute_aiota.o: compute_aiota.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c compute_aiota.f90
	mv compute_aiota.o OBJS
OBJS/magfie.o: magfie.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c magfie.f90
	mv magfie.o OBJS
OBJS/neo_mod.o: neo_mod.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c neo_mod.f90
	mv neo_mod.o OBJS
OBJS/polleg.o: polleg.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c polleg.f90
	mv polleg.o OBJS
OBJS/collision_operator.o: collision_operator.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c collision_operator.f90
	mv collision_operator.o OBJS
OBJS/lapack_band.o: lapack_band.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c lapack_band.f90
	mv lapack_band.o OBJS
OBJS/binarysplit_mod.o: binarysplit_mod.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c binarysplit_mod.f90
	mv binarysplit_mod.o OBJS
OBJS/propagator.o: propagator.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c propagator.f90
	mv propagator.o OBJS
OBJS/magnetics.o: magnetics.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c magnetics.f90
	mv magnetics.o OBJS
OBJS/device_mod.o: device_mod.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c device_mod.f90
	mv device_mod.o OBJS
OBJS/mag_interface.o: mag_interface.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c mag_interface.f90
	mv mag_interface.o OBJS
OBJS/fluxsplitter.o: fluxsplitter.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c fluxsplitter.f90
	mv fluxsplitter.o OBJS
OBJS/prodband.o: prodband.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c prodband.f90
	mv prodband.o OBJS
OBJS/neo2.o: neo2.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c neo2.f90
	mv neo2.o OBJS
OBJS/flint.o: flint.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c flint.f90
	mv flint.o OBJS
OBJS/mag.o: mag.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c mag.f90
	mv mag.o OBJS
OBJS/kin_allocate.o: kin_allocate.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c kin_allocate.f90
	mv kin_allocate.o OBJS
OBJS/rhs_kin.o: rhs_kin.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c rhs_kin.f90
	mv rhs_kin.o OBJS
OBJS/magdata_for_particles.o: magdata_for_particles.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c magdata_for_particles.f90
	mv magdata_for_particles.o OBJS
OBJS/rk4_kin.o: rk4_kin.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c rk4_kin.f90
	mv rk4_kin.o OBJS
OBJS/ripple_solver.o: ripple_solver.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c ripple_solver.f90
	mv ripple_solver.o OBJS
OBJS/join_ripples_int.o: join_ripples_int.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c join_ripples_int.f90
	mv join_ripples_int.o OBJS
OBJS/join_ripples.o: join_ripples.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c join_ripples.f90
	mv join_ripples.o OBJS
OBJS/join_ends.o: join_ends.f90 Neo_2.mk neo_mod.f90
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c join_ends.f90
	mv join_ends.o OBJS
OBJS/plot_distrf.o: plot_distrf.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c plot_distrf.f90
	mv plot_distrf.o OBJS
OBJS/vvn_tok.o: vvn_tok.f90 Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c vvn_tok.f90
	mv vvn_tok.o OBJS
OBJS/vvn_w7as.o: vvn_w7as.f Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c vvn_w7as.f
	mv vvn_w7as.o OBJS
OBJS/vvn_legendre.o: vvn_legendre.f Neo_2.mk
	$(FC) $(DEBUGFLAG) $(FFLAGS) -c vvn_legendre.f
	mv vvn_legendre.o OBJS
OBJS/c_fortran_dgssv.o: $(SUPERLU_F90)c_fortran_dgssv.c Neo_2.mk
	$(CC) $(CDEBUGFLAG) $(CFLAGS) -I$(SUPERLU_HDR) -c $(SUPERLU_F90)c_fortran_dgssv.c
	mv c_fortran_dgssv.o OBJS
OBJS/c_fortran_zgssv.o: $(SUPERLU_F90)c_fortran_zgssv.c Neo_2.mk
	$(CC) $(CDEBUGFLAG) $(CFLAGS) -I$(SUPERLU_HDR) -c $(SUPERLU_F90)c_fortran_zgssv.c
	mv c_fortran_zgssv.o OBJS
OBJS/umf4_f77wrapper64.o: $(SUITESPARSE_F90)umf4_f77wrapper.c Neo_2.mk
	$(CC) -I$(SUITESPARSE_HDR) -DDLONG -c $(SUITESPARSE_F90)umf4_f77wrapper.c -o umf4_f77wrapper64.o
	mv umf4_f77wrapper64.o OBJS
OBJS/umf4_f77zwrapper64.o: $(SUITESPARSE_F90)umf4_f77zwrapper.c Neo_2.mk
	$(CC) -I$(SUITESPARSE_HDR) -DZLONG -c $(SUITESPARSE_F90)umf4_f77zwrapper.c -o umf4_f77zwrapper64.o
	mv umf4_f77zwrapper64.o OBJS


#.c.o:
#	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)
