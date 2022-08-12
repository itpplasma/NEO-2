### Initialize CMake (some internal stuff)
cmake_minimum_required (VERSION 2.8.9)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMAKE_TOOLS")

project (neo-2-common)
enable_language (Fortran)
enable_testing()

### Macro for comparison of version strings
include(CMAKE_TOOLS/MacroCompareVersionStrings.cmake)

### SUPERLU ###
set(SUPERLU_SRC_FILES_VER_4_1
        c_fortran_dgssv.c
        c_fortran_zgssv.c
)
set(SUPERLU_FLAGS "${CFLAGS_DEBUG} ${CFLAGS} -I${SUPERLU_INCLUDES}")
set_source_files_properties(${SUPERLU_SRC_FILES_VER_4_1} PROPERTIES COMPILE_FLAGS ${SUPERLU_FLAGS})

### SUPERLU Version 5.0 ###
set(SUPERLU_SRC_FILES_VER_5_0
        c_fortran_dgssv_ver_5_0.c
        c_fortran_zgssv_ver_5_0.c
)
set(SUPERLU_FLAGS "${CFLAGS_DEBUG} ${CFLAGS} -I${SUPERLU_INCLUDES}")
set_source_files_properties(${SUPERLU_SRC_FILES_VER_5_0} PROPERTIES COMPILE_FLAGS ${SUPERLU_FLAGS})


### SUITESPARSE ###
set(SUITESPARSE_SRC_FILES_VER_3_6
        ../COMMON/umf4_f77wrapper.c
        ../COMMON/umf4_f77zwrapper.c
)
set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
set_source_files_properties(../COMMON/umf4_f77wrapper.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(../COMMON/umf4_f77zwrapper.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")

### SUITESPARSE Version 4.5 ###
set(SUITESPARSE_SRC_FILES_VER_4_5
        ../COMMON/umf4_f77wrapper_ver_4_5.c
        ../COMMON/umf4_f77zwrapper_ver_4_5.c
)
set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
set_source_files_properties(../COMMON/umf4_f77wrapper_ver_4_5.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(../COMMON/umf4_f77zwrapper_ver_4_5.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")

### FGSL user-defined interface ###
set(FGSL_USERINTERFACE_SOURCE_FILES
      gsl_roots_routines_mod.f90
      gsl_integration_routines_mod.f90
      gsl_bspline_routines_mod.f90
      gsl_specialfunctions_mod.f90
)
set(FGSL_USERINTERFACE_FLAGS "-I${FGSL_INC}")
set_source_files_properties(${FGSL_USERINTERFACE_SOURCE_FILES} PROPERTIES COMPILE_FLAGS ${FGSL_USERINTERFACE_FLAGS})

set(COMMON_FILES
  binarysplit_int.f90
  binarysplit_mod.f90
  collisionality_mod.f90
  collision_operator_mems.f90
  collop_bspline.f90
  collop_compute.f90
  collop_laguerre.f90
  collop_nbi.f90
  collop_polynomial.f90
  collop_spline.f90
  compute_aiota.f90
  development.f90
  device_mod.f90
  flint_mod.f90
  fluxsplitter.f90
  inter_interfaces.f90
  inter_precision.f90
  kin_allocate.f90
  lapack_band.f90
  magdata_for_particles.f90
  mag.f90
  magfie.f90
  mag_interface.f90
  magnetics.f90
  neo_actual_fluxs.f90
  neo_actual_spectra.f90
  neo_control.f90
  neo_eval_switch.f90
  neo_exchange.f90
  neo_input.f90
  neo_magfie.f90
  neo_parameters.f90
  neo_precision.f90
  neo_spline_b00.f90
  neo_spline_data.f90
  neo_spline.f90
  neo_sub.f90
  neo_support.f90
  neo_units.f90
  neo_van.f90
  neo_work.f90
  nrtype.f90
  nrutil.f90
  odeint_allroutines.f
  partpa_mod.f90
  plagrange_mod.f90
  plot_distrf.f90
  polleg.f90
  prodband.f90
  rel_kernels.f90
  rhs_kin.f90
  rk4_kin.f90
  rk4_kin_mod.f90
  rkstep_mod.f90
  size_mod.f90
  sizey_bo.f90
  sizey_cur.f90
  sizey_pla.f90
  solve_system.f90
  sparse_mod.f90
  sparsevec_mod.f90
  spline_cof.f90
  spline_int.f90
  spline_mod.f90
  test_function.f90
  vvn_legendre.f
  vvn_tok.f90
  vvn_w7as.f
  )
set_source_files_properties(${COMMON_FILES} PROPERTIES COMPILE_FLAGS "${DEBUGFLAG} ${FFLAG_DEBUG} ${FFLAGS} ${MPIFLAGS}")

### SuperLU
find_package(SuperLU REQUIRED)

### SuiteSparse and METIS - ignore if already found.
if (SUITESPARSEQR_FOUND)
else (SUITESPARSEQR_FOUND)
  find_package(SuiteSparse REQUIRED)
endif (SUITESPARSEQR_FOUND)

add_library(common STATIC
  ${COMMON_FILES}
  ${FGSL_USERINTERFACE_SOURCE_FILES}
  )

if ($ENV{COMMON_UNIT_TESTS}x STREQUAL "TRUEx")
  find_package(PFUNIT REQUIRED)
  enable_testing()
  add_subdirectory(unit_tests)
endif()
