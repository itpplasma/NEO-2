project(NEO-2-COMMON)

### FGSL user-defined interface
set(FGSL_USERINTERFACE_SOURCE_FILES
      gsl_roots_routines_mod.f90
      gsl_integration_routines_mod.f90
      gsl_bspline_routines_mod.f90
      gsl_specialfunctions_mod.f90
)

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
  sparse_types_mod.f90
  sparse_conversion_mod.f90
  sparse_io_mod.f90
  sparse_arithmetic_mod.f90
  sparse_solvers_mod.f90
  sparse_utils_mod.f90
  sparse_mod.f90
  sparsevec_mod.f90
  spline_cof.f90
  splinecof3_direct_sparse.f90
  spline_int.f90
  spline_mod.f90
  test_function.f90
  vvn_legendre.f
  vvn_tok.f90
  vvn_w7as.f
  )
set_source_files_properties(${COMMON_FILES} PROPERTIES COMPILE_FLAGS "${DEBUGFLAG} ${FFLAG_DEBUG} ${FFLAGS} ${MPIFLAGS}")


set(SUITESPARSE_SRC_FILES
  umf4_f77wrapper_ver_4_5.c
  umf4_f77zwrapper_ver_4_5.c
)
if(LEGACY_SUITESPARSE)
  set(SUITESPARSE_FLAGS "-I${SUITESPARSE_INCLUDE_DIRS}")
else()
  set(SUITESPARSE_FLAGS "")
endif()
set_source_files_properties(umf4_f77wrapper_ver_4_5.c  PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DDLONG")
set_source_files_properties(umf4_f77zwrapper_ver_4_5.c PROPERTIES COMPILE_FLAGS "${SUITESPARSE_FLAGS} -DZLONG")

add_library(common STATIC
  ${COMMON_FILES}
  ${SUITESPARSE_SRC_FILES}
  ${FGSL_USERINTERFACE_SOURCE_FILES}
  )
add_dependencies(common FGSL)
target_link_libraries(common PUBLIC
  ${LIBFGSL}
  GSL::gsl
  ${MPE_lib}
  ${MPI_Fortran_LIBRARIES}
  BLAS::BLAS
  LAPACK::LAPACK
  ${NETCDF_LIBS}
  -lz -lpthread
)
if(LEGACY_SUITESPARSE)
  target_link_libraries(common PUBLIC ${SUITESPARSE_LIBRARIES})
else()
  target_link_libraries(common PUBLIC SuiteSparse::umfpack)
endif()

target_link_libraries(common PUBLIC
  LIBNEO::neo
  LIBNEO::magfie
  LIBNEO::MyMPILib
  LIBNEO::hdf5_tools
)

target_link_libraries(common PUBLIC
  MPI::MPI_Fortran
)

if ($ENV{COMMON_UNIT_TESTS}x STREQUAL "TRUEx")
  find_package(PFUNIT REQUIRED)
  enable_testing()
  add_subdirectory(unit_tests)
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/cmake_version.f90.in"
  "${PROJECT_BINARY_DIR}/cmake_version.f90"
  )

add_custom_target(
    NEO-2-Versioning
    COMMAND ${PROJECT_SOURCE_DIR}/../ShellScripts/tag_version.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)
add_dependencies(common NEO-2-Versioning)
