# AMG preconditioner library for NEO-2

# Add AMG source files
set(AMG_SOURCES
  amg_types_mod.f90
  amg_smoothers_mod.f90
  amg_smoothed_aggregation_mod.f90
  amg_cycles_mod.f90
  amg_precond_mod.f90
)

# Create AMG library
add_library(neo2_amg STATIC ${AMG_SOURCES})

# Set Fortran module directory
set_target_properties(neo2_amg PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# Include directories for module files
target_include_directories(neo2_amg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include>
)

# Add dependency on sparse_types library
add_dependencies(neo2_amg sparse_types)

# Link dependencies
target_link_libraries(neo2_amg PUBLIC
  sparse_types
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
)

# Include the build directory where sparse_types_mod.mod will be
target_include_directories(neo2_amg PUBLIC
  ${CMAKE_BINARY_DIR}/COMMON
)

# Install targets
install(TARGETS neo2_amg
  EXPORT NEO2Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Install Fortran module files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
  DESTINATION include
  FILES_MATCHING PATTERN "*.mod"
)