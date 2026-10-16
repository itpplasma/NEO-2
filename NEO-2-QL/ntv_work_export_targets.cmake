add_library(ntv_work_export STATIC
  ${CMAKE_CURRENT_LIST_DIR}/ntv_work_export.f90)
set_target_properties(ntv_work_export PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/work_modules)
target_include_directories(ntv_work_export PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/work_modules)
