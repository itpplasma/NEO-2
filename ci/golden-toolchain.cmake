include("${CMAKE_CURRENT_LIST_DIR}/cluster-toolchain.cmake")
# Preserve the golden workflow's original project compiler flags. The unit
# tier's initial flags would change the checked legacy CMake_Flags metadata.
set(CMAKE_Fortran_FLAGS_INIT "")
