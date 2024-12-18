ExternalProject_Add(
    SUITESPARSE
    PREFIX ${CMAKE_BINARY_DIR}/suitesparse
    URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.8.3.tar.gz
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)
