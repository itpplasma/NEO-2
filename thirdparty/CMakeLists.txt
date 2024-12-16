cmake_minimum_required(VERSION 3.18)
project(FGSLProject)

include(ExternalProject)

ExternalProject_Add(
    FGSL
    URL https://github.com/reinh-bader/fgsl/archive/refs/tags/v1.6.0.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/fgsl
    CONFIGURE_COMMAND mkdir m4 && autoreconf -i && ./configure
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
