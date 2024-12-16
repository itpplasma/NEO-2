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

# Optionally, you can add include directories and link libraries
# to your project if you need to use FGSL in your code.
# For example:
# include_directories(${CMAKE_BINARY_DIR}/fgsl/src/FGSL/include)
# target_link_libraries(your_target ${CMAKE_BINARY_DIR}/fgsl/src/FGSL/.libs/libfgsl.a)
