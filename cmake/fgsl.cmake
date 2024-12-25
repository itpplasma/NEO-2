include(ExternalProject)

ExternalProject_Add(
    FGSL
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX ${CMAKE_BINARY_DIR}/fgsl
    URL https://github.com/reinh-bader/fgsl/archive/refs/tags/v1.6.0.tar.gz
    CONFIGURE_COMMAND mkdir m4 && autoreconf -i && ./configure
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS
       ${CMAKE_BINARY_DIR}/fgsl/src/FGSL/.libs/libfgsl${CMAKE_STATIC_LIBRARY_SUFFIX}
)

add_custom_target(external_fgsl DEPENDS FGSL)
add_dependencies(external_fgsl FGSL)
