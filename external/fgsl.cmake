FetchContent_Declare(
    FGSL
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL "https://github.com/reinh-bader/fgsl/archive/refs/tags/v1.6.0.tar.gz"
    CONFIGURE_COMMAND "mkdir m4 && autoreconf -i && ./configure"
    BUILD_COMMAND make
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/_deps/fgsl-build/src/FGSL/.libs/libfgsl${CMAKE_STATIC_LIBRARY_SUFFIX}
)

FetchContent_MakeAvailable(FGSL)
