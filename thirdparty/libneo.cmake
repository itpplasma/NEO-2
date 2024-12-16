ExternalProject_Add(
    LIBNEO
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG main
    PREFIX ${CMAKE_BINARY_DIR}/libneo
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
