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
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/libneo${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/libmagfie${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/src/MyMPILib/libMyMPILib${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libneo/src/LIBNEO/build/src/hdf5_tools/libhdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX}
)
