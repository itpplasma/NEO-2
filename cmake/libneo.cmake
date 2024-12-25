include(FetchContent)

FetchContent_Declare(
    LIBNEO
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG 14-user-friendly-installation
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/libneo${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/libmagfie${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/src/MyMPILib/libMyMPILib${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/hdf5_tools/libhdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX}
    OVERRIDE_FIND_PACKAGE TRUE
)

FetchContent_MakeAvailable(LIBNEO)
