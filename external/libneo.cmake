FetchContent_Declare(
    LIBNEO
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG 14-user-friendly-installation
    BUILD_BYPRODUCTS
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/src/LIBNEO/build/libneo${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/libneo/src/LIBNEO/build/libmagfie${CMAKE_SHARED_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/libneo/src/LIBNEO/build/src/MyMPILib/libMyMPILib${CMAKE_STATIC_LIBRARY_SUFFIX}
        ${CMAKE_BINARY_DIR}/_deps/libneo-build/libneo/src/LIBNEO/build/src/hdf5_tools/libhdf5_tools${CMAKE_STATIC_LIBRARY_SUFFIX}
)

FetchContent_MakeAvailable(LIBNEO)
