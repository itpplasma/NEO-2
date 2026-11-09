include(FetchContent)

FetchContent_Declare(
    SuiteSparse
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse.git
    GIT_TAG v7.6.1
    GIT_SHALLOW TRUE
    OVERRIDE_FIND_PACKAGE TRUE
)

# NEO-2 only links UMFPACK. SuiteSparse defaults to building every project,
# including the large GraphBLAS C/C++ tree, which is unrelated to this build
# and can introduce platform-specific toolchain failures.
set(SUITESPARSE_ENABLE_PROJECTS "umfpack" CACHE STRING
    "SuiteSparse projects to build")

FetchContent_MakeAvailable(SuiteSparse)
