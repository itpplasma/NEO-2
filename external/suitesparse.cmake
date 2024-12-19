FetchContent_Declare(
    SuiteSparse
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    GIT_REPOSITORY https://github.com/DrTimothyAldenDavis/SuiteSparse.git
    GIT_TAG v7.6.1
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(SuiteSparse)
