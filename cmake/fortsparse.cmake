include(FetchContent)

FetchContent_Declare(
  fortsparse
  GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
  GIT_TAG main
)

FetchContent_MakeAvailable(fortsparse)
