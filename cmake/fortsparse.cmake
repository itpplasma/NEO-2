include(FetchContent)

FetchContent_Declare(
  fortsparse
  GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
  GIT_TAG main
)

set(FORTSPARSE_ENABLE_SUPERLU OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(fortsparse)
