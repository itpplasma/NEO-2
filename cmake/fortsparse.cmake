include(FetchContent)

FetchContent_Declare(
  fortsparse
  GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
  GIT_TAG 23f3f6c62a1632dccd41ddde57957eb322a34f7d
)

set(FORTSPARSE_ENABLE_SUPERLU OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(fortsparse)
