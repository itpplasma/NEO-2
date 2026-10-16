cmake_minimum_required(VERSION 3.24)
project(neo2_field_codegen LANGUAGES Fortran C CXX)
include(FetchContent)
set(FORTSYM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fortsym
    GIT_REPOSITORY https://github.com/lazy-fortran/fortsym.git
    GIT_TAG 7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6)
FetchContent_MakeAvailable(fortsym)
add_executable(gen_covariant_curl gen_covariant_curl.f90)
target_link_libraries(gen_covariant_curl PRIVATE fortsym)
