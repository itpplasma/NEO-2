cmake_minimum_required(VERSION 3.22)
project(neo2_krook_codegen LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(gen_krook gen_krook.f90)
target_link_libraries(gen_krook PRIVATE fortsym::fortsym)

add_executable(gen_krook_geometry gen_krook_geometry.f90)
target_link_libraries(gen_krook_geometry PRIVATE fortsym::fortsym)
