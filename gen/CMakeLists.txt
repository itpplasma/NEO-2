cmake_minimum_required(VERSION 3.22)
# Standalone generator project: NEO-2 itself builds from the committed output
# and never needs fortsym.
project(neo2_codegen LANGUAGES Fortran)
find_package(fortsym CONFIG REQUIRED)
add_executable(gen_phi_ma gen_phi_ma.f90)
target_link_libraries(gen_phi_ma PRIVATE fortsym::fortsym)
