# Can also run without the full NEO-2 dependency tree:
# cmake -S TEST/moment_reduction -B build/moment -DCMAKE_BUILD_TYPE=Debug
cmake_minimum_required(VERSION 3.18)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(NEO2_MOMENT_TESTS LANGUAGES Fortran)
    enable_testing()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()
get_filename_component(NEO2_MOMENT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(MOMENT_MOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules")
file(MAKE_DIRECTORY "${MOMENT_MOD_DIR}")

# Do not duplicate the legacy implementation in a fixture: test the current
# production routine verbatim, including its actual signs and boundary terms.
set(NATIVE_SOURCE "${NEO2_MOMENT_ROOT}/NEO-2-QL/ripple_solver_ArnoldiOrder2_test.f90")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${NATIVE_SOURCE}")
file(READ "${NATIVE_SOURCE}" NATIVE_TEXT)
set(BEGIN_MARKER "    subroutine integral_part(vec_in, vec_out)")
set(END_MARKER "    end subroutine integral_part")
string(FIND "${NATIVE_TEXT}" "${BEGIN_MARKER}" BEGIN_INDEX)
if(BEGIN_INDEX LESS 0)
    message(FATAL_ERROR "Cannot extract native integral_part: review the oracle adapter")
endif()
string(SUBSTRING "${NATIVE_TEXT}" ${BEGIN_INDEX} -1 NATIVE_TAIL)
string(FIND "${NATIVE_TAIL}" "${END_MARKER}" END_INDEX)
if(END_INDEX LESS 0)
    message(FATAL_ERROR "Cannot find end of native integral_part")
endif()
string(LENGTH "${END_MARKER}" END_LENGTH)
math(EXPR NATIVE_LENGTH "${END_INDEX}+${END_LENGTH}")
string(SUBSTRING "${NATIVE_TAIL}" 0 ${NATIVE_LENGTH} NATIVE_ROUTINE)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/native_integral_part.inc" "${NATIVE_ROUTINE}\n")
string(SHA256 NATIVE_SHA "${NATIVE_ROUTINE}\n")
message(STATUS "Native integral_part oracle SHA256: ${NATIVE_SHA}")

# Separate module directory prevents .mod collisions with the full QL build.
add_library(neo2_moment_test_core STATIC
    "${NEO2_MOMENT_ROOT}/NEO-2-QL/collision_moment_map.f90"
    "${NEO2_MOMENT_ROOT}/NEO-2-QL/moment_schur_solver.f90"
    test_support.f90)
set_target_properties(neo2_moment_test_core PROPERTIES
    Fortran_MODULE_DIRECTORY "${MOMENT_MOD_DIR}")
target_include_directories(neo2_moment_test_core PUBLIC "${MOMENT_MOD_DIR}")
target_link_libraries(neo2_moment_test_core PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(test_collision_moment_map native_integral_oracle.f90 test_collision_moment_map.f90)
add_executable(test_moment_schur test_moment_schur.f90)
foreach(TARGET test_collision_moment_map test_moment_schur)
    target_link_libraries(${TARGET} PRIVATE neo2_moment_test_core)
    set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY "${MOMENT_MOD_DIR}")
    target_include_directories(${TARGET} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(TARGET neo2_moment_test_core test_collision_moment_map test_moment_schur)
        target_compile_options(${TARGET} PRIVATE -ffree-line-length-none
            "$<$<CONFIG:Debug>:-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>")
    endforeach()
endif()
add_test(NAME ql_moment_factorization COMMAND test_collision_moment_map)
add_test(NAME ql_moment_schur COMMAND test_moment_schur)
set_tests_properties(ql_moment_factorization ql_moment_schur PROPERTIES
    PASS_REGULAR_EXPRESSION "All tests passed!"
    FAIL_REGULAR_EXPRESSION "FAIL"
    TIMEOUT 60
    ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")

# Opt-in CPU microbenchmark (NOT CTest: timed tests are not CI correctness gates).
add_executable(bench_ql_moment_vs_full bench_moment_vs_full.f90)
target_link_libraries(bench_ql_moment_vs_full PRIVATE neo2_moment_test_core)
set_target_properties(bench_ql_moment_vs_full PROPERTIES Fortran_MODULE_DIRECTORY "${MOMENT_MOD_DIR}")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(bench_ql_moment_vs_full PRIVATE -ffree-line-length-none)
endif()
