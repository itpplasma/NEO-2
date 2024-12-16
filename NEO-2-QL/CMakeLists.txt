cmake_minimum_required (VERSION 3.18)

# Enable project
project(NEO-2-QL)
enable_testing()

### Check Build-Type and set Build-Name extension

set(EXE_NAME "neo_2_ql.x")

# Print the current compiler
message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

### Define flags for Fortran/C compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -ffpe-summary=invalid,zero,overflow")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  message(STATUS "Setting Intel flags.")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -cpp -mkl -heap-arrays -stand f08")
endif ()

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(MY_DEBUG_FLAG_LIST
  " -O0 -cpp -g"
 " -ffpe-summary=invalid,zero,overflow "
 " -fcheck=all "
 " -fbacktrace "
 " -Wall "
 " -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-label "
 " -Wshadow "
  )
  string(REPLACE ";" "" MY_DEBUG_FLAG  ${MY_DEBUG_FLAG_LIST})
  set(CMAKE_Fortran_FLAGS_DEBUG "${MY_DEBUG_FLAG}")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG "-cpp -mkl -g -warn all -stand f08 ")
endif ()
set(CFLAGS "-O2 -DDOUBLE_APPEND_FORTRAN")

### Define container for *.o and *.mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

### Source files
include(CMakeSources.in)

### Define executable
add_executable(${EXE_NAME}
  ${NEO2_SRC_FILES}
  ${COLLOPLIB_FILES}
)

### Add libraries linked with executable
target_link_libraries(${EXE_NAME} common)
target_include_directories(${EXE_NAME} PUBLIC "$<TARGET_PROPERTY:common,Fortran_MODULE_DIRECTORY>")

### Specify command "make copy"
add_custom_target(copy
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/${EXE_NAME}
  ${PROJECT_SOURCE_DIR}/../../TEST_RUNS/${EXE_NAME}
  )

### Testing section
ADD_TEST(Test_Build_Code "make")
ADD_TEST(Extended_QL_Multispecies  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies" "2" "QL")
ADD_TEST(Extended_QL_Multispecies_Read  ${PROJECT_SOURCE_DIR}/../ShellScripts/test_full_reconstruct.sh "Extended_QL_Multispecies_Read" "2" "QL")
#ADD_TEST(Extended_QL_Multispecies  ${CMAKE_SOURCE_DIR}/../ShellScripts/test_scan.sh "Extended_QL_Scan_Kink" "QL")
