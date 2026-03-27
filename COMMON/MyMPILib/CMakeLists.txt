find_package(MPI REQUIRED COMPONENTS Fortran)

message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -pthread")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
endif ()

if (MPE_SUPPORT)
  message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
endif ()

add_library(NEO2MyMPILib
  ${CMAKE_CURRENT_SOURCE_DIR}/Specific/packBuffer_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Specific/mpiprovider_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/list_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/packable_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/matrix_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/clientStatus_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Generic/workunit_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/wuDataRequester_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/intList_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/wuList_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/packableList_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/commandline_parser_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Generic/scheduler_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/mpelog_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/myLog_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Tools/configFile_parser_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/wuMergeWorkunit_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Internal/wuMergeChunk_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Generic/initWorkunit_module.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Generic/genericWorkunit_module.f90
)
set_target_properties(NEO2MyMPILib PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
  OUTPUT_NAME MyMPILib
)
target_compile_definitions(NEO2MyMPILib PRIVATE $<$<BOOL:${MPE_SUPPORT}>:MPE_SUPPORT>)
target_include_directories(NEO2MyMPILib
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/modules
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Generic
    ${CMAKE_CURRENT_SOURCE_DIR}/Specific
    ${CMAKE_CURRENT_SOURCE_DIR}/Internal
    ${CMAKE_CURRENT_SOURCE_DIR}/Tools
    $<$<BOOL:${MPE_SUPPORT}>:${MPE_PATH}/include>
)
target_link_libraries(NEO2MyMPILib PUBLIC MPI::MPI_Fortran)

add_custom_target(
  MyMPILib-Versioning
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/do_versioning.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_dependencies(NEO2MyMPILib MyMPILib-Versioning)
