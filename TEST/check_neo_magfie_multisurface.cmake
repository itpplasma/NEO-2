if(NOT DEFINED BINARY OR NOT DEFINED FIXTURE OR NOT DEFINED WORKROOT OR NOT DEFINED FIELD)
  message(FATAL_ERROR "BINARY, FIXTURE, FIELD and WORKROOT are required")
endif()
file(MAKE_DIRECTORY "${WORKROOT}")
file(COPY "${FIXTURE}/neo.in" DESTINATION "${WORKROOT}")
configure_file("${FIELD}" "${WORKROOT}/axi.bc" COPYONLY)
foreach(mode one_first one_second both compare)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    "${BINARY}" "${mode}" WORKING_DIRECTORY "${WORKROOT}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Native ${mode} failed (${status}): ${output} ${error}")
  endif()
endforeach()
message("All tests passed!")
