execute_process(COMMAND "${EXECUTABLE}" "${ARGUMENT}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 5)
if(status EQUAL 0 OR status MATCHES "timeout")
    message(FATAL_ERROR "Expected immediate sparse configuration rejection: ${status}")
endif()
set(receipt "${output}\n${errors}")
if(NOT receipt MATCHES "${EXPECTED}")
    message(FATAL_ERROR "Wrong sparse configuration diagnostic: ${receipt}")
endif()
