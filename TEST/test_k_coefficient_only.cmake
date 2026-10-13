# Exercise the actual binary with the existing axisymmetric geometry fixture.
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(READ "${CONTROL_TEMPLATE}" control)
foreach(setting "in_file|boozmn_test.nc" "inp_swi|10" "theta_n|32"
                "phi_n|16" "write_progress|0")
  string(REPLACE "|" ";" fields "${setting}")
  list(GET fields 0 key)
  list(GET fields 1 value)
  string(REGEX REPLACE "[^\n]*[ \t]${key}[ \t][^\n]*" "${value} ${key}"
         control "${control}")
endforeach()

set(input [=[
&settings
 nstep=80, nperiod=2, mag_nperiod_min=1,
 mag_magfield=1, magnetic_device=0, mag_coordinates=1, boozer_s=0.5d0,
 mag_close_fieldline=2, sparse_solve_method=3,
 epserr_iter=5d-2, niter=100, epserr_sink=1d-5
/
&collision
 k_coefficient_only=@K_ONLY@, isw_axisymm=1,
 conl_over_mfp=1d-2, lag=2, leg=2, legmax=2, z_eff=0d0,
 isw_lorentz=0, isw_integral=1, isw_energy=1
/
&binsplit
 eta_part=16, bin_split_mode=0, phi_split_mode=2, max_solver_try=10
/
&propagator
 prop_fileformat=1, prop_reconstruct=0, prop_ripple_plot=0, prop_write=2,
 prop_diagnostic=0, prop_diagphys=0, mag_talk=.false., mag_infotalk=.false.
/
&plotting /
]=])

function(run_neo name config expected_error)
  set(directory "${WORK_DIR}/${name}")
  set(threads 1)
  if(ARGC GREATER 3)
    set(threads "${ARGV3}")
  endif()
  file(MAKE_DIRECTORY "${directory}")
  file(WRITE "${directory}/neo2.in" "${config}")
  if(expected_error STREQUAL "")
    file(WRITE "${directory}/neo.in" "${control}")
    file(CREATE_LINK "${FIXTURE}" "${directory}/boozmn_test.nc" SYMBOLIC)
  endif()
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env "OMP_NUM_THREADS=${threads}"
                  OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "${NEO2_EXE}"
                  WORKING_DIRECTORY "${directory}" TIMEOUT 20
                  RESULT_VARIABLE result OUTPUT_FILE "${directory}/run.log"
                  ERROR_FILE "${directory}/run.log")
  if(expected_error STREQUAL "")
    if(NOT result STREQUAL "0" OR NOT EXISTS "${directory}/fulltransp.h5")
      message(FATAL_ERROR "${name} failed (${result}); see ${directory}/run.log")
    endif()
  else()
    file(READ "${directory}/run.log" log)
    if(result STREQUAL "0" OR NOT log MATCHES "${expected_error}")
      message(FATAL_ERROR "${name} was not rejected before geometry; see ${directory}/run.log")
    endif()
  endif()
endfunction()

string(REPLACE "@K_ONLY@" ".false." full_input "${input}")
string(REPLACE "@K_ONLY@" ".true." only_input "${input}")
run_neo(full "${full_input}" "")
run_neo(only "${only_input}" "")
run_neo(full2 "${full_input}" "" 2)
run_neo(only2 "${only_input}" "" 2)
file(GLOB incomplete "${WORK_DIR}/only/efinal.*" "${WORK_DIR}/only/evolve_*.h5"
                     "${WORK_DIR}/only/qf*.dat" "${WORK_DIR}/only/eta_[pm].dat"
                     "${WORK_DIR}/only2/efinal.*" "${WORK_DIR}/only2/evolve_*.h5"
                     "${WORK_DIR}/only2/qf*.dat" "${WORK_DIR}/only2/eta_[pm].dat")
if(incomplete)
  message(FATAL_ERROR "K-only emitted transport/distribution files: ${incomplete}")
endif()
execute_process(COMMAND "${CHECK_EXE}" "${WORK_DIR}" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "K-only numerical comparison failed (${result})")
endif()
execute_process(COMMAND "${CHECK_EXE}" "${WORK_DIR}" "2" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Two-thread K-only numerical comparison failed (${result})")
endif()
execute_process(COMMAND "${CHECK_EXE}" "${WORK_DIR}" "2" "2" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Two-thread default-mode numerical comparison failed (${result})")
endif()

# These cases deliberately have no neo.in or geometry file.
function(reject name collision propagator message)
  run_neo("${name}" "&settings magnetic_device=0 /\n&collision k_coefficient_only=.true., ${collision} /\n&binsplit /\n&propagator ${propagator} /\n&plotting /\n" "${message}")
endfunction()
reject(axisymmetric "isw_axisymm=0" "prop_fileformat=1" "single-species axisymmetric")
reject(multispecies "isw_axisymm=1, lsw_multispecies=.true." "prop_fileformat=1"
       "single-species axisymmetric")
reject(momentum "isw_axisymm=1, isw_momentum=1" "prop_fileformat=1" "Laguerre momentum")
reject(ascii "isw_axisymm=1" "prop_fileformat=0" "requires HDF5")
reject(reconstruction "isw_axisymm=1" "prop_fileformat=1, prop_reconstruct=1" "requires HDF5")
reject(plot "isw_axisymm=1" "prop_fileformat=1, prop_ripple_plot=1" "does not support")
