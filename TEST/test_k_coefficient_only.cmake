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
  set(output "fulltransp.h5")
  if(ARGC GREATER 3)
    set(output "${ARGV3}")
  endif()
  set(geometry_required FALSE)
  if(ARGC GREATER 4)
    set(geometry_required "${ARGV4}")
  endif()
  file(MAKE_DIRECTORY "${directory}")
  file(WRITE "${directory}/neo2.in" "${config}")
  if(expected_error STREQUAL "" OR geometry_required)
    file(WRITE "${directory}/neo.in" "${control}")
    file(CREATE_LINK "${FIXTURE}" "${directory}/boozmn_test.nc" SYMBOLIC)
  endif()
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env OMP_NUM_THREADS=1
                  OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "${NEO2_EXE}"
                  WORKING_DIRECTORY "${directory}" TIMEOUT 20
                  RESULT_VARIABLE result OUTPUT_FILE "${directory}/run.log"
                  ERROR_FILE "${directory}/run.log")
  if(expected_error STREQUAL "")
    if(NOT result STREQUAL "0" OR NOT EXISTS "${directory}/${output}")
      message(FATAL_ERROR "${name} failed (${result}); see ${directory}/run.log")
    endif()
  else()
    file(READ "${directory}/run.log" log)
    if(result STREQUAL "0" OR NOT log MATCHES "${expected_error}")
      message(FATAL_ERROR "${name} was not rejected; see ${directory}/run.log")
    endif()
    if(EXISTS "${directory}/${output}")
      message(FATAL_ERROR "${name} published an unconverged K result")
    endif()
  endif()
endfunction()

string(REPLACE "@K_ONLY@" ".false." full_input "${input}")
string(REPLACE "@K_ONLY@" ".true." only_input "${input}")
run_neo(full "${full_input}" "")
run_neo(only "${only_input}" "")
string(REPLACE "niter=100" "niter=3" last_input "${only_input}")
run_neo(laststep "${last_input}" "")
foreach(name only laststep)
  file(READ "${WORK_DIR}/${name}/run.log" convergence_log)
  foreach(source 1 3)
    if(NOT convergence_log MATCHES "Source: +${source} +Number of iterations: +3")
      message(FATAL_ERROR "${name} source ${source} did not converge at iteration 3")
    endif()
  endforeach()
endforeach()
string(REPLACE "niter=100" "niter=1" exhausted_input "${only_input}")
run_neo(exhausted "${exhausted_input}" "integral iteration did not converge"
        "fulltransp.h5" TRUE)
string(REPLACE "niter=100" "niter=0" zero_input "${only_input}")
run_neo(zero_budget "${zero_input}" "requires niter >= 1")
string(REPLACE "niter=100" "niter=-1" negative_input "${only_input}")
run_neo(negative_budget "${negative_input}" "requires niter >= 1")
string(REPLACE "boozer_s=0.5d0" "boozer_s=0.49d0" second_input "${only_input}")
string(REPLACE "conl_over_mfp=1d-2" "conl_over_mfp=2d-2" second_input "${second_input}")
run_neo(second "${second_input}" "")
string(REPLACE "k_coefficient_only=.true." "k_coefficient_only=.true., k_iter_solver=1"
       gmres_input "${only_input}")
string(REPLACE "niter=100" "niter=3" gmres_last_input "${gmres_input}")
run_neo(gmres "${gmres_last_input}" "")
string(REPLACE "k_coefficient_only=.true." "k_coefficient_only=.true., k_iter_solver=1"
       gmres_second_input "${second_input}")
run_neo(gmres_second "${gmres_second_input}" "")
string(REPLACE "niter=100" "niter=2" gmres_exhausted_input "${gmres_input}")
string(REPLACE "epserr_iter=5d-2" "epserr_iter=1d-14"
       gmres_exhausted_input "${gmres_exhausted_input}")
run_neo(exhaustion "${gmres_exhausted_input}" "failed its true-residual"
        "fulltransp.h5" TRUE)
file(GLOB incomplete "${WORK_DIR}/only/efinal.*" "${WORK_DIR}/only/evolve_*.h5"
                     "${WORK_DIR}/only/qf*.dat" "${WORK_DIR}/only/eta_[pm].dat")
if(incomplete)
  message(FATAL_ERROR "K-only emitted transport/distribution files: ${incomplete}")
endif()

string(REPLACE "&settings" "&settings\n k_surface_file='surfaces.dat',"
       batch_input "${only_input}")
set(surface_a "0.5 0 0 0.01\n")
set(surface_b "0.49 0 0 0.02\n")
function(run_batch name rows settings output)
  file(MAKE_DIRECTORY "${WORK_DIR}/${name}")
  file(WRITE "${WORK_DIR}/${name}/surfaces.dat" "${rows}")
  string(REPLACE "&settings" "&settings\n ${settings}" config "${batch_input}")
  set(expected_error "")
  if(ARGC GREATER 4)
    set(expected_error "${ARGV4}")
  endif()
  run_neo("${name}" "${config}" "${expected_error}" "${output}" TRUE)
  if(NOT expected_error STREQUAL "")
    return()
  endif()
  file(GLOB outputs RELATIVE "${WORK_DIR}/${name}" "${WORK_DIR}/${name}/*")
  list(REMOVE_ITEM outputs neo.in neo2.in surfaces.dat run.log boozmn_test.nc)
  if(NOT outputs STREQUAL output)
    message(FATAL_ERROR "${name} emitted unexpected solver output: ${outputs}")
  endif()
endfunction()
run_batch(aba "${surface_a}${surface_b}${surface_a}" "" "k_coefficients.h5")
run_batch(bab "${surface_b}${surface_a}${surface_b}" "" "k_coefficients.h5")
run_batch(chunk "${surface_a}${surface_b}${surface_a}"
          "surface_first=2, surface_last=2, k_output_file='chunk.h5'," "chunk.h5")
set(batch_default "${batch_input}")
string(REPLACE "niter=100" "niter=3" batch_input "${batch_default}")
run_batch(batch_laststep "${surface_a}" "" "k_coefficients.h5")
string(REPLACE "niter=100" "niter=1" batch_input "${batch_default}")
run_batch(batch_exhausted "${surface_a}" "" "k_coefficients.h5")
string(REPLACE "k_coefficient_only=.true." "k_coefficient_only=.true., k_iter_solver=1"
       gmres_batch_input "${batch_default}")
set(batch_input "${gmres_batch_input}")
run_batch(gmres_aba "${surface_a}${surface_b}${surface_a}" "" "k_coefficients.h5")
string(REPLACE "niter=100" "niter=3" batch_input "${gmres_batch_input}")
run_batch(gmres_laststep "${surface_a}" "" "k_coefficients.h5")
file(READ "${WORK_DIR}/gmres_laststep/run.log" gmres_last_log)
foreach(source 1 3)
  if(NOT gmres_last_log MATCHES "K source: +${source} +applications: +3")
    message(FATAL_ERROR "Batch GMRES source ${source} did not use the final application")
  endif()
endforeach()
string(REPLACE "niter=100" "niter=2" batch_input "${gmres_batch_input}")
string(REPLACE "epserr_iter=5d-2" "epserr_iter=1d-14" batch_input "${batch_input}")
run_batch(gmres_exhausted "${surface_a}" "" "k_coefficients.h5"
          "failed its true-residual")
set(batch_input "${batch_default}")
execute_process(COMMAND "${CHECK_EXE}" "${WORK_DIR}" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "K-only numerical comparison failed (${result})")
endif()
execute_process(COMMAND "${CHECK_EXE}" "${WORK_DIR}" "laststep" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Last-step K-only numerical comparison failed (${result})")
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

function(reject_table name rows settings message)
  file(MAKE_DIRECTORY "${WORK_DIR}/${name}")
  file(WRITE "${WORK_DIR}/${name}/surfaces.dat" "${rows}")
  string(REPLACE "&settings" "&settings\n ${settings}" config "${batch_input}")
  run_neo("${name}" "${config}" "${message}")
endfunction()
reject_table(empty "# no surfaces\n" "" "Empty k_surface_file")
reject_table(malformed "0.5 0 0\n" "" "Expected s, R, Z, conl_over_mfp")
reject_table(nonfinite "NaN 0 0 0.01\n" "" "Non-finite surface input")
reject_table(outside "1.1 0 0 0.01\n" "" "Invalid surface s")
reject_table(zero_collision "0.5 0 0 0\n" "" "zero conl_over_mfp")
reject_table(first_zero "${surface_a}" "surface_first=0," "Invalid surface_first")
reject_table(last_zero "${surface_a}" "surface_last=0," "Invalid surface_first")
reject_table(first_past_end "${surface_a}" "surface_first=2," "Invalid surface_first")
reject(iter_solver "k_iter_solver=2" "" "k_iter_solver must be")
reject(iter_solver_negative "k_iter_solver=-1" "" "k_iter_solver must be")
reject(no_integral "k_iter_solver=1, isw_integral=0" ""
       "requires the full integral")
reject(lorentz_only "k_iter_solver=1, isw_integral=1, isw_lorentz=1" ""
       "requires the full integral")
reject(precompute_only "k_iter_solver=1, collop_only_precompute=.true." ""
       "requires the full integral")
foreach(tolerance -1d0 NaN)
  run_neo("tolerance_${tolerance}"
          "&settings epserr_iter=${tolerance} /\n&collision k_coefficient_only=.true., k_iter_solver=1, isw_integral=1, isw_lorentz=0 /\n&binsplit /\n&propagator /\n&plotting /\n"
          "requires (epserr_iter>=0|niter>=2 and a finite tolerance)")
endforeach()
run_neo(no_k_mode
        "&settings /\n&collision k_iter_solver=1 /\n&binsplit /\n&propagator /\n&plotting /\n"
        "requires k_coefficient_only")
run_neo(iter_budget
        "&settings niter=1 /\n&collision k_coefficient_only=.true., k_iter_solver=1, isw_integral=1, isw_lorentz=0 /\n&binsplit /\n&propagator /\n&plotting /\n"
        "requires niter>=2")

function(reject_backend name backend k_only iter_solver method message)
  string(CONCAT config
    "&settings magnetic_device=0, sparse_solve_method=${method} /\n"
    "&collision k_coefficient_only=${k_only}, k_iter_solver=${iter_solver}, "
    "k_sparse_backend=${backend}, isw_axisymm=1, isw_integral=1, isw_lorentz=0 /\n"
    "&binsplit /\n&propagator prop_fileformat=1 /\n&plotting /\n")
  run_neo("${name}" "${config}" "${message}")
  file(GLOB output "${WORK_DIR}/${name}/*.h5")
  if(output)
    message(FATAL_ERROR "${name} wrote output before backend validation: ${output}")
  endif()
endfunction()
reject_backend(backend_negative -1 .true. 1 3 "k_sparse_backend must be")
reject_backend(backend_range 3 .true. 1 3 "k_sparse_backend must be")
reject_backend(backend_no_k 1 .false. 0 3 "requires k_coefficient_only")
reject_backend(backend_no_gmres 1 .true. 0 3 "requires k_iter_solver=1")
reject_backend(backend_method 1 .true. 1 2 "requires sparse_solve_method=3")
foreach(backend 1 2)
  if(NOT NEO2_USE_FORTSPARSE OR
     (backend EQUAL 1 AND NOT FORTSPARSE_ENABLE_SUPERLU) OR
     (backend EQUAL 2 AND NOT FORTSPARSE_ENABLE_MUMPS))
    reject_backend(backend_unavailable_${backend} ${backend} .true. 1 3
                   "k_sparse_backend=${backend} requires fortsparse")
  endif()
endforeach()
