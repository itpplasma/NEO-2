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
  set(geometry_required FALSE)
  if(ARGC GREATER 3)
    set(geometry_required "${ARGV3}")
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
    if(NOT result STREQUAL "0" OR NOT EXISTS "${directory}/fulltransp.h5")
      message(FATAL_ERROR "${name} failed (${result}); see ${directory}/run.log")
    endif()
  else()
    file(READ "${directory}/run.log" log)
    if(result STREQUAL "0" OR NOT log MATCHES "${expected_error}")
      message(FATAL_ERROR "${name} was not rejected; see ${directory}/run.log")
    endif()
    if(EXISTS "${directory}/fulltransp.h5")
      message(FATAL_ERROR "${name} published an unconverged K result")
    endif()
  endif()
endfunction()

string(REPLACE "@K_ONLY@" ".false." full_input "${input}")
string(REPLACE "@K_ONLY@" ".true." only_input "${input}")
run_neo(full "${full_input}" "")
run_neo(only "${only_input}" "")
# Rerun with the converged iteration count as budget: convergence on the last
# allowed step must be accepted.
file(READ "${WORK_DIR}/only/run.log" only_log)
string(REGEX MATCHALL "Number of iterations: +[0-9]+" counts "${only_log}")
set(last_iter 0)
foreach(count IN LISTS counts)
  string(REGEX REPLACE ".* " "" count "${count}")
  if(count GREATER last_iter)
    set(last_iter ${count})
  endif()
endforeach()
if(last_iter LESS 2)
  message(FATAL_ERROR "Fixture must need at least two integral iterations")
endif()
string(REPLACE "niter=100" "niter=${last_iter}" last_input "${only_input}")
run_neo(laststep "${last_input}" "")
string(REPLACE "niter=100" "niter=1" exhausted_input "${only_input}")
run_neo(exhausted "${exhausted_input}" "integral iteration did not converge" TRUE)
string(REPLACE "k_coefficient_only=.true." "k_coefficient_only=.true., k_iter_solver=1"
       gmres_input "${only_input}")
# GMRES needs a platform-dependent number of operator applications on this
# ill-conditioned fixture. Measure it with a free budget, then rerun with the
# budget that makes convergence fall on the final allowed application.
run_neo(gmres_free "${gmres_input}" "")
file(READ "${WORK_DIR}/gmres_free/run.log" gmres_log)
string(REGEX MATCHALL "applications: +[0-9]+" counts "${gmres_log}")
set(free_apps 0)
foreach(count IN LISTS counts)
  string(REGEX REPLACE ".* " "" count "${count}")
  if(count GREATER free_apps)
    set(free_apps ${count})
  endif()
endforeach()
if(free_apps LESS 3)
  message(FATAL_ERROR "Free GMRES run must take at least one Arnoldi step")
endif()
# The free run spends one extra application confirming the converged iterate.
math(EXPR gmres_budget "${free_apps} - 1")
string(REPLACE "niter=100" "niter=${gmres_budget}" gmres_last_input "${gmres_input}")
run_neo(gmres "${gmres_last_input}" "")
string(REPLACE "niter=100" "niter=2" gmres_exhausted_input "${gmres_input}")
string(REPLACE "epserr_iter=5d-2" "epserr_iter=1d-14"
       gmres_exhausted_input "${gmres_exhausted_input}")
run_neo(exhaustion "${gmres_exhausted_input}" "failed its true-residual" TRUE)
if(NEO2_USE_FORTSPARSE)
  string(REPLACE "k_iter_solver=1" "k_iter_solver=1, k_sparse_backend=2"
         mumps_input "${gmres_input}")
  run_neo(mumps "${mumps_input}" "")
endif()
file(GLOB incomplete "${WORK_DIR}/only/efinal.*" "${WORK_DIR}/only/evolve_*.h5"
                     "${WORK_DIR}/only/qf*.dat" "${WORK_DIR}/only/eta_[pm].dat")
if(incomplete)
  message(FATAL_ERROR "K-only emitted transport/distribution files: ${incomplete}")
endif()
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
reject(iter_solver "k_iter_solver=2" "" "k_iter_solver must be")
reject(no_integral "k_iter_solver=1, isw_integral=0" ""
       "requires the full integral")
reject(backend_range "k_iter_solver=1, k_sparse_backend=1" "" "k_sparse_backend must be")
reject(backend_no_gmres "k_sparse_backend=2" "" "requires k_iter_solver=1")
if(NOT NEO2_USE_FORTSPARSE)
  reject(backend_unavailable "k_iter_solver=1, k_sparse_backend=2" ""
         "requires NEO2_USE_FORTSPARSE")
endif()
run_neo(no_k_mode
        "&settings /\n&collision k_iter_solver=1 /\n&binsplit /\n&propagator /\n&plotting /\n"
        "requires k_coefficient_only")
run_neo(iter_budget
        "&settings niter=1 /\n&collision k_coefficient_only=.true., k_iter_solver=1, isw_integral=1, isw_lorentz=0 /\n&binsplit /\n&propagator /\n&plotting /\n"
        "requires niter >= 2")
