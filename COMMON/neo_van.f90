module neo_van
  use nrtype, only : dp, input_array_size

  real(kind=dp)                        :: v_phi0, v_theta0
  real(kind=dp)                        :: bmin_tol
  integer                              :: v_nper, v_steps
  integer                              :: v_num_mm
  integer                              :: no_minima
  integer, dimension(input_array_size) :: li_minima
  integer                              :: no_gamma
  integer                              :: tau_num
  integer                              :: tau_max_iter
  real(kind=dp)                        :: lambda_fac
  real(kind=dp)                        :: temp_e
  real(kind=dp)                        :: gamma_eps
  real(kind=dp)                        :: phi_eps
end module neo_van
