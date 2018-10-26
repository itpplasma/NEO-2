!> Control parameters from input file
module neo_control
  use nrtype, only : dp, input_array_size

  character(20)                      :: in_file
  integer                            :: theta_n
  integer                            :: phi_n
  integer                            :: s_ind_in
  integer                            :: write_progress
  integer                            :: write_output_files
  integer                            :: calc_fourier
  integer                            :: spline_test
  integer                            :: max_m_mode, max_n_mode
  integer                            :: lab_swi, inp_swi, ref_swi, eout_swi
  integer                            :: chk_swi
  integer                            :: fluxs_interp
  integer                            :: s_num
  real(kind=dp)                      :: s_start, s_end
  integer                            :: g11_swi
  integer                            :: eval_mode
  integer                            :: no_fluxs, no_fluxs_s
  integer, dimension(input_array_size) :: fluxs_arr
end module neo_control
