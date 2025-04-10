!> Control parameters from input file
MODULE neo_control
  use nrtype, only : dp
  CHARACTER(20)                      :: in_file
  INTEGER                            :: theta_n
  INTEGER                            :: phi_n
  INTEGER                            :: s_ind_in
  INTEGER                            :: write_progress
  INTEGER                            :: write_output_files
  INTEGER                            :: calc_fourier
  INTEGER                            :: spline_test
  INTEGER                            :: max_m_mode, max_n_mode
  INTEGER                            :: lab_swi, inp_swi, ref_swi, eout_swi
  INTEGER                            :: chk_swi
  INTEGER                            :: fluxs_interp
  INTEGER                            :: s_num
  REAL(kind=dp)                      :: s_start, s_end
  INTEGER                            :: g11_swi
  INTEGER                            :: eval_mode
  INTEGER                            :: no_fluxs, no_fluxs_s
  INTEGER, DIMENSION(:), ALLOCATABLE :: fluxs_arr

  !> Controls setting of rt0 and bmref in neo_sub.
  logical :: set_rt0_from_rmnc_for_zero_mode = .true.
END MODULE neo_control
